use super::bvp_02::{
    CP_J_PER_KG_K, DENSITY, LAMBDA, LENGTH_M, MASS_FLUX, MOLAR_MASS_KG, T_SCALE_K, T0_K, matrix_exp,
};
use crate::SimpleReactorsBVP::SimpleReactorBVP::{
    DensityModel, FastElemReact, NrbvpHandoffConfig, R_G, SimpleReactorTask,
};
use crate::SimpleReactorsBVP::reactor_BVP_utils::ScalingConfig;
use crate::SimpleReactorsBVP::solver_backend::{ReactorBvpSolverConfig, ReactorBvpSymbolicBackend};
use RustedSciThe::numerical::BVP_Damp::NR_Damp_solver_damped::{AdaptiveGridConfig, SolverParams};
use RustedSciThe::numerical::BVP_Damp::NR_Damp_solver_frozen::FrozenSolverOptions;
use RustedSciThe::numerical::BVP_Damp::grid_api::GridRefinementMethod;
use RustedSciThe::numerical::BVP_api::BvpApiTelemetryMode;
use nalgebra::{DMatrix, DVector};
use std::collections::HashMap;

const R_ARRHENIUS: f64 = 8.314;
const T_REF_K: f64 = 330.0;

#[derive(Clone, Debug)]
struct Reaction {
    equation: String,
    reactant: usize,
    stoich: Vec<f64>,
    k_ref: f64,
    ea_j_per_mol: f64,
}

#[derive(Clone, Debug)]
struct Network {
    names: Vec<String>,
    formulas: Vec<HashMap<String, usize>>,
    molar_masses_kg: Vec<f64>,
    inlet: Vec<f64>,
    diffusion: Vec<f64>,
    reactions: Vec<Reaction>,
}

impl Network {
    fn small() -> Self {
        let edges = [
            (0, 1, 0.25, 12.0),
            (1, 2, 0.60, 18.0),
            (2, 3, 1.00, 24.0),
            (0, 2, 0.35, 30.0),
            (1, 3, 0.45, 36.0),
        ];
        let mut reactions = Vec::with_capacity(edges.len());
        for (from, to, k_ref, ea_kj) in edges {
            let mut stoich = vec![0.0; 4];
            stoich[from] = -1.0;
            stoich[to] = 1.0;
            reactions.push(Reaction {
                equation: format!("X{from}=>X{to}"),
                reactant: from,
                stoich,
                k_ref,
                ea_j_per_mol: ea_kj * 1000.0,
            });
        }
        Self {
            names: (0..4).map(|i| format!("X{i}")).collect(),
            formulas: vec![
                HashMap::from([
                    ("C".to_string(), 2),
                    ("H".to_string(), 4),
                    ("O".to_string(), 2)
                ]);
                4
            ],
            molar_masses_kg: vec![MOLAR_MASS_KG; 4],
            inlet: vec![0.50, 0.25, 0.15, 0.10],
            diffusion: vec![1e-4; 4],
            reactions,
        }
    }

    fn ozone() -> Self {
        let mut first = vec![0.0; 4];
        first[0] = -1.0;
        first[1] = 1.0;
        first[2] = 1.0;
        let mut second = vec![0.0; 4];
        second[1] = -1.0;
        second[2] = 2.0;
        Self {
            names: vec!["O3".into(), "O2".into(), "O".into(), "N2".into()],
            formulas: vec![
                HashMap::from([("O".to_string(), 3)]),
                HashMap::from([("O".to_string(), 2)]),
                HashMap::from([("O".to_string(), 1)]),
                HashMap::from([("N".to_string(), 2)]),
            ],
            molar_masses_kg: vec![0.047997, 0.031998, 0.015999, 0.028014],
            inlet: vec![0.35, 0.25, 0.10, 0.30],
            diffusion: vec![1e-4; 4],
            reactions: vec![
                Reaction {
                    equation: "O3=>O2+O".into(),
                    reactant: 0,
                    stoich: first,
                    k_ref: 0.35,
                    ea_j_per_mol: 24_000.0,
                },
                Reaction {
                    equation: "O2=>2O".into(),
                    reactant: 1,
                    stoich: second,
                    k_ref: 0.75,
                    ea_j_per_mol: 32_000.0,
                },
            ],
        }
    }

    fn species_count(&self) -> usize {
        self.names.len()
    }

    fn arrhenius_a(reaction: &Reaction) -> f64 {
        reaction.k_ref * (reaction.ea_j_per_mol / (R_ARRHENIUS * T_REF_K)).exp()
    }

    fn k(&self, reaction_index: usize, temperature_k: f64) -> f64 {
        let reaction = &self.reactions[reaction_index];
        Self::arrhenius_a(reaction) * (-reaction.ea_j_per_mol / (R_ARRHENIUS * temperature_k)).exp()
    }

    fn source_matrix(&self, temperature_k: f64, density: f64) -> DMatrix<f64> {
        let n = self.species_count();
        let mut source = DMatrix::zeros(n, n);
        for reaction in &self.reactions {
            let rate_factor = density * Self::arrhenius_a(reaction)
                / self.molar_masses_kg[reaction.reactant]
                * (-reaction.ea_j_per_mol / (R_ARRHENIUS * temperature_k)).exp();
            for species in 0..n {
                source[(species, reaction.reactant)] +=
                    self.molar_masses_kg[species] * reaction.stoich[species] * rate_factor;
            }
        }
        source
    }

    fn fixture(&self, q_right_w_m2: f64) -> SimpleReactorTask {
        let reactions = self
            .reactions
            .iter()
            .map(|reaction| FastElemReact {
                eq: reaction.equation.clone(),
                A: Self::arrhenius_a(reaction),
                n: 0.0,
                E: reaction.ea_j_per_mol,
                Q: 0.0,
            })
            .collect();
        let mut reactor = SimpleReactorTask::new();
        reactor.set_density_model(DensityModel::Constant);
        reactor.fast_react_set(reactions).expect("BVP-04 reactions");
        reactor.kindata.substances = self.names.clone();
        reactor.kindata.groups = Some(
            self.names
                .iter()
                .cloned()
                .zip(self.formulas.iter().cloned())
                .collect(),
        );
        let mean_molar_mass = 1.0
            / self
                .inlet
                .iter()
                .zip(&self.molar_masses_kg)
                .map(|(fraction, mass)| fraction / mass)
                .sum::<f64>();
        let pressure = DENSITY * R_G * T0_K / mean_molar_mass;
        let mut boundary =
            HashMap::from([("T".to_string(), T0_K), ("q".to_string(), q_right_w_m2)]);
        for (name, fraction) in self.names.iter().zip(&self.inlet) {
            boundary.insert(name.clone(), *fraction);
            boundary.insert(format!("{name}_flux"), 0.0);
        }
        let diffusion = self
            .names
            .iter()
            .cloned()
            .zip(self.diffusion.iter().copied())
            .collect();
        reactor.set_parameters(
            vec![0.0; self.reactions.len()],
            pressure,
            T0_K,
            CP_J_PER_KG_K,
            boundary,
            LAMBDA,
            diffusion,
            MASS_FLUX,
            ScalingConfig::new(T0_K, LENGTH_M, T_SCALE_K),
        );
        reactor.M = mean_molar_mass;
        reactor
    }
}

fn expected_operator(network: &Network, temperature_k: f64, density: f64) -> DMatrix<f64> {
    let n = network.species_count();
    let mut operator = DMatrix::zeros(2 * n, 2 * n);
    let source = network.source_matrix(temperature_k, density);
    for i in 0..n {
        let rho_d = density * network.diffusion[i];
        operator[(i, n + i)] = 1.0 / rho_d;
        operator[(n + i, n + i)] = MASS_FLUX * LENGTH_M / rho_d;
        for j in 0..n {
            operator[(n + i, j)] = -LENGTH_M.powi(2) * source[(i, j)];
        }
    }
    operator
}

fn mixed_boundary_reference(network: &Network, operator: &DMatrix<f64>) -> (DVector<f64>, f64) {
    let n = network.species_count();
    let transition = matrix_exp(operator);
    let jj = transition.view((n, n), (n, n)).into_owned();
    let jc = transition.view((n, 0), (n, n)).into_owned();
    let inlet = DVector::from_vec(network.inlet.clone());
    let rhs = -jc * inlet.clone();
    let initial_flux = jj
        .clone()
        .lu()
        .solve(&rhs)
        .expect("BVP-04 mixed boundary system");
    let relative_residual = (&jj * &initial_flux - rhs).norm() / (1.0 + initial_flux.norm());
    let mut initial = DVector::zeros(2 * n);
    initial.rows_mut(0, n).copy_from(&inlet);
    initial.rows_mut(n, n).copy_from(&initial_flux);
    (initial, relative_residual)
}

#[test]
fn bvp04_arrhenius_preprocessing_uses_physical_scaled_temperature_and_j_per_mol() {
    let network = Network::small();
    let mut reactor = network.fixture(0.0);
    reactor.setup_bvp().expect("BVP-04 symbolic setup");
    let rate_constants = reactor
        .kindata
        .K_sym_vec
        .as_ref()
        .expect("symbolic Arrhenius constants");
    for temperature_k in [T0_K, T_REF_K, 361.36205588] {
        let theta = (temperature_k - T0_K) / T_SCALE_K;
        let state = [theta];
        let temperature_state_names = ["Teta"];
        for (index, reaction) in network.reactions.iter().enumerate() {
            let actual = rate_constants[index].eval_expression(&temperature_state_names, &state);
            let expected = network.k(index, temperature_k);
            assert!(
                (actual - expected).abs() <= 2e-12 * expected.abs().max(1.0),
                "k[{index}]({temperature_k} K)={actual:.16e}, expected {expected:.16e}"
            );
            assert!(reaction.ea_j_per_mol > 0.0);
            if (temperature_k - T_REF_K).abs() < f64::EPSILON {
                assert!((actual - reaction.k_ref).abs() < 1e-12);
            }
        }
    }
    assert!(
        network
            .reactions
            .iter()
            .enumerate()
            .all(|(i, _)| network.k(i, 361.36205588) > network.k(i, T_REF_K)
                && network.k(i, T_REF_K) > network.k(i, T0_K))
    );
}

#[test]
fn bvp04_variable_temperature_source_and_full_symbolic_jacobian_match_analytic_blocks() {
    for network in [Network::small(), Network::ozone()] {
        let mut reactor = network.fixture(1_000.0);
        reactor.setup_bvp().expect("BVP-04 derivative setup");
        let n = network.species_count();
        let temperature_k = 330.0;
        let theta = (temperature_k - T0_K) / T_SCALE_K;
        let concentrations: Vec<f64> = (0..n).map(|i| 0.12 + 0.03 * i as f64).collect();
        let mut state = vec![0.0; 2 + 2 * n];
        state[0] = theta;
        state[1] = 0.021;
        for i in 0..n {
            state[2 + 2 * i] = concentrations[i];
            state[3 + 2 * i] = -0.004 + 0.002 * i as f64;
        }
        let density = reactor.ideal_gas_density().unwrap();
        let source = network.source_matrix(temperature_k, density);
        let mut source_derivative: DMatrix<f64> = DMatrix::zeros(n, n);
        for (reaction_index, reaction) in network.reactions.iter().enumerate() {
            let dk_d_t = network.k(reaction_index, temperature_k) * reaction.ea_j_per_mol
                / (R_ARRHENIUS * temperature_k.powi(2));
            let rate_derivative = density * dk_d_t / network.molar_masses_kg[reaction.reactant];
            for i in 0..n {
                source_derivative[(i, reaction.reactant)] +=
                    network.molar_masses_kg[i] * reaction.stoich[i] * rate_derivative;
            }
        }
        let pe_q = MASS_FLUX * CP_J_PER_KG_K * LENGTH_M / LAMBDA;
        let variables: Vec<&str> = reactor.solver.unknowns.iter().map(String::as_str).collect();
        for (row, equation) in reactor.solver.eq_system.iter().enumerate() {
            let actual_rhs = equation.eval_expression(&variables, &state);
            let expected_rhs = if row == 0 {
                state[1] / LAMBDA
            } else if row == 1 {
                pe_q * state[1]
            } else {
                let species = (row - 2) / 2;
                if row % 2 == 0 {
                    state[3 + 2 * species] / (density * network.diffusion[species])
                } else {
                    MASS_FLUX * LENGTH_M / (density * network.diffusion[species])
                        * state[3 + 2 * species]
                        - LENGTH_M.powi(2)
                            * (0..n)
                                .map(|j| source[(species, j)] * concentrations[j])
                                .sum::<f64>()
                }
            };
            assert!(
                (actual_rhs - expected_rhs).abs() < 2e-10 * (1.0 + expected_rhs.abs()),
                "{} RHS[{row}]={actual_rhs:e}, expected {expected_rhs:e}",
                network.names.join(",")
            );
            for (column, variable) in variables.iter().enumerate() {
                let actual = equation.diff(variable).eval_expression(&variables, &state);
                let expected = if row == 0 {
                    if column == 1 { 1.0 / LAMBDA } else { 0.0 }
                } else if row == 1 {
                    if column == 1 { pe_q } else { 0.0 }
                } else {
                    let species = (row - 2) / 2;
                    if row % 2 == 0 {
                        if column == 3 + 2 * species {
                            1.0 / (density * network.diffusion[species])
                        } else {
                            0.0
                        }
                    } else if column == 0 {
                        -LENGTH_M.powi(2)
                            * T_SCALE_K
                            * (0..n)
                                .map(|j| source_derivative[(species, j)] * concentrations[j])
                                .sum::<f64>()
                    } else if column >= 2 && column % 2 == 0 {
                        -LENGTH_M.powi(2) * source[(species, (column - 2) / 2)]
                    } else if column == 3 + 2 * species {
                        MASS_FLUX * LENGTH_M / (density * network.diffusion[species])
                    } else {
                        0.0
                    }
                };
                assert!(
                    (actual - expected).abs() <= 2e-9 * (1.0 + expected.abs()),
                    "{} Jacobian[{row},{column}]={actual:e}, expected {expected:e}",
                    network.names.join(",")
                );
                if column == 0 || column == 2 || column == 3 + 2 * (n / 2) {
                    let step = 1e-7;
                    let mut shifted = state.clone();
                    shifted[column] += step;
                    let finite_difference =
                        (equation.eval_expression(&variables, &shifted) - actual_rhs) / step;
                    assert!(
                        (actual - finite_difference).abs() < 5e-5 * (1.0 + actual.abs()),
                        "{} finite-difference Jacobian[{row},{column}]={finite_difference:e}, symbolic={actual:e}",
                        network.names.join(",")
                    );
                }
            }
        }
    }
}

#[derive(Clone, Copy, Debug)]
enum Route {
    Damped,
    Frozen,
    Sci,
}

#[derive(Debug, Default)]
struct ProfileError {
    c: f64,
    j: f64,
    t: f64,
    q: f64,
    mass_invariant: f64,
    element_invariant: f64,
}

fn parse_reference(csv: &str, species_count: usize) -> Vec<Vec<f64>> {
    let rows: Vec<Vec<f64>> = csv
        .lines()
        .skip(1)
        .map(|line| {
            line.split(',')
                .map(|field| field.parse::<f64>().expect("numeric frozen reference"))
                .collect()
        })
        .collect();
    assert_eq!(rows.len(), 2001);
    assert!(rows.iter().all(|row| row.len() == 2 * species_count + 3));
    rows
}

fn interpolate_reference(rows: &[Vec<f64>], z: f64) -> Vec<f64> {
    let position = z.clamp(0.0, 1.0) * (rows.len() - 1) as f64;
    let left = position.floor() as usize;
    let right = (left + 1).min(rows.len() - 1);
    let fraction = position - left as f64;
    rows[left]
        .iter()
        .zip(&rows[right])
        .map(|(a, b)| a + fraction * (b - a))
        .collect()
}

fn prescribed_temperature_k(z: f64) -> f64 {
    let pe_q = MASS_FLUX * CP_J_PER_KG_K * LENGTH_M / LAMBDA;
    T0_K + (1_000.0 * LENGTH_M / LAMBDA) * ((pe_q * (z - 1.0)).exp() - (-pe_q).exp()) / pe_q
}

fn prescribed_heat_flux_w_m2(z: f64) -> f64 {
    let pe_q = MASS_FLUX * CP_J_PER_KG_K * LENGTH_M / LAMBDA;
    1_000.0 * (pe_q * (z - 1.0)).exp()
}

fn compare_profile(
    network: &Network,
    reactor: &SimpleReactorTask,
    reference_csv: &str,
) -> ProfileError {
    let n = network.species_count();
    let reference = parse_reference(reference_csv, n);
    let solution = reactor
        .solver
        .get_solution()
        .expect("published BVP-04 result");
    let mesh = reactor
        .solver
        .x_mesh
        .as_ref()
        .expect("published BVP-04 mesh");
    let mut error = ProfileError::default();
    let mut invariant_reference = None;
    let mut invariant_scale = 0.0_f64;
    let elements: std::collections::BTreeSet<String> = network
        .formulas
        .iter()
        .flat_map(|formula| formula.keys().cloned())
        .collect();
    let mut element_baselines = HashMap::<String, f64>::new();
    let mut element_scales = HashMap::<String, f64>::new();
    for (node, x) in mesh.iter().copied().enumerate() {
        let z = x / LENGTH_M;
        let expected = interpolate_reference(&reference, z);
        let mut c = 0.0;
        let mut j_internal = 0.0;
        for i in 0..n {
            let actual_c = solution[(node, 2 + 2 * i)];
            let actual_j = solution[(node, 3 + 2 * i)];
            assert!(actual_c.is_finite() && actual_c >= -1e-8);
            assert!(actual_j.is_finite());
            error.c = error.c.max((actual_c - expected[1 + i]).abs());
            error.j = error
                .j
                .max((actual_j - expected[1 + n + i] / LENGTH_M).abs());
            c += actual_c;
            j_internal += actual_j * LENGTH_M;
        }
        error.t = error
            .t
            .max((solution[(node, 0)] - prescribed_temperature_k(z)).abs());
        error.q = error
            .q
            .max((solution[(node, 1)] - prescribed_heat_flux_w_m2(z)).abs());
        let invariant = j_internal - MASS_FLUX * LENGTH_M * c;
        let baseline = *invariant_reference.get_or_insert(invariant);
        error.mass_invariant = error.mass_invariant.max((invariant - baseline).abs());
        invariant_scale = invariant_scale.max(j_internal.abs() + MASS_FLUX * LENGTH_M * c.abs());
        for element in &elements {
            let mut value = 0.0;
            let mut scale = 0.0;
            for i in 0..n {
                let count = network.formulas[i].get(element).copied().unwrap_or(0) as f64;
                let coefficient = count / network.molar_masses_kg[i];
                let term = coefficient
                    * (solution[(node, 3 + 2 * i)] * LENGTH_M
                        - MASS_FLUX * LENGTH_M * solution[(node, 2 + 2 * i)]);
                value += term;
                scale += term.abs();
            }
            let baseline = *element_baselines.entry(element.clone()).or_insert(value);
            error.element_invariant = error.element_invariant.max((value - baseline).abs());
            element_scales
                .entry(element.clone())
                .and_modify(|current| *current = current.max(scale))
                .or_insert(scale);
        }
    }
    error.mass_invariant /= invariant_scale.max(1e-12);
    let element_scale = element_scales.values().copied().fold(0.0, f64::max);
    error.element_invariant /= element_scale.max(1e-12);
    error
}

fn solve_route(
    network: &Network,
    route: Route,
    q_right_w_m2: f64,
    sci_initial_nodes: usize,
) -> SimpleReactorTask {
    solve_route_with_frontend(
        network,
        route,
        q_right_w_m2,
        sci_initial_nodes,
        ReactorBvpSymbolicBackend::ExprLegacy,
    )
}

fn solve_route_with_frontend(
    network: &Network,
    route: Route,
    q_right_w_m2: f64,
    sci_initial_nodes: usize,
    symbolic_backend: ReactorBvpSymbolicBackend,
) -> SimpleReactorTask {
    solve_route_with_grid(
        network,
        route,
        q_right_w_m2,
        sci_initial_nodes,
        symbolic_backend,
        None,
    )
}

fn solve_route_with_grid(
    network: &Network,
    route: Route,
    q_right_w_m2: f64,
    initial_nodes: usize,
    symbolic_backend: ReactorBvpSymbolicBackend,
    damped_grid: Option<(GridRefinementMethod, usize)>,
) -> SimpleReactorTask {
    let nodes = match route {
        Route::Damped if damped_grid.is_some() => initial_nodes,
        Route::Damped => 256,
        Route::Frozen => 2049,
        Route::Sci => initial_nodes,
    };
    let mut reactor = network.fixture(q_right_w_m2);
    reactor.setup_bvp().expect("BVP-04a setup");
    reactor.solver.arg_name = "z".into();
    reactor.solver.x_range = (0.0, 1.0);
    reactor.solver = reactor
        .solver
        .with_bvp_telemetry_mode(BvpApiTelemetryMode::Counters);
    let n = network.species_count();
    let guess = DMatrix::from_fn(2 + 2 * n, nodes, |row, _| match row {
        0 | 1 => 0.0,
        _ if row % 2 == 0 => network.inlet[(row - 2) / 2],
        _ => 0.0,
    });
    match route {
        Route::Damped => {
            let (grid_method, max_refinements) =
                damped_grid.unwrap_or((GridRefinementMethod::DoublePoints, 2));
            let params = SolverParams {
                adaptive: Some(AdaptiveGridConfig {
                    version: 1,
                    max_refinements,
                    grid_method,
                }),
                ..SolverParams::default()
            };
            let config = NrbvpHandoffConfig::new(
                guess,
                0.0,
                1.0,
                nodes,
                "forward".into(),
                "Damped".into(),
                Some(params),
                None,
                "Banded".into(),
                1e-8,
                Some(
                    reactor
                        .solver
                        .unknowns
                        .iter()
                        .map(|name| (name.clone(), 1e-8))
                        .collect(),
                ),
                200,
                Some(
                    reactor
                        .solver
                        .unknowns
                        .iter()
                        .map(|name| {
                            let bounds = if name == "Teta" {
                                (-5.0, 5.0)
                            } else if name == "q" || name.starts_with('J') {
                                (-10.0, 10.0)
                            } else {
                                (-1.0, 2.0)
                            };
                            (name.clone(), bounds)
                        })
                        .collect(),
                ),
                Some("off".into()),
            )
            .with_solver_backend_config(ReactorBvpSolverConfig {
                symbolic_backend,
                ..ReactorBvpSolverConfig::default_lambdify()
            });
            reactor
                .solver
                .solve_with_bvp_api(config)
                .expect("BVP-04a Damped solve");
        }
        Route::Frozen => reactor
            .solver
            .solve_BVPfrozen(
                guess,
                FrozenSolverOptions::default()
                    .with_tolerance(1e-9)
                    .with_max_iterations(500),
            )
            .expect("BVP-04a Frozen solve"),
        Route::Sci => reactor
            .solver
            .solve_BVPsci(guess, nodes, 512, 1e-8)
            .expect("BVP-04a BVP_sci solve"),
    }
    reactor.postprocessing().expect("BVP-04a publication");
    reactor
}

#[test]
fn bvp04a_isothermal_arrhenius_limit_matches_constant_operator_reference() {
    let network = Network::small();
    let source = network.source_matrix(T0_K, DENSITY);
    let operator = expected_operator(&network, T0_K, DENSITY);
    let (initial, boundary_residual) = mixed_boundary_reference(&network, &operator);
    assert!(boundary_residual < 1e-12);
    let end = &matrix_exp(&operator) * &initial;
    assert!(
        end.rows(network.species_count(), network.species_count())
            .norm()
            < 1e-10
    );
    let reactor = network.fixture(0.0);
    assert!(reactor.thermal_effects.iter().all(|value| *value == 0.0));
    for column in 0..source.ncols() {
        let sum: f64 = (0..source.nrows()).map(|row| source[(row, column)]).sum();
        assert!(sum.abs() < 1e-12);
    }
    assert!(operator.iter().all(|value| value.is_finite()));
    for route in [Route::Sci, Route::Frozen, Route::Damped] {
        let solution = solve_route(&network, route, 0.0, 33);
        let profile = solution.solver.get_solution().expect("BVP-04a solution");
        let mesh = solution.solver.x_mesh.as_ref().expect("BVP-04a mesh");
        let mut max_c = 0.0_f64;
        let mut max_j = 0.0_f64;
        let mut max_t = 0.0_f64;
        let mut max_q = 0.0_f64;
        for (node, x) in mesh.iter().copied().enumerate() {
            let z = x / LENGTH_M;
            let exact =
                crate::SimpleReactorsIVP::combustion_validation::matrix_exponential::expm_action(
                    &(operator.clone() * z),
                    &initial,
                );
            for i in 0..network.species_count() {
                max_c = max_c.max((profile[(node, 2 + 2 * i)] - exact[i]).abs());
                max_j = max_j.max(
                    (profile[(node, 3 + 2 * i)] - exact[network.species_count() + i] / LENGTH_M)
                        .abs(),
                );
            }
            max_t = max_t.max((profile[(node, 0)] - T0_K).abs());
            max_q = max_q.max(profile[(node, 1)].abs());
        }
        println!(
            "BVP-04a isothermal route={route:?} nodes={} dC={max_c:.3e} dJ={max_j:.3e} dT={max_t:.3e} dq={max_q:.3e}",
            mesh.len()
        );
        let (c_tol, j_tol) = match route {
            Route::Damped => (3e-4, 1e-5),
            Route::Frozen => (1e-4, 1e-5),
            Route::Sci => (2e-5, 2e-6),
        };
        assert!(
            max_c < c_tol,
            "{route:?} BVP-04a composition error {max_c:e}"
        );
        assert!(max_j < j_tol, "{route:?} BVP-04a flux error {max_j:e}");
        assert!(max_t < 1e-8 && max_q < 1e-8);
    }
}

#[test]
fn bvp04b_small_nonautonomous_network_matches_frozen_dop853_profiles() {
    const REFERENCE: &str = include_str!("fixtures/bvp04/bvp04b_4_5.csv");
    let network = Network::small();
    assert!((prescribed_temperature_k(1.0) - 361.36205588).abs() < 2e-8);
    for route in [Route::Sci, Route::Frozen, Route::Damped] {
        let initial_nodes = if matches!(route, Route::Sci) { 33 } else { 64 };
        let reactor = solve_route(&network, route, 1_000.0, initial_nodes);
        let solution = reactor.solver.get_solution().unwrap();
        let last = solution.nrows() - 1;
        let expected_outlet = interpolate_reference(&parse_reference(REFERENCE, 4), 1.0);
        println!(
            "BVP-04b outlet route={route:?} actual_C={:?} reference_C={:?}",
            (0..4)
                .map(|i| solution[(last, 2 + 2 * i)])
                .collect::<Vec<_>>(),
            &expected_outlet[1..5]
        );
        let errors = compare_profile(&network, &reactor, REFERENCE);
        println!(
            "BVP-04b route={route:?} nodes={} dC={:.3e} dJ={:.3e} dT={:.3e}K dq={:.3e}W/m2 mass_invariant={:.3e} element_invariant={:.3e}",
            reactor.solver.x_mesh.as_ref().unwrap().len(),
            errors.c,
            errors.j,
            errors.t,
            errors.q,
            errors.mass_invariant,
            errors.element_invariant
        );
        let (c_tol, j_tol, t_tol, q_tol) = match route {
            Route::Damped => (3e-4, 1e-5, 8e-2, 8e-1),
            Route::Frozen => (1e-4, 1e-5, 3e-2, 3e-1),
            Route::Sci => (2e-5, 2e-6, 2e-4, 2e-3),
        };
        assert!(errors.c < c_tol, "{route:?}: {errors:?}");
        assert!(errors.j < j_tol, "{route:?}: {errors:?}");
        assert!(errors.t < t_tol, "{route:?}: {errors:?}");
        assert!(errors.q < q_tol, "{route:?}: {errors:?}");
        assert!(errors.mass_invariant < 2e-5, "{route:?}: {errors:?}");
        assert!(errors.element_invariant < 2e-5, "{route:?}: {errors:?}");
    }
}

#[test]
fn bvp04b_damped_refinement_parameter_sweep_measures_nonautonomous_accuracy() {
    const REFERENCE: &str = include_str!("fixtures/bvp04/bvp04b_4_5.csv");
    let network = Network::small();
    let seed_nodes = 65;
    let baseline = solve_route_with_grid(
        &network,
        Route::Damped,
        1_000.0,
        seed_nodes,
        ReactorBvpSymbolicBackend::ExprLegacy,
        Some((GridRefinementMethod::DoublePoints, 0)),
    );
    let baseline_errors = compare_profile(&network, &baseline, REFERENCE);
    let candidates = [
        ("Pearson-1e-2", GridRefinementMethod::Pearson(1e-2, 1.4)),
        (
            "Grcar-1e-2-5e-2",
            GridRefinementMethod::GrcarSmooke(1e-2, 5e-2, 1.4),
        ),
        (
            "TwoPoint-2e-2-C1.4",
            GridRefinementMethod::TwoPoint(2e-2, 2e-2, 1.4),
        ),
    ];
    let mut measured = Vec::new();
    for (name, grid_method) in candidates {
        let reactor = solve_route_with_grid(
            &network,
            Route::Damped,
            1_000.0,
            seed_nodes,
            ReactorBvpSymbolicBackend::ExprLegacy,
            Some((grid_method, 3)),
        );
        let errors = compare_profile(&network, &reactor, REFERENCE);
        let nodes = reactor.solver.x_mesh.as_ref().unwrap().len();
        println!(
            "BVP-04b tuned grid={name} nodes={nodes} dC={:.3e} dJ={:.3e} dT={:.3e}K dq={:.3e}W/m2 mass={:.3e} element={:.3e}",
            errors.c, errors.j, errors.t, errors.q, errors.mass_invariant, errors.element_invariant,
        );
        measured.push((name, nodes, errors));
    }
    println!(
        "BVP-04b fixed seed nodes={seed_nodes} dC={:.3e} dJ={:.3e} dT={:.3e}K dq={:.3e}W/m2",
        baseline_errors.c, baseline_errors.j, baseline_errors.t, baseline_errors.q,
    );
    let grcar = measured
        .iter()
        .find(|row| row.0 == "Grcar-1e-2-5e-2")
        .expect("the Grcar nonautonomous candidate must be measured");
    assert_eq!(grcar.1, 173);
    assert!(
        grcar.2.c < baseline_errors.c
            && grcar.2.j < baseline_errors.j
            && grcar.2.t < baseline_errors.t
            && grcar.2.q < baseline_errors.q
            && grcar.2.mass_invariant < 1e-12
            && grcar.2.element_invariant < 1e-12,
        "selected Grcar settings must improve all fields and preserve invariants: baseline={baseline_errors:?}, candidate={grcar:?}"
    );
    assert!(
        measured.iter().any(|(_, _, errors)| {
            errors.c < baseline_errors.c
                && errors.j < baseline_errors.j
                && errors.t < baseline_errors.t
                && errors.q < baseline_errors.q
        }),
        "no tuned strategy improved all nonautonomous physical fields: baseline={baseline_errors:?}, candidates={measured:?}"
    );
}

#[test]
fn bvp04b_small_network_symbolic_frontends_match_the_same_profile() {
    let network = Network::small();
    let expr = solve_route_with_frontend(
        &network,
        Route::Damped,
        1_000.0,
        64,
        ReactorBvpSymbolicBackend::ExprLegacy,
    );
    let atom = solve_route_with_frontend(
        &network,
        Route::Damped,
        1_000.0,
        64,
        ReactorBvpSymbolicBackend::AtomView,
    );
    let expr_profile = expr.solver.get_solution().unwrap();
    let atom_profile = atom.solver.get_solution().unwrap();
    assert_eq!(expr_profile.shape(), atom_profile.shape());
    let delta = expr_profile
        .iter()
        .zip(atom_profile.iter())
        .map(|(left, right)| (left - right).abs())
        .fold(0.0_f64, f64::max);
    println!(
        "BVP-04b Damped ExprLegacy/AtomView nodes={} max_profile_delta={delta:.3e}",
        expr_profile.nrows()
    );
    assert!(delta < 2e-8);
}

#[test]
#[ignore = "12/20 BVP_sci reference gate is intentionally expensive in Debug; run explicitly in Release"]
fn bvp04c_large_nonautonomous_network_matches_frozen_dop853_profile() {
    const REFERENCE: &str = include_str!("fixtures/bvp04/bvp04c_12_20.csv");
    let mut network = Network::small();
    let large = super::bvp_02::Network::large();
    network.names = (0..large.species_count())
        .map(|i| format!("X{i}"))
        .collect();
    network.formulas = vec![network.formulas[0].clone(); large.species_count()];
    network.molar_masses_kg = vec![MOLAR_MASS_KG; large.species_count()];
    network.inlet = large.inlet.clone();
    network.diffusion = large.diffusion.clone();
    network.reactions = large
        .rates
        .iter()
        .enumerate()
        .map(|(j, edge)| {
            let mut stoich = vec![0.0; large.species_count()];
            stoich[edge.from] = -1.0;
            stoich[edge.to] = 1.0;
            Reaction {
                equation: format!("X{}=>X{}", edge.from, edge.to),
                reactant: edge.from,
                stoich,
                k_ref: edge.rate,
                ea_j_per_mol: (10.0 + 2.0 * j as f64) * 1000.0,
            }
        })
        .collect();
    let reactor = solve_route(&network, Route::Sci, 1_000.0, 17);
    let errors = compare_profile(&network, &reactor, REFERENCE);
    println!(
        "BVP-04c route=Sci nodes={} dC={:.3e} dJ={:.3e} dT={:.3e}K dq={:.3e}W/m2 mass_invariant={:.3e} element_invariant={:.3e}",
        reactor.solver.x_mesh.as_ref().unwrap().len(),
        errors.c,
        errors.j,
        errors.t,
        errors.q,
        errors.mass_invariant,
        errors.element_invariant
    );
    assert!(errors.c < 2e-5, "BVP_sci: {errors:?}");
    assert!(errors.j < 2e-6, "BVP_sci: {errors:?}");
    assert!(errors.t < 2e-4, "BVP_sci: {errors:?}");
    assert!(errors.q < 2e-3, "BVP_sci: {errors:?}");
    assert!(errors.mass_invariant < 2e-5, "BVP_sci: {errors:?}");
    assert!(errors.element_invariant < 2e-5, "BVP_sci: {errors:?}");
}

#[test]
fn bvp04d_unequal_molar_mass_network_matches_frozen_dop853_profiles() {
    const REFERENCE: &str = include_str!("fixtures/bvp04/bvp04d_ozone_4_2.csv");
    let network = Network::ozone();
    for route in [Route::Damped, Route::Frozen, Route::Sci] {
        let initial_nodes = if matches!(route, Route::Sci) { 33 } else { 64 };
        let reactor = solve_route(&network, route, 1_000.0, initial_nodes);
        let errors = compare_profile(&network, &reactor, REFERENCE);
        println!(
            "BVP-04d route={route:?} nodes={} dC={:.3e} dJ={:.3e} dT={:.3e}K dq={:.3e}W/m2 mass_invariant={:.3e} element_invariant={:.3e}",
            reactor.solver.x_mesh.as_ref().unwrap().len(),
            errors.c,
            errors.j,
            errors.t,
            errors.q,
            errors.mass_invariant,
            errors.element_invariant
        );
        let (c_tol, j_tol, t_tol, q_tol) = match route {
            Route::Damped => (3e-4, 1e-5, 8e-2, 8e-1),
            Route::Frozen => (1e-4, 1e-5, 3e-2, 3e-1),
            Route::Sci => (2e-5, 2e-6, 2e-4, 2e-3),
        };
        assert!(errors.c < c_tol, "{route:?}: {errors:?}");
        assert!(errors.j < j_tol, "{route:?}: {errors:?}");
        assert!(errors.t < t_tol, "{route:?}: {errors:?}");
        assert!(errors.q < q_tol, "{route:?}: {errors:?}");
        assert!(errors.mass_invariant < 2e-5, "{route:?}: {errors:?}");
        assert!(errors.element_invariant < 2e-5, "{route:?}: {errors:?}");
    }
}

#[test]
fn bvp04d_ozone_fixture_is_mass_and_element_balanced_with_distinct_molar_masses() {
    let network = Network::ozone();
    let mut reactor = network.fixture(0.0);
    reactor
        .setup_bvp()
        .expect("BVP-04d parser and preprocessing");
    let actual = reactor
        .kindata
        .stecheodata
        .vec_of_molmasses
        .as_ref()
        .expect("production molar masses");
    for (i, mass) in actual.iter().enumerate() {
        assert!(
            (mass / 1000.0 - network.molar_masses_kg[i]).abs() < 2e-6,
            "{} molar mass={} g/mol, expected {} g/mol",
            network.names[i],
            mass,
            network.molar_masses_kg[i] * 1000.0
        );
    }

    let source = network.source_matrix(330.0, DENSITY);
    for reaction_index in 0..network.reactions.len() {
        let reaction = &network.reactions[reaction_index];
        let mass_sum: f64 = (0..network.species_count())
            .map(|i| network.molar_masses_kg[i] * reaction.stoich[i])
            .sum();
        assert!(
            mass_sum.abs() < 1e-14,
            "mass balance reaction {reaction_index}: {mass_sum:e}"
        );
    }
    for column in 0..source.ncols() {
        let source_sum: f64 = (0..source.nrows()).map(|row| source[(row, column)]).sum();
        assert!(
            source_sum.abs() < 1e-12,
            "mass source column {column}: {source_sum:e}"
        );
    }
    let oxygen_atoms = [3.0, 2.0, 1.0, 0.0];
    let nitrogen_atoms = [0.0, 0.0, 0.0, 2.0];
    for (name, atoms) in [("oxygen", oxygen_atoms), ("nitrogen", nitrogen_atoms)] {
        for (j, reaction) in network.reactions.iter().enumerate() {
            let defect: f64 = atoms
                .iter()
                .zip(&reaction.stoich)
                .map(|(count, nu)| count * nu)
                .sum();
            assert!(
                defect.abs() < 1e-14,
                "{name} balance in reaction {j}: {defect}"
            );
        }
    }
    assert!(
        network
            .molar_masses_kg
            .windows(2)
            .any(|pair| pair[0] != pair[1])
    );
}
