use super::super::SimpleReactorBVP::{
    DensityModel, FastElemReact, NrbvpHandoffConfig, R_G, SimpleReactorTask,
};
use super::super::reactor_BVP_utils::ScalingConfig;
use super::super::solver_backend::{ReactorBvpSolverConfig, ReactorBvpSymbolicBackend};
use RustedSciThe::numerical::BVP_CorrectionNewton::{
    CorrectionNewtonMatrixLayout, CorrectionNewtonOptions,
};
use RustedSciThe::numerical::BVP_Damp::NR_Damp_solver_damped::{AdaptiveGridConfig, SolverParams};
use RustedSciThe::numerical::BVP_Damp::NR_Damp_solver_frozen::FrozenSolverOptions;
use RustedSciThe::numerical::BVP_Damp::grid_api::GridRefinementMethod;
use RustedSciThe::numerical::BVP_api::{
    BvpApiTelemetryMode, BvpSolverSettings, BvpSolverTelemetrySnapshot,
};
use RustedSciThe::numerical::logging::SolverLogLevel;
use nalgebra::{DMatrix, DVector};
use std::collections::HashMap;

pub(super) const T0_K: f64 = 298.15;
pub(super) const T_SCALE_K: f64 = 100.0;
pub(super) const LENGTH_M: f64 = 0.01;
pub(super) const CP_J_PER_KG_K: f64 = 1000.0;
pub(super) const MASS_FLUX: f64 = 0.01;
pub(super) const DENSITY: f64 = 1.0;
pub(super) const LAMBDA: f64 = 0.1;
pub(super) const MOLAR_MASS_KG: f64 = 0.060052;
pub(super) const MOLAR_MASS_G: f64 = MOLAR_MASS_KG * 1000.0;

#[derive(Clone, Copy, Debug)]
pub(super) struct Reaction {
    pub(super) from: usize,
    pub(super) to: usize,
    pub(super) rate: f64,
}

#[derive(Clone, Debug)]
pub(super) struct Network {
    pub(super) rates: Vec<Reaction>,
    pub(super) inlet: Vec<f64>,
    pub(super) diffusion: Vec<f64>,
    pub(super) kinetic_scale: f64,
}

impl Network {
    pub(super) fn small() -> Self {
        Self {
            rates: vec![
                Reaction {
                    from: 0,
                    to: 1,
                    rate: 0.25,
                },
                Reaction {
                    from: 1,
                    to: 2,
                    rate: 0.60,
                },
                Reaction {
                    from: 2,
                    to: 3,
                    rate: 1.00,
                },
                Reaction {
                    from: 0,
                    to: 2,
                    rate: 0.35,
                },
                Reaction {
                    from: 1,
                    to: 3,
                    rate: 0.45,
                },
            ],
            inlet: vec![0.50, 0.25, 0.15, 0.10],
            diffusion: vec![1e-4; 4],
            kinetic_scale: 1.0,
        }
    }

    pub(super) fn large() -> Self {
        let edges = [
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
        Self {
            rates: edges
                .into_iter()
                .map(|(from, to, rate)| Reaction { from, to, rate })
                .collect(),
            inlet: vec![
                0.12, 0.08, 0.07, 0.09, 0.10, 0.06, 0.11, 0.05, 0.07, 0.08, 0.09, 0.08,
            ],
            diffusion: vec![1e-4; 12],
            kinetic_scale: 1.0,
        }
    }

    pub(super) fn species_count(&self) -> usize {
        self.inlet.len()
    }

    pub(super) fn source_matrix(&self) -> DMatrix<f64> {
        let n = self.species_count();
        let mut w = DMatrix::zeros(n, n);
        for reaction in &self.rates {
            let k = reaction.rate * self.kinetic_scale;
            w[(reaction.from, reaction.from)] -= DENSITY * k;
            w[(reaction.to, reaction.from)] += DENSITY * k;
        }
        w
    }

    fn transport_matrix(&self) -> DMatrix<f64> {
        let n = self.species_count();
        let mut a = DMatrix::zeros(2 * n, 2 * n);
        let w = self.source_matrix();
        for i in 0..n {
            let k_i = DENSITY * self.diffusion[i];
            let pe_i = MASS_FLUX * LENGTH_M / k_i;
            a[(i, n + i)] = 1.0 / k_i;
            a[(n + i, n + i)] = pe_i;
            for j in 0..n {
                a[(n + i, j)] = -LENGTH_M.powi(2) * w[(i, j)];
            }
        }
        a
    }
}

#[derive(Debug)]
struct ReferenceProfile {
    concentrations: DVector<f64>,
    internal_fluxes: DVector<f64>,
    temperature_k: f64,
    heat_flux_w_m2: f64,
}

pub(super) fn matrix_exp(matrix: &DMatrix<f64>) -> DMatrix<f64> {
    let n = matrix.nrows();
    let mut result = DMatrix::zeros(n, n);
    for col in 0..n {
        let mut basis = DVector::zeros(n);
        basis[col] = 1.0;
        result.set_column(
            col,
            &crate::SimpleReactorsIVP::combustion_validation::matrix_exponential::expm_action(
                matrix, &basis,
            ),
        );
    }
    result
}

fn reference_initial_state(network: &Network) -> DVector<f64> {
    let n = network.species_count();
    let a = network.transport_matrix();
    let e1 = matrix_exp(&a);
    let e_jc = e1.view((n, 0), (n, n)).into_owned();
    let e_jj = e1.view((n, n), (n, n)).into_owned();
    let singular_values = e_jj.clone().svd(false, false).singular_values;
    let sigma_max = singular_values.iter().copied().fold(0.0_f64, f64::max);
    let sigma_min = singular_values
        .iter()
        .copied()
        .fold(f64::INFINITY, f64::min);
    assert!(
        sigma_min > 0.0 && sigma_min.is_finite(),
        "mixed-BC matrix is numerically singular"
    );
    let condition_estimate = sigma_max / sigma_min;
    assert!(
        condition_estimate < 1e10,
        "mixed-BC matrix condition estimate={condition_estimate:.3e}"
    );
    println!("BVP-02 reference species={n} cond2(E_JJ)={condition_estimate:.6e}");

    let inlet = DVector::from_vec(network.inlet.clone());
    let right_flux = DVector::zeros(n);
    let rhs = right_flux - e_jc * &inlet;
    let initial_flux = e_jj
        .clone()
        .lu()
        .solve(&rhs)
        .expect("mixed boundary linear system should be solvable");
    let linear_residual = &e_jj * &initial_flux - &rhs;
    let relative_residual = linear_residual.norm() / (1.0 + rhs.norm());
    assert!(
        relative_residual < 1e-12,
        "mixed-BC linear residual={relative_residual:.3e}"
    );

    let mut initial = DVector::zeros(2 * n);
    initial.rows_mut(0, n).copy_from(&inlet);
    initial.rows_mut(n, n).copy_from(&initial_flux);
    initial
}

fn reference_profile(network: &Network, z: f64, initial: &DVector<f64>) -> ReferenceProfile {
    let n = network.species_count();
    let state = crate::SimpleReactorsIVP::combustion_validation::matrix_exponential::expm_action(
        &(network.transport_matrix() * z),
        initial,
    );
    let pe_q = MASS_FLUX * CP_J_PER_KG_K * LENGTH_M / LAMBDA;
    let q_right = 1000.0;
    let heat_flux_w_m2 = q_right * (pe_q * (z - 1.0)).exp();
    let integral = ((pe_q * (z - 1.0)).exp() - (-pe_q).exp()) / pe_q;
    let temperature_k = T0_K + q_right * LENGTH_M / LAMBDA * integral;
    ReferenceProfile {
        concentrations: state.rows(0, n).into_owned(),
        internal_fluxes: state.rows(n, n).into_owned(),
        temperature_k,
        heat_flux_w_m2,
    }
}

pub(super) fn fixture(network: &Network) -> SimpleReactorTask {
    let names: Vec<String> = (0..network.species_count())
        .map(|i| format!("X{i}"))
        .collect();
    let reactions = network
        .rates
        .iter()
        .map(|reaction| FastElemReact {
            eq: format!("X{}=>X{}", reaction.from, reaction.to),
            A: reaction.rate * network.kinetic_scale,
            n: 0.0,
            E: 0.0,
            Q: 0.0,
        })
        .collect();
    let mut reactor = SimpleReactorTask::new();
    reactor.set_density_model(DensityModel::Constant);
    reactor
        .fast_react_set(reactions)
        .expect("BVP-02 reaction table should be valid");
    reactor.kindata.substances = names.clone();
    reactor.kindata.groups = Some(
        names
            .iter()
            .map(|name| {
                (
                    name.clone(),
                    HashMap::from([
                        ("C".to_string(), 2),
                        ("H".to_string(), 4),
                        ("O".to_string(), 2),
                    ]),
                )
            })
            .collect(),
    );

    let mut boundary = HashMap::from([("T".to_string(), T0_K), ("q".to_string(), 1000.0)]);
    for (i, name) in names.iter().enumerate() {
        boundary.insert(name.clone(), network.inlet[i]);
        boundary.insert(format!("{name}_flux"), 0.0);
    }
    let diffusion = names
        .iter()
        .enumerate()
        .map(|(i, name)| (name.clone(), network.diffusion[i]))
        .collect();
    let pressure = DENSITY * R_G * T0_K / MOLAR_MASS_KG;
    reactor.set_parameters(
        vec![0.0; network.rates.len()],
        pressure,
        T0_K,
        CP_J_PER_KG_K,
        boundary,
        LAMBDA,
        diffusion,
        MASS_FLUX,
        ScalingConfig::new(T0_K, LENGTH_M, T_SCALE_K),
    );
    reactor.M = MOLAR_MASS_KG;
    reactor
}

fn initial_guess(nodes: usize, network: &Network) -> DMatrix<f64> {
    DMatrix::from_fn(2 + 2 * network.species_count(), nodes, |row, _| match row {
        0 => 0.0,
        1 => 0.0,
        _ if row % 2 == 0 => network.inlet[(row - 2) / 2],
        _ => 0.0,
    })
}

#[derive(Clone, Copy, Debug)]
enum Route {
    Damped,
    Frozen,
    Sci,
}

fn solve(network: &Network, route: Route) -> SimpleReactorTask {
    let nodes = match route {
        Route::Damped => 64,
        Route::Frozen => 1025,
        Route::Sci if network.species_count() > 4 => 17,
        Route::Sci => 33,
    };
    solve_with_nodes_and_frontend(network, route, nodes, ReactorBvpSymbolicBackend::ExprLegacy)
}

fn solve_with_nodes(network: &Network, route: Route, nodes: usize) -> SimpleReactorTask {
    solve_with_nodes_and_frontend(network, route, nodes, ReactorBvpSymbolicBackend::ExprLegacy)
}

fn solve_correction_newton(network: &Network, nodes: usize) -> SimpleReactorTask {
    let mut reactor = fixture(network);
    reactor
        .setup_bvp()
        .expect("BVP-02 CorrectionNewton fixture setup");
    reactor.solver.arg_name = "z".into();
    reactor.solver.x_range = (0.0, 1.0);
    reactor.solver = reactor
        .solver
        .with_bvp_telemetry_mode(BvpApiTelemetryMode::Detailed);
    let mut options = CorrectionNewtonOptions {
        matrix_layout: CorrectionNewtonMatrixLayout::Sparse,
        absolute_tolerance: 1e-9,
        default_relative_tolerance: 1e-8,
        max_iterations: 100,
        logging_level: SolverLogLevel::Off,
        ..CorrectionNewtonOptions::default()
    };
    options.relative_tolerances = reactor
        .solver
        .unknowns
        .iter()
        .map(|name| (name.clone(), 1e-8))
        .collect();
    options.bounds = reactor
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
        .collect();
    reactor
        .solver
        .solve_with_hybrid_settings(
            BvpSolverSettings::CorrectionNewton(options),
            initial_guess(nodes, network),
            0.0,
            1.0,
            nodes,
        )
        .expect("BVP-02 CorrectionNewton should solve the analytic network");
    reactor
}

fn solve_with_nodes_and_frontend(
    network: &Network,
    route: Route,
    nodes: usize,
    symbolic_backend: ReactorBvpSymbolicBackend,
) -> SimpleReactorTask {
    solve_with_nodes_frontend_and_grid(
        network,
        route,
        nodes,
        symbolic_backend,
        GridRefinementMethod::DoublePoints,
        2,
    )
}

fn solve_with_nodes_frontend_and_grid(
    network: &Network,
    route: Route,
    nodes: usize,
    symbolic_backend: ReactorBvpSymbolicBackend,
    grid_method: GridRefinementMethod,
    max_refinements: usize,
) -> SimpleReactorTask {
    let mut reactor = fixture(network);
    reactor
        .setup_bvp()
        .expect("BVP-02 fixture should prepare via production setup");
    reactor.solver.arg_name = "z".to_string();
    reactor.solver.x_range = (0.0, 1.0);
    reactor.solver = reactor
        .solver
        .with_bvp_telemetry_mode(BvpApiTelemetryMode::Counters);
    match route {
        Route::Damped => {
            let strategy = SolverParams {
                adaptive: Some(AdaptiveGridConfig {
                    version: 1,
                    max_refinements,
                    grid_method,
                }),
                ..SolverParams::default()
            };
            let config = NrbvpHandoffConfig::new(
                initial_guess(nodes, network),
                0.0,
                1.0,
                nodes,
                "forward".to_string(),
                "Damped".to_string(),
                Some(strategy),
                None,
                "Banded".to_string(),
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
                Some("off".to_string()),
            )
            .with_solver_backend_config(ReactorBvpSolverConfig {
                symbolic_backend,
                ..ReactorBvpSolverConfig::default_lambdify()
            });
            reactor
                .solver
                .solve_with_bvp_api(config)
                .expect("Damped BVP-02 solve");
        }
        Route::Frozen => reactor
            .solver
            .solve_BVPfrozen(
                initial_guess(nodes, network),
                FrozenSolverOptions::default()
                    .with_tolerance(1e-9)
                    .with_max_iterations(500),
            )
            .expect("Frozen BVP-02 solve"),
        Route::Sci => reactor
            .solver
            .solve_BVPsci(
                initial_guess(nodes, network),
                nodes,
                if network.species_count() > 4 {
                    128
                } else {
                    512
                },
                1e-8,
            )
            .expect("BVP_sci BVP-02 solve"),
    }
    reactor
}

#[derive(Debug, Default)]
struct ProfileErrors {
    max_mass_fraction: f64,
    max_physical_flux: f64,
    max_temperature_k: f64,
    max_heat_flux: f64,
    max_mass_invariant_drift: f64,
    max_total_fraction_drift: Option<f64>,
    max_total_physical_flux: Option<f64>,
}

fn compare_published_profile(network: &Network, reactor: &mut SimpleReactorTask) -> ProfileErrors {
    reactor
        .postprocessing()
        .expect("BVP-02 dimensional publication");
    let solution = reactor
        .solver
        .get_solution()
        .expect("published BVP-02 solution");
    let mesh = reactor
        .solver
        .x_mesh
        .as_ref()
        .expect("published BVP-02 mesh");
    assert!(solution.iter().all(|value| value.is_finite()));
    let initial = reference_initial_state(network);
    let mut errors = ProfileErrors::default();
    let mut mass_invariant_initial = None;
    for (node, x) in mesh.iter().enumerate() {
        let reference = reference_profile(network, x / LENGTH_M, &initial);
        errors.max_temperature_k = errors
            .max_temperature_k
            .max((solution[(node, 0)] - reference.temperature_k).abs());
        errors.max_heat_flux = errors
            .max_heat_flux
            .max((solution[(node, 1)] - reference.heat_flux_w_m2).abs());
        for i in 0..network.species_count() {
            let mass_fraction = solution[(node, 2 + 2 * i)];
            assert!(
                mass_fraction >= -1e-8,
                "negative BVP-02 mass fraction at node={node} species=X{i}: {mass_fraction:.6e}"
            );
            errors.max_mass_fraction = errors
                .max_mass_fraction
                .max((mass_fraction - reference.concentrations[i]).abs());
            errors.max_physical_flux = errors
                .max_physical_flux
                .max((solution[(node, 3 + 2 * i)] - reference.internal_fluxes[i] / LENGTH_M).abs());
        }
        let total_fraction: f64 = (0..network.species_count())
            .map(|i| solution[(node, 2 + 2 * i)])
            .sum();
        let total_physical_flux: f64 = (0..network.species_count())
            .map(|i| solution[(node, 3 + 2 * i)])
            .sum();
        let invariant = LENGTH_M * total_physical_flux - MASS_FLUX * LENGTH_M * total_fraction;
        let initial_invariant = *mass_invariant_initial.get_or_insert(invariant);
        errors.max_mass_invariant_drift = errors
            .max_mass_invariant_drift
            .max((invariant - initial_invariant).abs());
        if network
            .diffusion
            .iter()
            .all(|value| *value == network.diffusion[0])
        {
            errors.max_total_fraction_drift = Some(
                errors
                    .max_total_fraction_drift
                    .unwrap_or_default()
                    .max((total_fraction - 1.0).abs()),
            );
            errors.max_total_physical_flux = Some(
                errors
                    .max_total_physical_flux
                    .unwrap_or_default()
                    .max(total_physical_flux.abs()),
            );
        }
    }
    errors
}

fn assert_source_matrix_and_production_rhs(network: &Network) {
    let mut reactor = fixture(network);
    reactor.setup_bvp().expect("BVP-02 RHS fixture setup");
    let n = network.species_count();
    let source = network.source_matrix();
    for col in 0..n {
        let mass_sum: f64 = (0..n).map(|row| source[(row, col)]).sum();
        assert!(
            mass_sum.abs() < 1e-12,
            "mass source invariant column {col}: {mass_sum}"
        );
    }
    let molar_masses = reactor
        .kindata
        .stecheodata
        .vec_of_molmasses
        .as_ref()
        .expect("BVP-02 kinetic preprocessing should publish molar masses");
    for (index, mass) in molar_masses.iter().enumerate() {
        assert!(
            (mass - MOLAR_MASS_G).abs() < 2e-3,
            "X{index} parsed molar mass={mass}, expected {MOLAR_MASS_G} g/mol"
        );
    }
    for reaction in &reactor.kindata.stecheodata.stecheo_matrx {
        let net_stoichiometry: f64 = reaction.iter().sum();
        assert!(
            net_stoichiometry.abs() < 1e-14,
            "mass/element stoichiometry: {reaction:?}"
        );
    }

    let test_fractions: Vec<f64> = (0..n).map(|i| 0.02 + 0.01 * i as f64).collect();
    let mut state = vec![0.0; 2 + 2 * n];
    state[0] = 0.1;
    state[1] = 0.025;
    for i in 0..n {
        state[2 + 2 * i] = test_fractions[i];
        state[3 + 2 * i] = 1e-5 * (i as f64 - 2.0);
    }
    let variables: Vec<&str> = reactor.solver.unknowns.iter().map(String::as_str).collect();
    let expected_source = &source * DVector::from_vec(test_fractions);
    let mut expected_rhs = vec![0.0; state.len()];
    let mut expected_jacobian = DMatrix::zeros(state.len(), state.len());
    let pe_q = MASS_FLUX * CP_J_PER_KG_K * LENGTH_M / LAMBDA;
    expected_rhs[0] = state[1] / LAMBDA;
    expected_rhs[1] = pe_q * state[1];
    expected_jacobian[(0, 1)] = 1.0 / LAMBDA;
    expected_jacobian[(1, 1)] = pe_q;
    for i in 0..n {
        let k_i = DENSITY * network.diffusion[i];
        let pe_i = MASS_FLUX * LENGTH_M / k_i;
        expected_rhs[2 + 2 * i] = state[3 + 2 * i] / k_i;
        expected_rhs[3 + 2 * i] = pe_i * state[3 + 2 * i] - LENGTH_M.powi(2) * expected_source[i];
        expected_jacobian[(2 + 2 * i, 3 + 2 * i)] = 1.0 / k_i;
        expected_jacobian[(3 + 2 * i, 3 + 2 * i)] = pe_i;
        for j in 0..n {
            expected_jacobian[(3 + 2 * i, 2 + 2 * j)] = -LENGTH_M.powi(2) * source[(i, j)];
        }
    }
    for (row, equation) in reactor.solver.eq_system.iter().enumerate() {
        let actual = equation.eval_expression(&variables, &state);
        assert!(
            (actual - expected_rhs[row]).abs() <= 2e-11 * (1.0 + expected_rhs[row].abs()),
            "RHS row={row} actual={actual:.16e} expected={:.16e}",
            expected_rhs[row]
        );
        for (col, variable) in variables.iter().enumerate() {
            let derivative = equation.diff(variable).eval_expression(&variables, &state);
            assert!(
                (derivative - expected_jacobian[(row, col)]).abs()
                    <= 1e-10 * (1.0 + expected_jacobian[(row, col)].abs()),
                "Jacobian row={row} col={col} symbolic={derivative:.8e} analytic={:.8e}",
                expected_jacobian[(row, col)]
            );
            let step = 1e-7;
            let mut shifted = state.clone();
            shifted[col] += step;
            let numerical = (equation.eval_expression(&variables, &shifted) - actual) / step;
            assert!(
                (derivative - numerical).abs() <= 2e-6 * (1.0 + derivative.abs()),
                "Jacobian row={row} col={col} symbolic={derivative:.8e} finite_difference={numerical:.8e}"
            );
        }
    }
}

fn assert_route_matches_reference(network: &Network, route: Route) -> ProfileErrors {
    let mut reactor = solve(network, route);
    let errors = compare_published_profile(network, &mut reactor);
    let mesh = reactor.solver.x_mesh.as_ref().unwrap();
    let total_c_drift = errors
        .max_total_fraction_drift
        .map(|value| format!("{value:.3e}"))
        .unwrap_or_else(|| "n/a".to_string());
    let total_flux = errors
        .max_total_physical_flux
        .map(|value| format!("{value:.3e}"))
        .unwrap_or_else(|| "n/a".to_string());
    println!(
        "BVP-02 route={route:?} species={} reactions={} nodes={} max_dC={:.3e} max_dJ={:.3e} max_dT={:.3e} K max_dq={:.3e} W/m2 mass_invariant={:.3e} total_C_drift={total_c_drift} total_J={total_flux}",
        network.species_count(),
        network.rates.len(),
        mesh.len(),
        errors.max_mass_fraction,
        errors.max_physical_flux,
        errors.max_temperature_k,
        errors.max_heat_flux,
        errors.max_mass_invariant_drift
    );
    let (c_tol, j_tol, t_tol, q_tol) = match route {
        Route::Damped => (3e-4, 1e-5, 8e-2, 8e-1),
        Route::Frozen => (1e-4, 1e-5, 3e-2, 3e-1),
        Route::Sci => (2e-5, 2e-6, 2e-4, 2e-3),
    };
    assert!(
        errors.max_mass_fraction < c_tol,
        "{route:?} mass-fraction error: {errors:?}"
    );
    assert!(
        errors.max_physical_flux < j_tol,
        "{route:?} physical-flux error: {errors:?}"
    );
    assert!(
        errors.max_temperature_k < t_tol,
        "{route:?} temperature error: {errors:?}"
    );
    assert!(
        errors.max_heat_flux < q_tol,
        "{route:?} heat-flux error: {errors:?}"
    );
    assert!(
        errors.max_mass_invariant_drift < 1e-8,
        "{route:?} mass invariant: {errors:?}"
    );
    if network
        .diffusion
        .iter()
        .all(|value| *value == network.diffusion[0])
    {
        assert!(
            errors.max_total_fraction_drift.unwrap() < 1e-3,
            "{route:?} total composition: {errors:?}"
        );
        assert!(
            errors.max_total_physical_flux.unwrap() < 1e-6,
            "{route:?} total flux: {errors:?}"
        );
    }
    errors
}

#[test]
fn bvp02a_damped_expr_and_atom_frontends_match_the_same_reference() {
    let network = Network::small();
    let mut expr = solve_with_nodes_and_frontend(
        &network,
        Route::Damped,
        64,
        ReactorBvpSymbolicBackend::ExprLegacy,
    );
    let mut atom = solve_with_nodes_and_frontend(
        &network,
        Route::Damped,
        64,
        ReactorBvpSymbolicBackend::AtomView,
    );
    let expr_errors = compare_published_profile(&network, &mut expr);
    let atom_errors = compare_published_profile(&network, &mut atom);
    assert!(
        expr_errors.max_mass_fraction < 3e-4,
        "ExprLegacy: {expr_errors:?}"
    );
    assert!(
        atom_errors.max_mass_fraction < 3e-4,
        "AtomView: {atom_errors:?}"
    );
    assert!(
        atom_errors.max_physical_flux < 1e-5,
        "AtomView: {atom_errors:?}"
    );
    assert!(
        atom_errors.max_temperature_k < 8e-2,
        "AtomView: {atom_errors:?}"
    );
    assert!(
        atom_errors.max_heat_flux < 8e-1,
        "AtomView: {atom_errors:?}"
    );
    let expr_solution = expr.solver.get_solution().unwrap();
    let atom_solution = atom.solver.get_solution().unwrap();
    assert_eq!(expr_solution.shape(), atom_solution.shape());
    let max_frontend_delta = expr_solution
        .iter()
        .zip(atom_solution.iter())
        .map(|(left, right)| (left - right).abs())
        .fold(0.0_f64, f64::max);
    println!(
        "BVP-02a Damped frontend parity ExprLegacy/AtomView nodes={} max_profile_delta={max_frontend_delta:.3e}",
        expr_solution.nrows()
    );
    assert!(max_frontend_delta < 2e-8);
}

#[test]
fn bvp02a_damped_and_frozen_mesh_refinement_reduces_reference_error() {
    let network = Network::small();
    for (route, mesh_sizes) in [
        (Route::Damped, vec![32, 64, 128]),
        (Route::Frozen, vec![129, 257, 513, 1025]),
    ] {
        let mut errors = Vec::new();
        for nodes in mesh_sizes {
            let mut reactor = solve_with_nodes(&network, route, nodes);
            let error = compare_published_profile(&network, &mut reactor);
            println!(
                "BVP-02a mesh convergence route={route:?} seed={nodes} published_nodes={} dC={:.3e} dJ={:.3e} dT={:.3e} K dq={:.3e} W/m2",
                reactor.solver.x_mesh.as_ref().unwrap().len(),
                error.max_mass_fraction,
                error.max_physical_flux,
                error.max_temperature_k,
                error.max_heat_flux
            );
            errors.push(error);
        }
        for pair in errors.windows(2) {
            assert!(
                pair[1].max_mass_fraction < pair[0].max_mass_fraction
                    && pair[1].max_physical_flux < pair[0].max_physical_flux
                    && pair[1].max_temperature_k < pair[0].max_temperature_k
                    && pair[1].max_heat_flux < pair[0].max_heat_flux,
                "{route:?} mesh refinement should reduce each physical-field error: {errors:?}"
            );
        }
    }
}

#[test]
fn bvp02a_damped_refinement_parameter_sweep_measures_physical_accuracy() {
    let network = Network::small();
    let seed_nodes = 33;
    let candidates = [
        (
            "TwoPoint-2e-2-C1.4",
            GridRefinementMethod::TwoPoint(2e-2, 2e-2, 1.4),
        ),
        (
            "TwoPoint-5e-2-C2.5",
            GridRefinementMethod::TwoPoint(5e-2, 5e-2, 2.5),
        ),
        ("Pearson-1e-2", GridRefinementMethod::Pearson(1e-2, 1.4)),
        (
            "Grcar-1e-2-5e-2",
            GridRefinementMethod::GrcarSmooke(1e-2, 5e-2, 1.4),
        ),
    ];
    let mut baseline = solve_with_nodes_frontend_and_grid(
        &network,
        Route::Damped,
        seed_nodes,
        ReactorBvpSymbolicBackend::ExprLegacy,
        GridRefinementMethod::DoublePoints,
        0,
    );
    let baseline_errors = compare_published_profile(&network, &mut baseline);
    let mut measured = Vec::new();

    for (name, grid_method) in candidates {
        let mut reactor = solve_with_nodes_frontend_and_grid(
            &network,
            Route::Damped,
            seed_nodes,
            ReactorBvpSymbolicBackend::ExprLegacy,
            grid_method,
            3,
        );
        let errors = compare_published_profile(&network, &mut reactor);
        let nodes = reactor.solver.x_mesh.as_ref().unwrap().len();
        let BvpSolverTelemetrySnapshot::Damped(telemetry) = reactor
            .solver
            .last_telemetry
            .as_ref()
            .expect("adaptive Damped solve should preserve telemetry")
        else {
            panic!("BVP-02 parameter sweep must use Damped telemetry");
        };
        let refinements = telemetry.counters.grid_refinements;
        assert!(nodes >= seed_nodes);
        assert!(refinements <= 3);
        println!(
            "BVP-02a tuned grid={name} nodes={nodes} refinements={refinements} dY={:.3e} dJ={:.3e} dT={:.3e}K dq={:.3e}W/m2",
            errors.max_mass_fraction,
            errors.max_physical_flux,
            errors.max_temperature_k,
            errors.max_heat_flux,
        );
        measured.push((name, nodes, refinements, errors));
    }
    println!(
        "BVP-02a fixed seed nodes={seed_nodes} dY={:.3e} dJ={:.3e} dT={:.3e}K dq={:.3e}W/m2",
        baseline_errors.max_mass_fraction,
        baseline_errors.max_physical_flux,
        baseline_errors.max_temperature_k,
        baseline_errors.max_heat_flux,
    );
    let pearson = measured
        .iter()
        .find(|row| row.0 == "Pearson-1e-2")
        .expect("the Pearson candidate must be measured");
    assert_eq!(pearson.1, 166);
    assert_eq!(pearson.2, 3);
    assert!(
        pearson.3.max_mass_fraction < baseline_errors.max_mass_fraction
            && pearson.3.max_physical_flux < baseline_errors.max_physical_flux
            && pearson.3.max_temperature_k < baseline_errors.max_temperature_k
            && pearson.3.max_heat_flux < baseline_errors.max_heat_flux,
        "the measured Pearson settings must improve every physical field: baseline={baseline_errors:?}, candidate={pearson:?}"
    );
    assert!(
        measured.iter().any(|(_, _, _, errors)| {
            errors.max_mass_fraction < baseline_errors.max_mass_fraction
                || errors.max_physical_flux < baseline_errors.max_physical_flux
                || errors.max_temperature_k < baseline_errors.max_temperature_k
                || errors.max_heat_flux < baseline_errors.max_heat_flux
        }),
        "no parameterized adaptive strategy improved any physical field: baseline={baseline_errors:?}, candidates={measured:?}"
    );
}

#[test]
fn bvp02a_correction_newton_matches_matrix_exponential_and_retains_telemetry() {
    let network = Network::small();
    let mut coarse = solve_correction_newton(&network, 129);
    let coarse_errors = compare_published_profile(&network, &mut coarse);
    let mut fine = solve_correction_newton(&network, 257);
    let fine_errors = compare_published_profile(&network, &mut fine);
    println!(
        "BVP-02a CorrectionNewton fixed mesh coarse=129 dY={:.3e} dJ={:.3e} dT={:.3e}K dq={:.3e}W/m2 fine=257 dY={:.3e} dJ={:.3e} dT={:.3e}K dq={:.3e}W/m2",
        coarse_errors.max_mass_fraction,
        coarse_errors.max_physical_flux,
        coarse_errors.max_temperature_k,
        coarse_errors.max_heat_flux,
        fine_errors.max_mass_fraction,
        fine_errors.max_physical_flux,
        fine_errors.max_temperature_k,
        fine_errors.max_heat_flux,
    );
    assert!(
        fine_errors.max_mass_fraction < coarse_errors.max_mass_fraction
            && fine_errors.max_physical_flux < coarse_errors.max_physical_flux
            && fine_errors.max_temperature_k < coarse_errors.max_temperature_k
            && fine_errors.max_heat_flux < coarse_errors.max_heat_flux,
        "CorrectionNewton fixed-mesh refinement should improve all oracle fields: coarse={coarse_errors:?}, fine={fine_errors:?}"
    );
    let telemetry = fine
        .solver
        .last_telemetry
        .as_ref()
        .expect("KiThe should retain CorrectionNewton telemetry");
    let BvpSolverTelemetrySnapshot::CorrectionNewton(Some(telemetry)) = telemetry else {
        panic!("CorrectionNewton telemetry should retain its typed snapshot");
    };
    assert!(telemetry.residual_evaluations > 0);
    assert!(telemetry.jacobian_evaluations > 0);
    assert!(telemetry.linear_solves > 0);
    assert!(telemetry.residual_callback_time > std::time::Duration::ZERO);
    assert!(telemetry.jacobian_callback_time > std::time::Duration::ZERO);
    assert!(telemetry.linear_solve_time > std::time::Duration::ZERO);
}

#[test]
fn bvp02a_small_network_checks_independent_source_rhs_and_solver_routes() {
    let network = Network::small();
    assert!((network.inlet.iter().sum::<f64>() - 1.0).abs() < 1e-14);
    let expected_w = DMatrix::from_row_slice(
        4,
        4,
        &[
            -0.60, 0.0, 0.0, 0.0, 0.25, -1.05, 0.0, 0.0, 0.35, 0.60, -1.0, 0.0, 0.0, 0.45, 1.0, 0.0,
        ],
    );
    assert!((&network.source_matrix() - expected_w).norm() < 1e-14);
    assert_source_matrix_and_production_rhs(&network);
    for route in [Route::Damped, Route::Frozen, Route::Sci] {
        assert_route_matches_reference(&network, route);
    }
}

#[test]
fn bvp02b_large_network_preserves_mass_and_matches_matrix_exponential() {
    let network = Network::large();
    assert_eq!(network.species_count(), 12);
    assert_eq!(network.rates.len(), 20);
    assert!((network.inlet.iter().sum::<f64>() - 1.0).abs() < 1e-14);
    assert_source_matrix_and_production_rhs(&network);
    for route in [Route::Damped, Route::Frozen, Route::Sci] {
        assert_route_matches_reference(&network, route);
    }
}

#[test]
fn bvp02c_unequal_diffusion_uses_species_keyed_transport_without_normalization() {
    let mut network = Network::large();
    network.diffusion = (0..network.species_count())
        .map(|i| 1e-4 * (1.0 + i as f64 / 10.0))
        .collect();
    assert!(network.diffusion.windows(2).all(|pair| pair[0] < pair[1]));
    assert_source_matrix_and_production_rhs(&network);
    assert_route_matches_reference(&network, Route::Damped);
}

#[test]
fn bvp02d_zero_kinetics_recovers_transport_only_and_varies_continuously() {
    let base = Network::large();
    let mut zero = base.clone();
    zero.kinetic_scale = 0.0;
    assert_eq!(zero.source_matrix().norm(), 0.0);
    assert_route_matches_reference(&zero, Route::Damped);
    let mut unit = base.clone();
    unit.kinetic_scale = 1.0;
    let mut previous_reference = Some(reference_initial_state(&zero));
    for alpha in [0.01, 0.1, 1.0] {
        let mut scaled = base.clone();
        scaled.kinetic_scale = alpha;
        assert_route_matches_reference(&scaled, Route::Damped);
        let initial = reference_initial_state(&scaled);
        let reference_delta = (&initial - previous_reference.as_ref().unwrap()).norm();
        assert!(reference_delta.is_finite() && reference_delta < 2e-5);
        println!("BVP-02d alpha={alpha} reference_initial_step_delta={reference_delta:.6e}");
        previous_reference = Some(initial);
    }
    let zero_profile = reference_initial_state(&zero);
    let unit_profile = reference_initial_state(&unit);
    assert!((&unit_profile - &zero_profile).norm() > 0.0);
}
