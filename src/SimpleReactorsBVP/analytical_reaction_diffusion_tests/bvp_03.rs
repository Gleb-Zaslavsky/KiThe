use super::bvp_02::{
    CP_J_PER_KG_K, DENSITY, LAMBDA, LENGTH_M, MASS_FLUX, MOLAR_MASS_KG, Network, T_SCALE_K, T0_K,
    fixture as base_fixture, matrix_exp,
};
use crate::SimpleReactorsBVP::SimpleReactorBVP::{NrbvpHandoffConfig, SimpleReactorTask};
use crate::SimpleReactorsBVP::solver_backend::{ReactorBvpSolverConfig, ReactorBvpSymbolicBackend};
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

const H_BASE_J_PER_MOL: [f64; 12] = [
    0.0, -40.0, -25.0, -65.0, -15.0, -80.0, -35.0, -100.0, -50.0, -90.0, -60.0, -115.0,
];
const GAMMA: f64 = 300.0;
const Q_RIGHT_PHYSICAL: f64 = 1000.0;

#[derive(Clone, Copy, Debug)]
enum Route {
    Damped,
    Frozen,
    Sci,
}

#[derive(Debug)]
struct ThermoReference {
    stoich: DMatrix<f64>,
    reaction_rates: DMatrix<f64>,
    reaction_heats: DVector<f64>,
    source: DMatrix<f64>,
    operator: DMatrix<f64>,
    initial: DVector<f64>,
    gamma: f64,
    condition_estimate: f64,
}

#[derive(Debug, Default)]
struct Errors {
    c: f64,
    c_z: f64,
    c_relative: f64,
    j: f64,
    j_z: f64,
    j_relative: f64,
    t: f64,
    t_z: f64,
    t_relative: f64,
    chemical_temperature_relative: f64,
    q: f64,
    q_z: f64,
    q_relative: f64,
    composition_sum_deviation: f64,
    energy_abs: f64,
    energy_rel: f64,
    energy_z: f64,
    heat_integral_abs: f64,
    heat_integral_rel: f64,
}

fn enthalpies(network: &Network, gamma: f64) -> DVector<f64> {
    DVector::from_iterator(
        network.species_count(),
        H_BASE_J_PER_MOL[..network.species_count()]
            .iter()
            .map(|h| gamma * h),
    )
}

fn stoichiometry(network: &Network) -> DMatrix<f64> {
    let n = network.species_count();
    let mut s = DMatrix::zeros(n, network.rates.len());
    for (j, reaction) in network.rates.iter().enumerate() {
        s[(reaction.from, j)] -= 1.0;
        s[(reaction.to, j)] += 1.0;
    }
    s
}

fn thermochemical_reference(network: &Network, gamma: f64) -> ThermoReference {
    let n = network.species_count();
    let reaction_count = network.rates.len();
    let s = stoichiometry(network);
    let h = enthalpies(network, gamma);
    let q_reaction = -s.transpose() * &h;
    let mut r_c = DMatrix::zeros(reaction_count, n);
    for (j, reaction) in network.rates.iter().enumerate() {
        r_c[(j, reaction.from)] = reaction.rate * network.kinetic_scale * DENSITY / MOLAR_MASS_KG;
    }

    let molar_masses = DVector::from_element(n, MOLAR_MASS_KG);
    let source = DMatrix::from_diagonal(&molar_masses) * &s * &r_c;
    let mut operator = DMatrix::zeros(2 * n + 2, 2 * n + 2);
    let c_start = 2;
    let j_start = 2 + n;
    let pe_q = MASS_FLUX * CP_J_PER_KG_K * LENGTH_M / LAMBDA;
    operator[(0, 1)] = 1.0 / LAMBDA;
    operator[(1, 1)] = pe_q;
    for i in 0..n {
        let q_coefficient: f64 = (0..reaction_count)
            .map(|j| q_reaction[j] * r_c[(j, i)])
            .sum();
        operator[(1, c_start + i)] = -LENGTH_M.powi(2) / T_SCALE_K * q_coefficient;

        let k_i = DENSITY * network.diffusion[i];
        let pe_i = MASS_FLUX * LENGTH_M / k_i;
        operator[(c_start + i, j_start + i)] = 1.0 / k_i;
        operator[(j_start + i, j_start + i)] = pe_i;
        for j in 0..n {
            operator[(j_start + i, c_start + j)] = -LENGTH_M.powi(2) * source[(i, j)];
        }
    }

    let transition = matrix_exp(&operator);
    let mut known_initial = DVector::zeros(2 * n + 2);
    for (i, fraction) in network.inlet.iter().enumerate() {
        known_initial[c_start + i] = *fraction;
    }
    let unknown_indices: Vec<usize> = std::iter::once(1).chain(j_start..j_start + n).collect();
    let known_indices = unknown_indices.clone();
    let mut boundary_matrix = DMatrix::zeros(n + 1, n + 1);
    let mut boundary_rhs = DVector::zeros(n + 1);
    let q_right_internal = Q_RIGHT_PHYSICAL * LENGTH_M / T_SCALE_K;
    for (row, &state_row) in known_indices.iter().enumerate() {
        boundary_rhs[row] = if row == 0 { q_right_internal } else { 0.0 }
            - (0..operator.nrows())
                .map(|col| transition[(state_row, col)] * known_initial[col])
                .sum::<f64>();
        for (col, &state_col) in unknown_indices.iter().enumerate() {
            boundary_matrix[(row, col)] = transition[(state_row, state_col)];
        }
    }

    let singular_values = boundary_matrix.clone().svd(false, false).singular_values;
    let sigma_max = singular_values.iter().copied().fold(0.0_f64, f64::max);
    let sigma_min = singular_values
        .iter()
        .copied()
        .fold(f64::INFINITY, f64::min);
    assert!(sigma_min.is_finite() && sigma_min > 0.0);
    let condition_estimate = sigma_max / sigma_min;
    assert!(condition_estimate.is_finite());
    let unknown = boundary_matrix
        .clone()
        .lu()
        .solve(&boundary_rhs)
        .expect("thermochemical mixed-boundary system should be solvable");
    let linear_residual = &boundary_matrix * &unknown - &boundary_rhs;
    let relative_residual = linear_residual.norm() / (1.0 + boundary_rhs.norm());
    assert!(
        relative_residual < 1e-12,
        "boundary linear residual={relative_residual:e}"
    );

    let mut initial = known_initial;
    for (index, value) in unknown_indices.iter().zip(unknown.iter()) {
        initial[*index] = *value;
    }
    let right_state = &transition * &initial;
    assert!((initial[0]).abs() < 1e-14);
    assert!((right_state[1] - q_right_internal).abs() < 1e-10);
    assert!(right_state.rows(j_start, n).norm() < 1e-10);
    ThermoReference {
        stoich: s,
        reaction_rates: r_c,
        reaction_heats: q_reaction,
        source,
        operator,
        initial,
        gamma,
        condition_estimate,
    }
}

fn reference_state(reference: &ThermoReference, z: f64) -> DVector<f64> {
    crate::SimpleReactorsIVP::combustion_validation::matrix_exponential::expm_action(
        &(&reference.operator * z),
        &reference.initial,
    )
}

fn interleaved_state(block: &DVector<f64>, n: usize) -> DVector<f64> {
    let mut state = DVector::zeros(2 * n + 2);
    state[0] = block[0];
    state[1] = block[1];
    for i in 0..n {
        state[2 + 2 * i] = block[2 + i];
        state[3 + 2 * i] = block[2 + n + i];
    }
    state
}

fn block_state(interleaved: &DVector<f64>, n: usize) -> DVector<f64> {
    let mut state = DVector::zeros(2 * n + 2);
    state[0] = interleaved[0];
    state[1] = interleaved[1];
    for i in 0..n {
        state[2 + i] = interleaved[2 + 2 * i];
        state[2 + n + i] = interleaved[3 + 2 * i];
    }
    state
}

fn physical_energy_terms(
    temperature: f64,
    heat_flux: f64,
    fractions: &[f64],
    species_fluxes: &[f64],
    gamma: f64,
) -> (f64, f64) {
    let mut invariant = heat_flux - MASS_FLUX * CP_J_PER_KG_K * temperature;
    let mut scale = heat_flux.abs() + (MASS_FLUX * CP_J_PER_KG_K * temperature).abs();
    for (i, h) in H_BASE_J_PER_MOL.iter().take(fractions.len()).enumerate() {
        let g = gamma * h / MOLAR_MASS_KG;
        let term = g * (species_fluxes[i] - MASS_FLUX * fractions[i]);
        invariant += term;
        scale += term.abs();
    }
    (invariant, scale)
}

fn internal_energy_terms(state: &DVector<f64>, n: usize, gamma: f64) -> (f64, f64) {
    let convective = MASS_FLUX * CP_J_PER_KG_K * LENGTH_M * state[0];
    let mut invariant = state[1] - convective;
    let mut scale = state[1].abs() + convective.abs();
    for i in 0..n {
        let enthalpy_per_mass = gamma * H_BASE_J_PER_MOL[i] / MOLAR_MASS_KG;
        let term = enthalpy_per_mass * (state[2 + n + i] - MASS_FLUX * LENGTH_M * state[2 + i])
            / T_SCALE_K;
        invariant += term;
        scale += term.abs();
    }
    (invariant, scale)
}

fn make_reactor(network: &Network, gamma: f64) -> SimpleReactorTask {
    let mut reactor = base_fixture(network);
    let h = enthalpies(network, gamma);
    let s = stoichiometry(network);
    reactor.thermal_effects = (-s.transpose() * h).as_slice().to_vec();
    reactor
}

fn solve(network: &Network, gamma: f64, route: Route, nodes: usize) -> SimpleReactorTask {
    solve_with_frontend(
        network,
        gamma,
        route,
        nodes,
        ReactorBvpSymbolicBackend::ExprLegacy,
    )
}

fn solve_correction_newton(network: &Network, gamma: f64, nodes: usize) -> SimpleReactorTask {
    let mut reactor = make_reactor(network, gamma);
    reactor.setup_bvp().expect("BVP-03 CorrectionNewton setup");
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
    let guess = DMatrix::from_fn(2 + 2 * network.species_count(), nodes, |row, _| match row {
        0 | 1 => 0.0,
        _ if row % 2 == 0 => network.inlet[(row - 2) / 2],
        _ => 0.0,
    });
    reactor
        .solver
        .solve_with_hybrid_settings(
            BvpSolverSettings::CorrectionNewton(options),
            guess,
            0.0,
            1.0,
            nodes,
        )
        .expect("BVP-03 CorrectionNewton should solve thermochemical fixture");
    reactor.postprocessing().expect("BVP-03 publication");
    reactor
}

fn solve_with_frontend(
    network: &Network,
    gamma: f64,
    route: Route,
    nodes: usize,
    symbolic_backend: ReactorBvpSymbolicBackend,
) -> SimpleReactorTask {
    solve_with_grid(
        network,
        gamma,
        route,
        nodes,
        symbolic_backend,
        GridRefinementMethod::DoublePoints,
        2,
    )
}

fn solve_with_grid(
    network: &Network,
    gamma: f64,
    route: Route,
    nodes: usize,
    symbolic_backend: ReactorBvpSymbolicBackend,
    grid_method: GridRefinementMethod,
    max_refinements: usize,
) -> SimpleReactorTask {
    let mut reactor = make_reactor(network, gamma);
    reactor.setup_bvp().expect("BVP-03 production setup");
    reactor.solver.arg_name = "z".to_string();
    reactor.solver.x_range = (0.0, 1.0);
    reactor.solver = reactor
        .solver
        .with_bvp_telemetry_mode(BvpApiTelemetryMode::Counters);
    let guess = DMatrix::from_fn(2 + 2 * network.species_count(), nodes, |row, _| match row {
        0 => 0.0,
        1 => 0.0,
        _ if row % 2 == 0 => network.inlet[(row - 2) / 2],
        _ => 0.0,
    });
    match route {
        Route::Damped => {
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
                "forward".to_string(),
                "Damped".to_string(),
                Some(params),
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
                .expect("Damped BVP-03 solve");
        }
        Route::Frozen => reactor
            .solver
            .solve_BVPfrozen(
                guess,
                FrozenSolverOptions::default()
                    .with_tolerance(1e-9)
                    .with_max_iterations(500),
            )
            .expect("Frozen BVP-03 solve"),
        Route::Sci => reactor
            .solver
            .solve_BVPsci(
                guess,
                nodes,
                if network.species_count() > 4 {
                    128
                } else {
                    512
                },
                1e-8,
            )
            .expect("BVP_sci BVP-03 solve"),
    }
    reactor.postprocessing().expect("BVP-03 publication");
    let solution = reactor.solver.get_solution().expect("published solution");
    let last = solution.nrows() - 1;
    assert!((solution[(0, 0)] - T0_K).abs() < 1e-7);
    for (i, inlet) in network.inlet.iter().enumerate() {
        assert!((solution[(0, 2 + 2 * i)] - inlet).abs() < 1e-7);
        assert!(solution[(last, 3 + 2 * i)].abs() < 1e-7);
    }
    assert!((solution[(last, 1)] - Q_RIGHT_PHYSICAL).abs() < 1e-6);
    reactor
}

fn profile_errors(
    network: &Network,
    reference: &ThermoReference,
    reactor: &SimpleReactorTask,
) -> Errors {
    let solution = reactor.solver.get_solution().expect("published solution");
    let mesh = reactor.solver.x_mesh.as_ref().expect("published mesh");
    let n = network.species_count();
    let zero_heat_reference = thermochemical_reference(network, 0.0);
    let mut errors = Errors::default();
    let mut invariant_initial = None;
    let mut quadrature = 0.0;
    let mut previous_z = None;
    let mut previous_qvol = None;

    for (node, &x) in mesh.iter().enumerate() {
        let z = x / LENGTH_M;
        let exact = reference_state(reference, z);
        let zero_heat = reference_state(&zero_heat_reference, z);
        let t_exact = T0_K + T_SCALE_K * exact[0];
        let q_exact = exact[1] * T_SCALE_K / LENGTH_M;
        let t_error = (solution[(node, 0)] - t_exact).abs();
        assert!(
            (0.0..1000.0).contains(&solution[(node, 0)]),
            "BVP-03 published temperature outside physical fixture envelope: {} K",
            solution[(node, 0)]
        );
        if t_error > errors.t {
            errors.t = t_error;
            errors.t_z = z;
        }
        errors.t_relative = errors
            .t_relative
            .max(t_error / (t_exact - T0_K).abs().max(1.0));
        let exact_chemical_shift = T_SCALE_K * (exact[0] - zero_heat[0]);
        let actual_chemical_shift = solution[(node, 0)] - (T0_K + T_SCALE_K * zero_heat[0]);
        errors.chemical_temperature_relative = errors.chemical_temperature_relative.max(
            (actual_chemical_shift - exact_chemical_shift).abs()
                / exact_chemical_shift.abs().max(1.0),
        );
        let q_error = (solution[(node, 1)] - q_exact).abs();
        if q_error > errors.q {
            errors.q = q_error;
            errors.q_z = z;
        }
        errors.q_relative = errors.q_relative.max(q_error / q_exact.abs().max(1.0));
        let mut fractions = Vec::with_capacity(n);
        let mut fluxes = Vec::with_capacity(n);
        let mut c_block = DVector::zeros(n);
        let mut total_fraction = 0.0;
        for i in 0..n {
            let c = solution[(node, 2 + 2 * i)];
            let j = solution[(node, 3 + 2 * i)];
            assert!(
                c.is_finite() && c >= -1e-8,
                "invalid C[{i}]={c} at node {node}"
            );
            assert!(j.is_finite());
            fractions.push(c);
            fluxes.push(j);
            total_fraction += c;
            c_block[i] = c;
            let c_error = (c - exact[2 + i]).abs();
            if c_error > errors.c {
                errors.c = c_error;
                errors.c_z = z;
            }
            errors.c_relative = errors
                .c_relative
                .max(c_error / exact[2 + i].abs().max(1e-6));
            let j_error = (j - exact[2 + n + i] / LENGTH_M).abs();
            if j_error > errors.j {
                errors.j = j_error;
                errors.j_z = z;
            }
            errors.j_relative = errors
                .j_relative
                .max(j_error / (exact[2 + n + i] / LENGTH_M).abs().max(1e-6));
        }
        errors.composition_sum_deviation = errors
            .composition_sum_deviation
            .max((total_fraction - 1.0).abs());

        let (energy, energy_scale) = physical_energy_terms(
            solution[(node, 0)],
            solution[(node, 1)],
            &fractions,
            &fluxes,
            reference.gamma,
        );
        let baseline = *invariant_initial.get_or_insert(energy);
        let residual = energy - baseline;
        if residual.abs() > errors.energy_abs {
            errors.energy_abs = residual.abs();
            errors.energy_z = z;
        }
        errors.energy_rel = errors
            .energy_rel
            .max(residual.abs() / energy_scale.max(1.0));

        let qvol = (&reference.reaction_rates * &c_block).dot(&reference.reaction_heats);
        if let (Some(old_z), Some(old_qvol)) = (previous_z, previous_qvol) {
            quadrature += 0.5 * (qvol + old_qvol) * (z - old_z);
        }
        previous_z = Some(z);
        previous_qvol = Some(qvol);
    }

    let left_flux_sum: f64 = (0..n)
        .map(|i| {
            let c = solution[(0, 2 + 2 * i)];
            let j = solution[(0, 3 + 2 * i)];
            (reference.gamma * H_BASE_J_PER_MOL[i] / MOLAR_MASS_KG) * (j - MASS_FLUX * c)
        })
        .sum();
    let last = solution.nrows() - 1;
    let right_flux_sum: f64 = (0..n)
        .map(|i| {
            let c = solution[(last, 2 + 2 * i)];
            let j = solution[(last, 3 + 2 * i)];
            (reference.gamma * H_BASE_J_PER_MOL[i] / MOLAR_MASS_KG) * (j - MASS_FLUX * c)
        })
        .sum();
    let flux_integral = (right_flux_sum - left_flux_sum) / LENGTH_M;
    let integral_residual = quadrature - flux_integral;
    errors.heat_integral_abs = integral_residual.abs();
    errors.heat_integral_rel =
        integral_residual.abs() / quadrature.abs().max(flux_integral.abs()).max(1.0);
    errors
}

fn max_chemical_temperature_shift(network: &Network, gamma: f64) -> f64 {
    let cold = thermochemical_reference(network, 0.0);
    let hot = thermochemical_reference(network, gamma);
    (0..=128)
        .map(|i| {
            let z = i as f64 / 128.0;
            T_SCALE_K * (reference_state(&hot, z)[0] - reference_state(&cold, z)[0]).abs()
        })
        .fold(0.0_f64, f64::max)
}

fn check_hess_and_operator(network: &Network, reference: &ThermoReference) {
    let h = enthalpies(network, reference.gamma);
    let hess_q = -reference.stoich.transpose() * h;
    assert!((&hess_q - &reference.reaction_heats).norm() < 1e-12);
    assert!(reference.reaction_heats.iter().any(|q| q.abs() > 0.0));
    assert!(
        (&reference.source - network.source_matrix()).norm() < 1e-11,
        "BVP-02/03 independently assembled species source matrices disagree"
    );
    println!(
        "BVP-03 Hess species={} gamma={} reaction_heats_J_per_mol={:?}",
        network.species_count(),
        reference.gamma,
        reference.reaction_heats
    );
    for col in 0..reference.source.ncols() {
        let stoich_sum: f64 = (0..reference.stoich.nrows())
            .map(|row| reference.stoich[(row, col)])
            .sum();
        assert!(
            stoich_sum.abs() < 1e-14,
            "reaction {col} atom-count surrogate={stoich_sum:e}"
        );
        let mass_sum: f64 = (0..reference.source.nrows())
            .map(|row| reference.source[(row, col)])
            .sum();
        assert!(
            mass_sum.abs() < 1e-12,
            "mass-source column {col}={mass_sum:e}"
        );
    }
    if network.species_count() == 12 {
        let cycle_heat: f64 = (0..12).map(|j| reference.reaction_heats[j]).sum();
        assert!(cycle_heat.abs() < 1e-12, "closed-cycle heat={cycle_heat:e}");
    }
}

fn assert_rhs_and_jacobian(network: &Network, gamma: f64, reference: &ThermoReference) {
    let mut reactor = make_reactor(network, gamma);
    reactor
        .setup_bvp()
        .expect("BVP-03 derivative fixture setup");
    let n = network.species_count();
    let block = reference_state(reference, 0.37);
    let state = interleaved_state(&block, n);
    let expected_block = &reference.operator * &block;
    let expected = interleaved_state(&expected_block, n);
    let mut permutation = DMatrix::zeros(state.len(), state.len());
    for col in 0..state.len() {
        let mut basis = DVector::zeros(state.len());
        basis[col] = 1.0;
        let block_basis = block_state(&basis, n);
        let block_derivative = &reference.operator * block_basis;
        let interleaved_derivative = interleaved_state(&block_derivative, n);
        permutation.set_column(col, &interleaved_derivative);
    }
    let variables: Vec<&str> = reactor.solver.unknowns.iter().map(String::as_str).collect();
    for (row, equation) in reactor.solver.eq_system.iter().enumerate() {
        let actual = equation.eval_expression(&variables, state.as_slice());
        assert!(
            (actual - expected[row]).abs() <= 5e-10 * (1.0 + expected[row].abs()),
            "BVP-03 RHS row={row}: actual={actual:e}, expected={:e}",
            expected[row]
        );
        for (col, variable) in variables.iter().enumerate() {
            let derivative = equation
                .diff(variable)
                .eval_expression(&variables, state.as_slice());
            assert!(
                (derivative - permutation[(row, col)]).abs()
                    <= 5e-10 * (1.0 + permutation[(row, col)].abs()),
                "BVP-03 Jacobian [{row},{col}] {derivative:e} != {:e}",
                permutation[(row, col)]
            );
            if col == 1 || col == 2 || col == 3 + 2 * (n / 2) {
                let step = 1e-7;
                let mut shifted = state.clone();
                shifted[col] += step;
                let numerical =
                    (equation.eval_expression(&variables, shifted.as_slice()) - actual) / step;
                assert!(
                    (derivative - numerical).abs() < 5e-5 * (1.0 + derivative.abs()),
                    "finite-difference Jacobian [{row},{col}] {derivative:e} != {numerical:e}"
                );
            }
        }
    }
}

fn assert_route(network: &Network, gamma: f64, route: Route) -> Errors {
    let nodes = match route {
        Route::Damped => 64,
        Route::Frozen => 1025,
        Route::Sci if network.species_count() > 4 => 17,
        Route::Sci => 33,
    };
    let reference = thermochemical_reference(network, gamma);
    let reactor = solve(network, gamma, route, nodes);
    let errors = profile_errors(network, &reference, &reactor);
    println!(
        "BVP-03 route={route:?} species={} reactions={} gamma={gamma} nodes={} cond2={:.3e} dC={:.3e}@{:.3}/rel={:.3e} dJ={:.3e}@{:.3}/rel={:.3e} dT={:.3e}K@{:.3}/rel={:.3e} dTchem_rel={:.3e} dq={:.3e}W/m2@{:.3}/rel={:.3e} sumC_deviation={:.3e} energy_abs={:.3e} energy_rel={:.3e} energy_z={:.3} heat_integral_abs={:.3e} heat_integral_rel={:.3e}",
        network.species_count(),
        network.rates.len(),
        reactor.solver.x_mesh.as_ref().unwrap().len(),
        reference.condition_estimate,
        errors.c,
        errors.c_z,
        errors.c_relative,
        errors.j,
        errors.j_z,
        errors.j_relative,
        errors.t,
        errors.t_z,
        errors.t_relative,
        errors.chemical_temperature_relative,
        errors.q,
        errors.q_z,
        errors.q_relative,
        errors.composition_sum_deviation,
        errors.energy_abs,
        errors.energy_rel,
        errors.energy_z,
        errors.heat_integral_abs,
        errors.heat_integral_rel
    );
    let (c_tol, j_tol, t_tol, q_tol) = match route {
        Route::Damped => (3e-4, 1e-5, 8e-2, 8e-1),
        Route::Frozen => (1e-4, 1e-5, 3e-2, 3e-1),
        Route::Sci => (2e-5, 2e-6, 2e-4, 2e-3),
    };
    assert!(errors.c < c_tol, "{route:?} composition error={errors:?}");
    assert!(errors.j < j_tol, "{route:?} species-flux error={errors:?}");
    assert!(errors.t < t_tol, "{route:?} temperature error={errors:?}");
    assert!(
        errors.chemical_temperature_relative < 2e-2,
        "{route:?} chemical temperature response error={errors:?}"
    );
    assert!(errors.q < q_tol, "{route:?} heat-flux error={errors:?}");
    assert!(
        errors.energy_rel < 2e-5,
        "{route:?} energy invariant={errors:?}"
    );
    assert!(
        errors.heat_integral_rel < 2e-3,
        "{route:?} heat integral={errors:?}"
    );
    if network
        .diffusion
        .iter()
        .all(|value| *value == network.diffusion[0])
    {
        assert!(
            errors.composition_sum_deviation < 1e-3,
            "composition sum={errors:?}"
        );
    }
    errors
}

#[test]
fn bvp03_thermal_scale_calibration_is_linear_and_well_conditioned() {
    for network in [Network::small(), Network::large()] {
        let mut previous = 0.0;
        let unit_response = max_chemical_temperature_shift(&network, 1.0);
        for gamma in [0.0, 10.0, 100.0, 300.0, 1_000.0, 10_000.0] {
            let delta_t = max_chemical_temperature_shift(&network, gamma);
            let reference = thermochemical_reference(&network, gamma);
            println!(
                "BVP-03 calibration species={} gamma={gamma:.1e} delta_T_chem={delta_t:.6e} K cond2={:.6e}",
                network.species_count(),
                reference.condition_estimate
            );
            assert!(delta_t.is_finite());
            assert!(delta_t + 1e-12 >= previous);
            assert!((delta_t - gamma * unit_response).abs() < 2e-10 * (1.0 + delta_t));
            previous = delta_t;
        }
        let selected = max_chemical_temperature_shift(&network, GAMMA);
        let selected_reference = thermochemical_reference(&network, GAMMA);
        assert!(
            (5.0..=30.0).contains(&selected),
            "selected gamma={GAMMA} gives chemical temperature shift {selected} K"
        );
        assert!(selected_reference.condition_estimate < 1e10);
    }
}

fn exact_reference_heat_quadrature_error(network: &Network, gamma: f64, intervals: usize) -> f64 {
    let reference = thermochemical_reference(network, gamma);
    let n = network.species_count();
    let weighted_transport = |state: &DVector<f64>| {
        (0..n)
            .map(|i| {
                gamma * H_BASE_J_PER_MOL[i] / MOLAR_MASS_KG
                    * (state[2 + n + i] - MASS_FLUX * LENGTH_M * state[2 + i])
            })
            .sum::<f64>()
    };
    let mut quadrature = 0.0;
    let mut previous: Option<(f64, f64)> = None;
    for index in 0..=intervals {
        let z = index as f64 / intervals as f64;
        let state = reference_state(&reference, z);
        let c = DVector::from_iterator(n, (0..n).map(|i| state[2 + i]));
        let q_vol = (&reference.reaction_rates * c).dot(&reference.reaction_heats);
        if let Some((previous_z, previous_q)) = previous {
            quadrature += 0.5 * (q_vol + previous_q) * (z - previous_z);
        }
        previous = Some((z, q_vol));
    }
    let left = weighted_transport(&reference_state(&reference, 0.0));
    let right = weighted_transport(&reference_state(&reference, 1.0));
    let flux_identity = (right - left) / LENGTH_M.powi(2);
    (quadrature - flux_identity).abs() / flux_identity.abs().max(1.0)
}

#[test]
fn bvp03_integrated_heat_quadrature_converges_to_the_flux_identity() {
    let network = Network::large();
    let errors: Vec<f64> = [32, 128, 512]
        .into_iter()
        .map(|intervals| {
            let error = exact_reference_heat_quadrature_error(&network, GAMMA, intervals);
            println!(
                "BVP-03 exact heat quadrature intervals={intervals} relative_error={error:.6e}"
            );
            error
        })
        .collect();
    assert!(errors.windows(2).all(|pair| pair[1] < pair[0]));
    assert!(*errors.last().unwrap() < 1e-6);
}

#[test]
fn bvp03a_damped_and_frozen_mesh_refinement_reduces_thermochemical_error() {
    let network = Network::small();
    let reference = thermochemical_reference(&network, GAMMA);
    for (route, meshes) in [
        (Route::Damped, vec![32, 64, 128]),
        (Route::Frozen, vec![129, 257, 513, 1025]),
    ] {
        let errors: Vec<Errors> = meshes
            .into_iter()
            .map(|nodes| {
                let reactor = solve(&network, GAMMA, route, nodes);
                let result = profile_errors(&network, &reference, &reactor);
                println!(
                    "BVP-03 mesh route={route:?} seed={nodes} published={} dC={:.3e} dJ={:.3e} dT={:.3e} K dq={:.3e} W/m2",
                    reactor.solver.x_mesh.as_ref().unwrap().len(),
                    result.c,
                    result.j,
                    result.t,
                    result.q
                );
                result
            })
            .collect();
        assert!(errors.windows(2).all(|pair| {
            pair[1].c < pair[0].c
                && pair[1].j < pair[0].j
                && pair[1].t < pair[0].t
                && pair[1].q < pair[0].q
        }));
    }
}

#[test]
fn bvp03a_damped_refinement_parameter_sweep_measures_thermochemical_accuracy() {
    let network = Network::small();
    let reference = thermochemical_reference(&network, GAMMA);
    let seed_nodes = 33;
    let baseline = solve_with_grid(
        &network,
        GAMMA,
        Route::Damped,
        seed_nodes,
        ReactorBvpSymbolicBackend::ExprLegacy,
        GridRefinementMethod::DoublePoints,
        0,
    );
    let baseline_errors = profile_errors(&network, &reference, &baseline);
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
        let reactor = solve_with_grid(
            &network,
            GAMMA,
            Route::Damped,
            seed_nodes,
            ReactorBvpSymbolicBackend::ExprLegacy,
            grid_method,
            3,
        );
        let errors = profile_errors(&network, &reference, &reactor);
        let nodes = reactor.solver.x_mesh.as_ref().unwrap().len();
        println!(
            "BVP-03a tuned grid={name} nodes={nodes} dC={:.3e} dJ={:.3e} dT={:.3e}K dq={:.3e}W/m2 energy_rel={:.3e}",
            errors.c, errors.j, errors.t, errors.q, errors.energy_rel,
        );
        measured.push((name, nodes, errors));
    }
    println!(
        "BVP-03a fixed seed nodes={seed_nodes} dC={:.3e} dJ={:.3e} dT={:.3e}K dq={:.3e}W/m2 energy_rel={:.3e}",
        baseline_errors.c,
        baseline_errors.j,
        baseline_errors.t,
        baseline_errors.q,
        baseline_errors.energy_rel,
    );
    let pearson = measured
        .iter()
        .find(|row| row.0 == "Pearson-1e-2")
        .expect("the Pearson thermochemical candidate must be measured");
    assert_eq!(pearson.1, 170);
    assert!(
        pearson.2.c < baseline_errors.c
            && pearson.2.j < baseline_errors.j
            && pearson.2.t < baseline_errors.t
            && pearson.2.q < baseline_errors.q,
        "selected Pearson settings must improve every thermochemical field: baseline={baseline_errors:?}, candidate={pearson:?}"
    );
    assert!(
        measured.iter().any(|(_, _, errors)| {
            errors.c < baseline_errors.c
                && errors.j < baseline_errors.j
                && errors.t < baseline_errors.t
                && errors.q < baseline_errors.q
        }),
        "no tuned strategy improved all thermochemical fields: baseline={baseline_errors:?}, candidates={measured:?}"
    );
}

#[test]
fn bvp03a_correction_newton_matches_thermochemical_oracle_and_energy_invariant() {
    let network = Network::small();
    let reference = thermochemical_reference(&network, GAMMA);
    let coarse = solve_correction_newton(&network, GAMMA, 129);
    let coarse_errors = profile_errors(&network, &reference, &coarse);
    let fine = solve_correction_newton(&network, GAMMA, 257);
    let fine_errors = profile_errors(&network, &reference, &fine);
    println!(
        "BVP-03a CorrectionNewton fixed mesh coarse=129 dC={:.3e} dJ={:.3e} dT={:.3e}K dq={:.3e}W/m2 energy_rel={:.3e} fine=257 dC={:.3e} dJ={:.3e} dT={:.3e}K dq={:.3e}W/m2 energy_rel={:.3e}",
        coarse_errors.c,
        coarse_errors.j,
        coarse_errors.t,
        coarse_errors.q,
        coarse_errors.energy_rel,
        fine_errors.c,
        fine_errors.j,
        fine_errors.t,
        fine_errors.q,
        fine_errors.energy_rel,
    );
    assert!(
        fine_errors.c < coarse_errors.c
            && fine_errors.j < coarse_errors.j
            && fine_errors.t < coarse_errors.t
            && fine_errors.q < coarse_errors.q,
        "CorrectionNewton refinement should improve all thermochemical profile errors: coarse={coarse_errors:?}, fine={fine_errors:?}"
    );
    assert!(fine_errors.energy_rel < 1e-12);
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
fn bvp03a_small_network_checks_thermochemistry_rhs_jacobian_and_routes() {
    let network = Network::small();
    let reference = thermochemical_reference(&network, GAMMA);
    assert_eq!(reference.operator.shape(), (10, 10));
    check_hess_and_operator(&network, &reference);
    let expected_heats = [40.0, -15.0, 40.0, 25.0, 25.0].map(|q| q * GAMMA);
    for (actual, expected) in reference.reaction_heats.iter().zip(expected_heats) {
        assert!((actual - expected).abs() < 1e-12);
    }
    assert_rhs_and_jacobian(&network, GAMMA, &reference);
    for route in [Route::Damped, Route::Frozen, Route::Sci] {
        assert_route(&network, GAMMA, route);
    }

    let expr = solve_with_frontend(
        &network,
        GAMMA,
        Route::Damped,
        64,
        ReactorBvpSymbolicBackend::ExprLegacy,
    );
    let atom = solve_with_frontend(
        &network,
        GAMMA,
        Route::Damped,
        64,
        ReactorBvpSymbolicBackend::AtomView,
    );
    let expr_solution = expr.solver.get_solution().unwrap();
    let atom_solution = atom.solver.get_solution().unwrap();
    assert_eq!(expr_solution.shape(), atom_solution.shape());
    let frontend_delta = expr_solution
        .iter()
        .zip(atom_solution.iter())
        .map(|(a, b)| (a - b).abs())
        .fold(0.0_f64, f64::max);
    println!("BVP-03a ExprLegacy/AtomView max profile delta={frontend_delta:.3e}");
    assert!(frontend_delta < 2e-8);
}

#[test]
fn bvp03b_large_network_checks_all_reaction_heats_and_solver_routes() {
    let network = Network::large();
    let reference = thermochemical_reference(&network, GAMMA);
    assert_eq!(reference.operator.shape(), (26, 26));
    assert_eq!(reference.reaction_heats.len(), 20);
    check_hess_and_operator(&network, &reference);
    assert_rhs_and_jacobian(&network, GAMMA, &reference);
    for route in [Route::Damped, Route::Frozen, Route::Sci] {
        assert_route(&network, GAMMA, route);
    }
}

#[test]
fn bvp03c_unequal_diffusion_preserves_energy_without_fraction_normalization() {
    let mut network = Network::large();
    network.diffusion = (0..network.species_count())
        .map(|i| 1e-4 * (1.0 + i as f64 / 10.0))
        .collect();
    let reference = thermochemical_reference(&network, GAMMA);
    check_hess_and_operator(&network, &reference);
    let errors = assert_route(&network, GAMMA, Route::Damped);
    assert!(errors.energy_rel < 2e-5);
    assert!(
        errors.composition_sum_deviation > 1e-6,
        "unequal-diffusion solution unexpectedly appears composition-normalized: {errors:?}"
    );
}

#[test]
fn bvp03d_zero_heat_limit_recovers_species_solution_and_temperature_response_is_linear() {
    let network = Network::large();
    let reference_zero = thermochemical_reference(&network, 0.0);
    let reference_one = thermochemical_reference(&network, 1.0);
    for gamma in [0.0, 10.0, 100.0] {
        let reference = thermochemical_reference(&network, gamma);
        for index in 0..reference.initial.len() {
            let is_species_coordinate = index >= 2;
            if is_species_coordinate {
                assert!(
                    (reference.initial[index] - reference_zero.initial[index]).abs() < 1e-12,
                    "species initial state changed with gamma at index={index}"
                );
            }
        }
        for sample in 0..=32 {
            let z = sample as f64 / 32.0;
            let state = reference_state(&reference, z);
            let zero = reference_state(&reference_zero, z);
            let unit = reference_state(&reference_one, z);
            for index in 2..state.len() {
                assert!((state[index] - zero[index]).abs() < 2e-11);
            }
            let predicted_theta = zero[0] + gamma * (unit[0] - zero[0]);
            let predicted_q = zero[1] + gamma * (unit[1] - zero[1]);
            assert!((state[0] - predicted_theta).abs() < 2e-10);
            assert!((state[1] - predicted_q).abs() < 2e-10);
        }
    }

    let hot = thermochemical_reference(&network, GAMMA);
    let n = network.species_count();
    let first_state = reference_state(&hot, 0.0);
    let (baseline_internal, _) = internal_energy_terms(&first_state, n, GAMMA);
    let first_c: Vec<f64> = (0..n).map(|i| first_state[2 + i]).collect();
    let first_j: Vec<f64> = (0..n).map(|i| first_state[2 + n + i] / LENGTH_M).collect();
    let (baseline_dimensional, _) = physical_energy_terms(
        T0_K + T_SCALE_K * first_state[0],
        first_state[1] * T_SCALE_K / LENGTH_M,
        &first_c,
        &first_j,
        GAMMA,
    );
    let mut max_internal_invariant = 0.0_f64;
    let mut max_internal_z = 0.0_f64;
    let mut max_internal_relative = 0.0_f64;
    let mut max_internal_relative_z = 0.0_f64;
    let mut max_dimensional_invariant = 0.0_f64;
    let mut max_dimensional_z = 0.0_f64;
    let mut max_dimensional_relative = 0.0_f64;
    let mut max_dimensional_relative_z = 0.0_f64;
    let mut max_internal_scale = 0.0_f64;
    let mut max_dimensional_scale = 0.0_f64;
    for sample in 1..=128 {
        let z = sample as f64 / 128.0;
        let state = reference_state(&hot, z);
        let (internal, internal_scale) = internal_energy_terms(&state, n, GAMMA);
        let c: Vec<f64> = (0..n).map(|i| state[2 + i]).collect();
        let internal_j: Vec<f64> = (0..n).map(|i| state[2 + n + i]).collect();
        let physical_j: Vec<f64> = internal_j.iter().map(|value| value / LENGTH_M).collect();
        let physical_t = T0_K + T_SCALE_K * state[0];
        let physical_q = state[1] * T_SCALE_K / LENGTH_M;
        let (dimensional, _) =
            physical_energy_terms(physical_t, physical_q, &c, &physical_j, GAMMA);
        let internal_abs = (internal - baseline_internal).abs();
        if internal_abs > max_internal_invariant {
            max_internal_invariant = internal_abs;
            max_internal_z = z;
        }
        let dimensional_abs = (dimensional - baseline_dimensional).abs();
        if dimensional_abs > max_dimensional_invariant {
            max_dimensional_invariant = dimensional_abs;
            max_dimensional_z = z;
        }
        max_internal_scale = max_internal_scale.max(internal_scale);
        let internal_relative = internal_abs / internal_scale.max(1e-12);
        if internal_relative > max_internal_relative {
            max_internal_relative = internal_relative;
            max_internal_relative_z = z;
        }
        let scale = physical_energy_terms(physical_t, physical_q, &c, &physical_j, GAMMA).1;
        max_dimensional_scale = max_dimensional_scale.max(scale);
        let dimensional_relative = dimensional_abs / scale.max(1.0);
        if dimensional_relative > max_dimensional_relative {
            max_dimensional_relative = dimensional_relative;
            max_dimensional_relative_z = z;
        }
    }
    println!(
        "BVP-03d analytic energy invariant internal_abs={max_internal_invariant:.3e}@{max_internal_z:.3} internal_rel={max_internal_relative:.3e}@{max_internal_relative_z:.3} internal_scale_max={max_internal_scale:.3e} dimensional_abs={max_dimensional_invariant:.3e}@{max_dimensional_z:.3} dimensional_rel={max_dimensional_relative:.3e}@{max_dimensional_relative_z:.3} dimensional_scale_max={max_dimensional_scale:.3e}"
    );
    assert!(max_internal_relative < 2e-12);
    assert!(max_dimensional_relative < 2e-12);

    let baseline_reactor = solve(&network, 0.0, Route::Damped, 64);
    let reactor = solve(&network, GAMMA, Route::Damped, 64);
    let solution = reactor.solver.get_solution().expect("published solution");
    let baseline_solution = baseline_reactor
        .solver
        .get_solution()
        .expect("published zero-heat solution");
    assert_eq!(solution.shape(), baseline_solution.shape());
    let mut max_species_delta = 0.0_f64;
    for row in 0..solution.nrows() {
        for i in 0..network.species_count() {
            max_species_delta = max_species_delta
                .max((solution[(row, 2 + 2 * i)] - baseline_solution[(row, 2 + 2 * i)]).abs())
                .max((solution[(row, 3 + 2 * i)] - baseline_solution[(row, 3 + 2 * i)]).abs());
        }
    }
    println!(
        "BVP-03d production gamma=0 vs gamma={GAMMA} max species-state delta={max_species_delta:.3e}"
    );
    assert!(max_species_delta < 2e-8);
    let mesh = reactor.solver.x_mesh.as_ref().unwrap();
    let mut max_solver_relative = 0.0_f64;
    let mut baseline = None;
    for (node, &x) in mesh.iter().enumerate() {
        let fractions: Vec<f64> = (0..network.species_count())
            .map(|i| solution[(node, 2 + 2 * i)])
            .collect();
        let fluxes: Vec<f64> = (0..network.species_count())
            .map(|i| solution[(node, 3 + 2 * i)])
            .collect();
        let (value, scale) = physical_energy_terms(
            solution[(node, 0)],
            solution[(node, 1)],
            &fractions,
            &fluxes,
            GAMMA,
        );
        let base = *baseline.get_or_insert(value);
        max_solver_relative = max_solver_relative.max((value - base).abs() / scale.max(1.0));
        let exact = reference_state(&hot, x / LENGTH_M);
        assert!((solution[(node, 0)] - (T0_K + T_SCALE_K * exact[0])).abs() < 8e-2);
    }
    assert!(
        max_solver_relative < 2e-5,
        "solver energy invariant={max_solver_relative:e}"
    );
}
