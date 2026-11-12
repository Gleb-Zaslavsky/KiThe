use crate::SimpleReactorsBVP::SimpleReactorBVP::{
    DensityModel, FastElemReact, NrbvpHandoffConfig, R_G, SimpleReactorTask,
};
use crate::SimpleReactorsBVP::reactor_BVP_utils::ScalingConfig;
use crate::SimpleReactorsBVP::solver_backend::ReactorBvpSolverConfig;
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
use nalgebra::DMatrix;
use std::collections::HashMap;

const R: f64 = 8.314;
const T0_K: f64 = 300.0;
const TRANSPORT_T_K: f64 = 298.15;
const RHO: f64 = 1.0;
const MASS_FLUX: f64 = 0.01;
const CP: f64 = 1000.0;
const LAMBDA: f64 = 0.1;
const DIFFUSION: f64 = 1e-4;
const TAU: f64 = 0.5;
const A0: f64 = 1.0 / TAU;
const MOLAR_MASS: f64 = 0.060052;
const Q_MASS: f64 = 10.0 * CP * T0_K;
const Q_MOLAR: f64 = Q_MASS * MOLAR_MASS;
const DELTA_T_SCALE: f64 = Q_MASS / CP;
const ACTIVATION_ENERGY: f64 = 7.0 * R * T0_K;

#[derive(Clone, Copy, Debug)]
enum Route {
    Sci,
    Damped,
    Frozen,
}

fn make_reactor(length_m: f64) -> SimpleReactorTask {
    let mut reactor = SimpleReactorTask::new();
    reactor.set_density_model(DensityModel::Constant);
    reactor
        .fast_react_set(vec![FastElemReact {
            eq: "A=>B".into(),
            A: A0,
            n: 0.0,
            E: ACTIVATION_ENERGY,
            Q: 0.0,
        }])
        .expect("ZZ-1962-B reaction fixture");
    reactor.kindata.substances = vec!["A".into(), "B".into()];
    let formula = HashMap::from([
        ("C".to_string(), 2),
        ("H".to_string(), 4),
        ("O".to_string(), 2),
    ]);
    reactor.kindata.groups = Some(HashMap::from([
        ("A".to_string(), formula.clone()),
        ("B".to_string(), formula),
    ]));

    let pressure = RHO * R_G * TRANSPORT_T_K / MOLAR_MASS;
    let boundary = HashMap::from([
        ("T".to_string(), T0_K),
        ("A".to_string(), 1.0),
        ("B".to_string(), 0.0),
        ("q".to_string(), 0.0),
        ("A_flux".to_string(), 0.0),
        ("B_flux".to_string(), 0.0),
    ]);
    reactor.set_parameters(
        vec![Q_MOLAR],
        pressure,
        TRANSPORT_T_K,
        CP,
        boundary,
        LAMBDA,
        HashMap::from([("A".to_string(), DIFFUSION), ("B".to_string(), DIFFUSION)]),
        MASS_FLUX,
        ScalingConfig::new(T0_K, length_m, DELTA_T_SCALE),
    );
    reactor.M = MOLAR_MASS;
    reactor
}

fn solve_full(
    length_m: f64,
    seed_profile: ReferenceBranch,
    route: Route,
) -> Result<SimpleReactorTask, String> {
    let nodes = match route {
        Route::Sci => 101,
        Route::Damped => 256,
        Route::Frozen => 2049,
    };
    let guess = reference_initial_guess(length_m, seed_profile, nodes);
    solve_full_with_guess(length_m, guess, route, seed_profile.name)
}

fn solve_correction_newton(
    length_m: f64,
    seed_profile: ReferenceBranch,
    nodes: usize,
) -> Result<SimpleReactorTask, String> {
    let mut reactor = make_reactor(length_m);
    reactor
        .setup_bvp()
        .map_err(|error| format!("CorrectionNewton setup: {error}"))?;
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
            let bounds = if name == "q" || name.starts_with('J') {
                (-10.0, 10.0)
            } else {
                (-0.1, 1.1)
            };
            (name.clone(), bounds)
        })
        .collect();
    reactor
        .solver
        .solve_with_hybrid_settings(
            BvpSolverSettings::CorrectionNewton(options),
            reference_initial_guess(length_m, seed_profile, nodes),
            0.0,
            1.0,
            nodes,
        )
        .map_err(|error| {
            format!(
                "CorrectionNewton branch seed={}: {error}",
                seed_profile.name
            )
        })?;
    reactor
        .postprocessing()
        .map_err(|error| format!("CorrectionNewton postprocessing: {error}"))?;
    Ok(reactor)
}

fn reference_initial_guess(
    length_m: f64,
    seed_profile: ReferenceBranch,
    nodes: usize,
) -> DMatrix<f64> {
    let eta_max = length_m / (MASS_FLUX / RHO * TAU);
    let reference = parse_reference(seed_profile.csv);
    DMatrix::from_fn(6, nodes, |row, column| {
        let eta = eta_max * column as f64 / (nodes - 1) as f64;
        let pos = eta / eta_max * (reference.len() - 1) as f64;
        let lo = pos.floor() as usize;
        let hi = (lo + 1).min(reference.len() - 1);
        let fraction = pos - lo as f64;
        let theta = reference[lo][1] + fraction * (reference[hi][1] - reference[lo][1]);
        let theta_eta = reference[lo][2] + fraction * (reference[hi][2] - reference[lo][2]);
        reduced_state_guess(row, eta_max, theta, theta_eta)
    })
}

fn reduced_state_guess(row: usize, eta_max: f64, theta: f64, theta_eta: f64) -> f64 {
    let theta_z = eta_max * theta_eta;
    match row {
        0 => theta,
        1 => LAMBDA * theta_z,
        2 => 1.0 - theta,
        3 => -RHO * DIFFUSION * theta_z,
        4 => theta,
        5 => RHO * DIFFUSION * theta_z,
        _ => unreachable!(),
    }
}

fn solve_full_with_guess(
    length_m: f64,
    guess: DMatrix<f64>,
    route: Route,
    seed_name: &str,
) -> Result<SimpleReactorTask, String> {
    solve_full_with_limits(length_m, guess, route, seed_name, 20_000, 1e-9)
}

fn solve_full_with_limits(
    length_m: f64,
    guess: DMatrix<f64>,
    route: Route,
    seed_name: &str,
    max_nodes: usize,
    tolerance: f64,
) -> Result<SimpleReactorTask, String> {
    solve_full_with_grid(
        length_m, guess, route, seed_name, max_nodes, tolerance, None,
    )
}

fn solve_full_with_grid(
    length_m: f64,
    guess: DMatrix<f64>,
    route: Route,
    seed_name: &str,
    max_nodes: usize,
    tolerance: f64,
    damped_grid: Option<(GridRefinementMethod, usize)>,
) -> Result<SimpleReactorTask, String> {
    let mut reactor = make_reactor(length_m);
    reactor.setup_bvp().map_err(|e| format!("setup: {e}"))?;
    let nodes = guess.ncols();
    match route {
        Route::Sci => reactor
            .solver
            .solve_BVPsci(guess, nodes, max_nodes, tolerance)
            .map_err(|e| format!("BVP_sci branch seed={seed_name}: {e}"))?,
        Route::Frozen => reactor
            .solver
            .solve_BVPfrozen(
                guess,
                FrozenSolverOptions::default()
                    .with_tolerance(1e-9)
                    .with_max_iterations(500),
            )
            .map_err(|e| format!("Frozen branch seed={seed_name}: {e}"))?,
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
            let tolerances = reactor
                .solver
                .unknowns
                .iter()
                .map(|name| (name.clone(), 1e-8))
                .collect();
            let bounds = reactor
                .solver
                .unknowns
                .iter()
                .map(|name| {
                    let range = if name == "q" || name.starts_with('J') {
                        (-10.0, 10.0)
                    } else {
                        (-0.1, 1.1)
                    };
                    (name.clone(), range)
                })
                .collect();
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
                Some(tolerances),
                200,
                Some(bounds),
                Some("off".into()),
            )
            .with_solver_backend_config(ReactorBvpSolverConfig::default_lambdify());
            reactor
                .solver
                .solve_with_bvp_api(config)
                .map_err(|e| format!("Damped branch seed={seed_name}: {e}"))?;
        }
    }
    reactor
        .postprocessing()
        .map_err(|e| format!("postprocessing: {e}"))?;
    Ok(reactor)
}

fn smooth_high_branch_guess(eta_max: f64, nodes: usize, theta_r: f64, a: f64) -> DMatrix<f64> {
    let denominator = 1.0 - (-a * eta_max).exp();
    DMatrix::from_fn(6, nodes, |row, column| {
        let eta = eta_max * column as f64 / (nodes - 1) as f64;
        let theta = theta_r * (1.0 - (-a * eta).exp()) / denominator;
        let theta_eta = theta_r * a * (-a * eta).exp() / denominator;
        reduced_state_guess(row, eta_max, theta, theta_eta)
    })
}

fn prolong_published_solution(
    previous: &SimpleReactorTask,
    previous_eta_max: f64,
    new_eta_max: f64,
    nodes: usize,
) -> DMatrix<f64> {
    let solution = previous
        .solver
        .get_solution()
        .expect("previous BVP solution");
    let mesh = previous.solver.x_mesh.as_ref().expect("previous BVP mesh");
    let old_length = previous.L;
    let new_length = new_eta_max * (MASS_FLUX / RHO) * TAU;
    let previous_domain_length = previous_eta_max * (MASS_FLUX / RHO) * TAU;
    DMatrix::from_fn(6, nodes, |row, column| {
        let z = column as f64 / (nodes - 1) as f64;
        let x = z * new_length;
        let position = (x / previous_domain_length).clamp(0.0, 1.0);
        let old_x = position * old_length;
        let upper = mesh
            .as_slice()
            .partition_point(|mesh_x| *mesh_x < old_x)
            .clamp(1, mesh.len() - 1);
        let lower = upper - 1;
        let span = mesh[upper] - mesh[lower];
        let fraction = if span > 0.0 {
            ((old_x - mesh[lower]) / span).clamp(0.0, 1.0)
        } else {
            0.0
        };
        let value =
            solution[(lower, row)] + fraction * (solution[(upper, row)] - solution[(lower, row)]);
        match row {
            0 => (value - T0_K) / DELTA_T_SCALE,
            1 => value * new_length / DELTA_T_SCALE,
            2 | 4 => value,
            3 | 5 => value * new_length,
            _ => unreachable!(),
        }
    })
}

#[derive(Clone, Copy)]
struct ReferenceBranch {
    name: &'static str,
    csv: &'static str,
}

fn branches_eta10() -> [ReferenceBranch; 3] {
    [
        ReferenceBranch {
            name: "low",
            csv: include_str!("fixtures/zz1962b_eta10_branch1.csv"),
        },
        ReferenceBranch {
            name: "middle",
            csv: include_str!("fixtures/zz1962b_eta10_branch2.csv"),
        },
        ReferenceBranch {
            name: "high",
            csv: include_str!("fixtures/zz1962b_eta10_branch3.csv"),
        },
    ]
}

fn branch_summary(csv: &str) -> (f64, f64, f64, f64, f64) {
    let rows = parse_reference(csv);
    let inlet_slope = rows[0][2];
    let last = rows.last().expect("nonempty scalar reference");
    let peak = rows
        .iter()
        .max_by(|left, right| left[6].total_cmp(&right[6]))
        .expect("reaction-rate values present");
    (
        last[4],
        last[3],
        (DIFFUSION / ((MASS_FLUX / RHO).powi(2) * TAU)) * inlet_slope.abs(),
        peak[6],
        peak[0],
    )
}

fn parse_reference(csv: &str) -> Vec<[f64; 7]> {
    csv.lines()
        .skip(1)
        .map(|line| {
            let values: Vec<f64> = line
                .split(',')
                .map(|field| field.parse().expect("numeric ZZ-1962-B reference"))
                .collect();
            values.try_into().expect("seven reference columns")
        })
        .collect()
}

fn interpolate(rows: &[[f64; 7]], z: f64) -> [f64; 7] {
    let pos = z.clamp(0.0, 1.0) * (rows.len() - 1) as f64;
    let lo = pos.floor() as usize;
    let hi = (lo + 1).min(rows.len() - 1);
    let f = pos - lo as f64;
    std::array::from_fn(|i| rows[lo][i] + f * (rows[hi][i] - rows[lo][i]))
}

fn compare_to_branch(
    reactor: &SimpleReactorTask,
    branch: ReferenceBranch,
) -> (f64, f64, f64, f64, f64) {
    let rows = parse_reference(branch.csv);
    let profile = reactor.solver.get_solution().expect("published profile");
    let mesh = reactor.solver.x_mesh.as_ref().expect("published mesh");
    let mut max_t = 0.0_f64;
    let mut max_y = 0.0_f64;
    let mut max_flux = 0.0_f64;
    let mut max_rate = 0.0_f64;
    let mut max_rate_eta = 0.0_f64;
    let mut peak_rate_value = f64::NEG_INFINITY;
    for (node, x) in mesh.iter().enumerate() {
        let z = *x / reactor.L;
        let reference = interpolate(&rows, z);
        max_t = max_t.max((profile[(node, 0)] - reference[3]).abs());
        max_y = max_y.max((profile[(node, 2)] - reference[4]).abs());
        max_y = max_y.max((profile[(node, 4)] - reference[5]).abs());
        let expected_flux_a = -RHO * DIFFUSION * reference[2] / (MASS_FLUX / RHO * TAU);
        max_flux = max_flux.max((profile[(node, 3)] - expected_flux_a).abs());
        let rate = profile[(node, 2)] * (-ACTIVATION_ENERGY / (R * profile[(node, 0)])).exp();
        max_rate = max_rate.max((rate - reference[6]).abs());
        if rate > peak_rate_value {
            peak_rate_value = rate;
            max_rate_eta = *x / (MASS_FLUX / RHO * TAU);
        }
    }
    let expected_peak_eta = rows
        .iter()
        .max_by(|left, right| left[6].total_cmp(&right[6]))
        .map(|row| row[0])
        .unwrap_or_default();
    (
        max_t,
        max_y,
        max_flux,
        max_rate,
        (max_rate_eta - expected_peak_eta).abs(),
    )
}

fn assert_high_branch_solution(reactor: &SimpleReactorTask, reference: ReferenceBranch) {
    assert_high_branch_solution_with_fraction_tolerance(reactor, reference, 2e-7);
}

fn assert_high_branch_solution_with_fraction_tolerance(
    reactor: &SimpleReactorTask,
    reference: ReferenceBranch,
    fraction_tolerance: f64,
) {
    let profile = reactor.solver.get_solution().expect("published profile");
    let mesh = reactor.solver.x_mesh.as_ref().expect("published mesh");
    let (dt, dy, dj, drate, peak_eta_error) = compare_to_branch(reactor, reference);
    let dq = max_heat_flux_error_to_branch(reactor, reference);
    let (peak_rate, peak_rate_x) = physical_rate_diagnostics(reactor);
    let (ode_scaled, ode_absolute, ode_row, ode_z) = production_ode_residual_diagnostic(reactor);
    println!(
        "ZZ-1962-B high profile diagnostic L={:.3}m nodes={} errors[T,Y,J,q,rate]=[{dt:.3e},{dy:.3e},{dj:.3e},{dq:.3e},{drate:.3e}] peak_eta_error={peak_eta_error:.3e} physical_rate_peak={peak_rate:.6e}mol/(m3 s) at x={peak_rate_x:.6e}m production_ODE_scaled_defect={ode_scaled:.3e} abs={ode_absolute:.3e} row={ode_row} z={ode_z:.6}",
        reactor.L,
        profile.nrows(),
    );
    let mut first_integral_error = 0.0_f64;
    let mut mass_fraction_sum_error = 0.0_f64;
    let mut flux_sum_error = 0.0_f64;
    let mut energy_min = f64::INFINITY;
    let mut energy_max = f64::NEG_INFINITY;
    let mut energy_scale = 0.0_f64;
    for node in 0..profile.nrows() {
        let temperature = profile[(node, 0)];
        let ya = profile[(node, 2)];
        let yb = profile[(node, 4)];
        assert!(temperature.is_finite() && temperature > 0.0);
        assert!(ya.is_finite() && (-1e-7..=1.0 + 1e-7).contains(&ya));
        assert!(yb.is_finite() && (-1e-7..=1.0 + 1e-7).contains(&yb));
        first_integral_error =
            first_integral_error.max((temperature - T0_K - DELTA_T_SCALE * yb).abs());
        mass_fraction_sum_error = mass_fraction_sum_error.max((ya + yb - 1.0).abs());
        flux_sum_error = flux_sum_error.max((profile[(node, 3)] + profile[(node, 5)]).abs());
        let energy = profile[(node, 1)] - MASS_FLUX * CP * temperature
            + Q_MASS * (profile[(node, 3)] - MASS_FLUX * ya);
        energy_min = energy_min.min(energy);
        energy_max = energy_max.max(energy);
        energy_scale = energy_scale.max(
            profile[(node, 1)].abs()
                + (MASS_FLUX * CP * temperature).abs()
                + (Q_MASS * (profile[(node, 3)] - MASS_FLUX * ya)).abs(),
        );
    }
    let energy_relative_drift = (energy_max - energy_min) / energy_scale.max(1e-12);
    assert!(
        first_integral_error < 2e-3,
        "first integral defect={first_integral_error:e} K"
    );
    assert!(mass_fraction_sum_error < 2e-5);
    assert!(flux_sum_error < 2e-5);
    assert!(dt < 2e-3, "high-branch temperature profile error={dt:e} K");
    assert!(
        dy < fraction_tolerance,
        "high-branch mass-fraction profile error={dy:e} at L={}m",
        reactor.L
    );
    assert!(dj < 2e-8, "high-branch flux profile error={dj:e}");
    assert!(drate < 2e-5, "high-branch rate profile error={drate:e}");
    assert!(peak_eta_error < 0.2);
    assert!(dq < 0.05, "high-branch heat-flux profile error={dq:e} W/m2");
    assert!(
        ode_scaled < 1e-4,
        "production ODE scaled defect={ode_scaled:e} row={ode_row} z={ode_z}"
    );
    assert!(energy_relative_drift < 2e-7);

    let first = 0;
    let last = profile.nrows() - 1;
    assert!((profile[(first, 0)] - T0_K).abs() < 2e-5);
    assert!((profile[(first, 2)] - 1.0).abs() < 2e-7);
    assert!(profile[(first, 4)].abs() < 2e-7);
    for flux_state in [1, 3, 5] {
        assert!(profile[(last, flux_state)].abs() < 2e-5);
    }
    assert_eq!(mesh.len(), profile.nrows());
}

fn max_profile_distance_to_high_branch(
    reactor: &SimpleReactorTask,
    branch: ReferenceBranch,
) -> f64 {
    let (dt, dy, dj, drate, _) = compare_to_branch(reactor, branch);
    dt + 3000.0 * dy + 100.0 * dj + 1000.0 * drate
}

fn max_heat_flux_error_to_branch(reactor: &SimpleReactorTask, branch: ReferenceBranch) -> f64 {
    let rows = parse_reference(branch.csv);
    let profile = reactor.solver.get_solution().expect("published profile");
    let mesh = reactor.solver.x_mesh.as_ref().expect("published mesh");
    let velocity = MASS_FLUX / RHO;
    let mut max_error = 0.0_f64;
    for (node, x) in mesh.iter().enumerate() {
        let reference = interpolate(&rows, *x / reactor.L);
        let expected_q = LAMBDA * DELTA_T_SCALE * reference[2] / (velocity * TAU);
        max_error = max_error.max((profile[(node, 1)] - expected_q).abs());
    }
    max_error
}

fn physical_rate_diagnostics(reactor: &SimpleReactorTask) -> (f64, f64) {
    let profile = reactor.solver.get_solution().expect("published profile");
    let mesh = reactor.solver.x_mesh.as_ref().expect("published mesh");
    let mut peak_rate = f64::NEG_INFINITY;
    let mut peak_x = 0.0;
    for (node, x) in mesh.iter().enumerate() {
        let rate = RHO / MOLAR_MASS
            * A0
            * profile[(node, 2)]
            * (-ACTIVATION_ENERGY / (R * profile[(node, 0)])).exp();
        if rate > peak_rate {
            peak_rate = rate;
            peak_x = *x;
        }
    }
    (peak_rate, peak_x)
}

fn nonuniform_derivative(values: &[f64], coordinates: &[f64], index: usize) -> f64 {
    if index == 0 {
        let h0 = coordinates[1] - coordinates[0];
        let h1 = coordinates[2] - coordinates[1];
        return -(2.0 * h0 + h1) / (h0 * (h0 + h1)) * values[0] + (h0 + h1) / (h0 * h1) * values[1]
            - h0 / (h1 * (h0 + h1)) * values[2];
    }
    if index + 1 == values.len() {
        let count = values.len();
        let h0 = coordinates[count - 2] - coordinates[count - 3];
        let h1 = coordinates[count - 1] - coordinates[count - 2];
        return h1 / (h0 * (h0 + h1)) * values[count - 3]
            - (h0 + h1) / (h0 * h1) * values[count - 2]
            + (h0 + 2.0 * h1) / (h1 * (h0 + h1)) * values[count - 1];
    }
    let h0 = coordinates[index] - coordinates[index - 1];
    let h1 = coordinates[index + 1] - coordinates[index];
    -h1 / (h0 * (h0 + h1)) * values[index - 1]
        + (h1 - h0) / (h0 * h1) * values[index]
        + h0 / (h1 * (h0 + h1)) * values[index + 1]
}

fn production_ode_residual_diagnostic(reactor: &SimpleReactorTask) -> (f64, f64, usize, f64) {
    let profile = reactor.solver.get_solution().expect("published profile");
    let mesh = reactor.solver.x_mesh.as_ref().expect("published mesh");
    let coordinates: Vec<f64> = mesh.iter().map(|x| *x / reactor.L).collect();
    let variables: Vec<&str> = reactor.solver.unknowns.iter().map(String::as_str).collect();
    let mut internal_states = vec![vec![0.0; reactor.solver.unknowns.len()]; profile.nrows()];
    for (state_index, name) in reactor.solver.unknowns.iter().enumerate() {
        for node in 0..profile.nrows() {
            let physical = profile[(node, state_index)];
            internal_states[node][state_index] = match name.as_str() {
                "Teta" => (physical - T0_K) / DELTA_T_SCALE,
                "q" => physical * reactor.L / DELTA_T_SCALE,
                name if name.starts_with('J') => physical * reactor.L,
                _ => physical,
            };
        }
    }
    let state_profiles: Vec<Vec<f64>> = (0..reactor.solver.unknowns.len())
        .map(|column| internal_states.iter().map(|state| state[column]).collect())
        .collect();

    let mut max_scaled_defect = 0.0_f64;
    let mut max_absolute_defect = 0.0_f64;
    let mut max_row = 0usize;
    let mut max_z = 0.0_f64;
    for (node, state) in internal_states.iter().enumerate() {
        for (row, equation) in reactor.solver.eq_system.iter().enumerate() {
            let state_derivative: Vec<f64> = state_profiles
                .iter()
                .map(|values| nonuniform_derivative(values, &coordinates, node))
                .collect();
            let rhs = equation.eval_expression(&variables, state);
            let absolute = (state_derivative[row] - rhs).abs();
            let scaled = absolute / (1.0 + rhs.abs());
            if scaled > max_scaled_defect {
                max_scaled_defect = scaled;
                max_absolute_defect = absolute;
                max_row = row;
                max_z = coordinates[node];
            }
        }
    }
    (max_scaled_defect, max_absolute_defect, max_row, max_z)
}

fn reaction_zone_width_10_percent(reactor: &SimpleReactorTask) -> f64 {
    let solution = reactor.solver.get_solution().expect("published profile");
    let mesh = reactor.solver.x_mesh.as_ref().expect("published mesh");
    let rates: Vec<f64> = (0..solution.nrows())
        .map(|node| solution[(node, 2)] * (-ACTIVATION_ENERGY / (R * solution[(node, 0)])).exp())
        .collect();
    let peak = rates.iter().copied().fold(0.0_f64, f64::max);
    let active_x: Vec<f64> = rates
        .iter()
        .zip(mesh.iter())
        .filter_map(|(rate, x)| (*rate >= 0.1 * peak).then_some(*x))
        .collect();
    match (active_x.first(), active_x.last()) {
        (Some(first), Some(last)) => last - first,
        _ => 0.0,
    }
}

fn report_coarse_route_branch_selection(
    reactor: &SimpleReactorTask,
    route_name: &str,
) -> &'static str {
    let profile = reactor.solver.get_solution().expect("published profile");
    let last = profile.nrows() - 1;
    let mut first_integral_error = 0.0_f64;
    let mut mass_fraction_sum_error = 0.0_f64;
    let mut flux_sum_error = 0.0_f64;
    let mut energy_min = f64::INFINITY;
    let mut energy_max = f64::NEG_INFINITY;
    let mut energy_scale = 0.0_f64;
    for node in 0..profile.nrows() {
        let temperature = profile[(node, 0)];
        let ya = profile[(node, 2)];
        let yb = profile[(node, 4)];
        assert!(temperature.is_finite() && temperature > 0.0);
        assert!((-1e-7..=1.0 + 1e-7).contains(&ya));
        assert!((-1e-7..=1.0 + 1e-7).contains(&yb));
        first_integral_error =
            first_integral_error.max((temperature - T0_K - DELTA_T_SCALE * yb).abs());
        mass_fraction_sum_error = mass_fraction_sum_error.max((ya + yb - 1.0).abs());
        flux_sum_error = flux_sum_error.max((profile[(node, 3)] + profile[(node, 5)]).abs());
        let energy = profile[(node, 1)] - MASS_FLUX * CP * temperature
            + Q_MASS * (profile[(node, 3)] - MASS_FLUX * ya);
        energy_min = energy_min.min(energy);
        energy_max = energy_max.max(energy);
        energy_scale = energy_scale.max(
            profile[(node, 1)].abs()
                + (MASS_FLUX * CP * temperature).abs()
                + (Q_MASS * (profile[(node, 3)] - MASS_FLUX * ya)).abs(),
        );
    }
    let energy_relative_drift = (energy_max - energy_min) / energy_scale.max(1e-12);
    assert!(profile.iter().all(|value| value.is_finite()));
    assert!(first_integral_error < 2e-3);
    assert!(mass_fraction_sum_error < 2e-5);
    assert!(flux_sum_error < 2e-5);
    assert!(energy_relative_drift < 2e-7);
    assert!(profile[(last, 0)] > 0.0);
    assert!((-1e-7..=1.0 + 1e-7).contains(&profile[(last, 2)]));
    assert!((-1e-7..=1.0 + 1e-7).contains(&profile[(last, 4)]));

    let mut candidates: Vec<_> = branches_eta10()
        .iter()
        .map(|branch| {
            let (dt, dy, dj, drate, peak_eta_error) = compare_to_branch(reactor, *branch);
            let score = dt + 3000.0 * dy + 100.0 * dj + 1000.0 * drate;
            (*branch, score, (dt, dy, dj, drate, peak_eta_error))
        })
        .collect();
    candidates.sort_by(|left, right| left.1.total_cmp(&right.1));
    let (nearest, _, (dt, dy, dj, drate, peak_eta_error)) = candidates[0];
    let selected_branch = nearest.name;
    let epsilon_inlet = profile[(0, 3)].abs() / MASS_FLUX;
    let dq = max_heat_flux_error_to_branch(reactor, nearest);
    let (peak_rate, peak_rate_x) = physical_rate_diagnostics(reactor);
    println!(
        "ZZ-1962-B branch selection route={} nodes={} nearest={} YAout={:.9} Tout={:.6}K epsilon_inlet={epsilon_inlet:.7e} nearest_errors[T,Y,J,q,rate]=[{dt:.3e},{dy:.3e},{dj:.3e},{dq:.3e},{drate:.3e}] peak_eta_error={peak_eta_error:.3e} invariant[T,sumY,sumJ,energy]=[{first_integral_error:.3e},{mass_fraction_sum_error:.3e},{flux_sum_error:.3e},{energy_relative_drift:.3e}] physical_rate_peak={peak_rate:.6e}mol/(m3 s) at x={peak_rate_x:.6e}m zone_width_10pct_x={:.6e}m",
        route_name,
        profile.nrows(),
        nearest.name,
        profile[(last, 2)],
        profile[(last, 0)],
        reaction_zone_width_10_percent(reactor),
    );
    for (branch, score, (dt, dy, dj, drate, peak_eta_error)) in candidates {
        println!(
            "ZZ-1962-B branch candidate route={} branch={} score={score:.6e} errors[T,Y,J,q,rate]=[{dt:.3e},{dy:.3e},{dj:.3e},{:.3e},{drate:.3e}] peak_eta_error={peak_eta_error:.3e}",
            route_name,
            branch.name,
            max_heat_flux_error_to_branch(reactor, branch),
        );
    }
    selected_branch
}

#[test]
fn zz1962b_parameters_map_to_the_claimed_dimensionless_groups() {
    let u = MASS_FLUX / RHO;
    let alpha = DIFFUSION / (u * u * TAU);
    let beta = ACTIVATION_ENERGY / (R * T0_K);
    let delta = Q_MASS / (CP * T0_K);
    let lewis = LAMBDA / (RHO * DIFFUSION * CP);
    let eta_max = 0.05 / (u * TAU);
    assert!((alpha - 2.0).abs() < 1e-14);
    assert!((beta - 7.0).abs() < 1e-14);
    assert!((delta - 10.0).abs() < 1e-14);
    assert!((lewis - 1.0).abs() < 1e-14);
    assert!((eta_max - 10.0).abs() < 1e-14);
    assert!((Q_MOLAR / MOLAR_MASS - Q_MASS).abs() < 1e-8);
    println!("ZZ-1962-B alpha={alpha} beta={beta} delta={delta} Le={lewis} eta_max={eta_max}");
}

#[test]
fn zz1962b_reduced_reference_records_multistart_branches_and_finite_domain_effect() {
    let branches = branches_eta10().map(|branch| branch_summary(branch.csv));
    assert!((branches[0].0 - 0.9881569499).abs() < 2e-8);
    assert!((branches[1].0 - 0.8509502583).abs() < 2e-8);
    assert!((branches[2].0 - 0.1479815712).abs() < 2e-8);
    assert!(branches[0].0 > branches[1].0 && branches[1].0 > branches[2].0);
    assert!((branches[2].4 - 4.555).abs() < 0.01);

    const ETA15_HIGH: &str = include_str!("fixtures/zz1962b_eta15_branch2.csv");
    const ETA20_HIGH: &str = include_str!("fixtures/zz1962b_eta20_branch1.csv");
    let high_eta10 = branches[2];
    let high_eta15 = branch_summary(ETA15_HIGH);
    let high_eta20 = branch_summary(ETA20_HIGH);
    assert!(high_eta10.0 > high_eta15.0 && high_eta15.0 > high_eta20.0);
    assert!(high_eta10.1 < high_eta15.1 && high_eta15.1 < high_eta20.1);
    assert!((high_eta10.2 - high_eta15.2).abs() < 4e-4);
    assert!((high_eta15.2 - high_eta20.2).abs() < 2e-5);
    assert!((high_eta10.3 - high_eta20.3).abs() < 2e-7);
    assert!((high_eta10.4 - high_eta20.4).abs() < 0.02);
    println!(
        "ZZ-1962-B reduced branches eta=10: low/middle/high Y_A,out={:.8}/{:.8}/{:.8}; high-branch etaMax 10/15/20 Y_A,out={:.8}/{:.8}/{:.8}, T_out={:.3}/{:.3}/{:.3} K, epsilon_in={:.7e}/{:.7e}/{:.7e}",
        branches[0].0,
        branches[1].0,
        branches[2].0,
        high_eta10.0,
        high_eta15.0,
        high_eta20.0,
        high_eta10.1,
        high_eta15.1,
        high_eta20.1,
        high_eta10.2,
        high_eta15.2,
        high_eta20.2,
    );
}

#[test]
fn zz1962b_full_bvp_routes_match_the_classified_reduced_branch() {
    let mut matched = Vec::<(&'static str, f64)>::new();
    for route in [Route::Sci, Route::Damped, Route::Frozen] {
        let seed_branch = branches_eta10()[0];
        let reactor = match solve_full(0.05, seed_branch, route) {
            Ok(reactor) => reactor,
            Err(error) => {
                panic!(
                    "ZZ-1962-B branch {} did not solve: {error}",
                    seed_branch.name
                );
            }
        };
        let profile = reactor.solver.get_solution().unwrap();
        let mesh = reactor.solver.x_mesh.as_ref().unwrap();
        let mut invariant_error = 0.0_f64;
        let mut energy_invariant_min = f64::INFINITY;
        let mut energy_invariant_max = f64::NEG_INFINITY;
        let mut energy_invariant_scale = 0.0_f64;
        let mut species_sum_error = 0.0_f64;
        let mut flux_sum_error = 0.0_f64;
        for node in 0..profile.nrows() {
            let t = profile[(node, 0)];
            let ya = profile[(node, 2)];
            let yb = profile[(node, 4)];
            invariant_error = invariant_error.max((t - T0_K - DELTA_T_SCALE * yb).abs());
            species_sum_error = species_sum_error.max((ya + yb - 1.0).abs());
            flux_sum_error = flux_sum_error.max((profile[(node, 3)] + profile[(node, 5)]).abs());
            let energy = profile[(node, 1)] - MASS_FLUX * CP * t
                + Q_MASS * (profile[(node, 3)] - MASS_FLUX * ya);
            energy_invariant_min = energy_invariant_min.min(energy);
            energy_invariant_max = energy_invariant_max.max(energy);
            energy_invariant_scale = energy_invariant_scale.max(
                profile[(node, 1)].abs()
                    + (MASS_FLUX * CP * t).abs()
                    + (Q_MASS * (profile[(node, 3)] - MASS_FLUX * ya)).abs(),
            );
        }
        let energy_invariant_relative_drift =
            (energy_invariant_max - energy_invariant_min) / energy_invariant_scale.max(1e-12);
        assert!(
            invariant_error < 2e-3,
            "first-integral defect={invariant_error:e}"
        );
        assert!(
            species_sum_error < 2e-5,
            "mass-fraction sum defect={species_sum_error:e}"
        );
        assert!(
            flux_sum_error < 2e-5,
            "species-flux sum defect={flux_sum_error:e}"
        );
        assert!(energy_invariant_relative_drift < 2e-7);

        let mut errors: Vec<_> = branches_eta10()
            .iter()
            .map(|branch| (*branch, compare_to_branch(&reactor, *branch)))
            .collect();
        errors.sort_by(|left, right| {
            (left.1.0 + left.1.1 * 3000.0 + left.1.2 * 100.0 + left.1.3 * 1000.0).total_cmp(
                &(right.1.0 + right.1.1 * 3000.0 + right.1.2 * 100.0 + right.1.3 * 1000.0),
            )
        });
        let (branch, (dt, dy, dj, drate, rate_peak_eta_error)) = errors[0];
        let outlet = profile.nrows() - 1;
        let inlet_flux = profile[(0, 3)];
        let epsilon_inlet = RHO * DIFFUSION / MASS_FLUX * (inlet_flux / (RHO * DIFFUSION)).abs();
        println!(
            "ZZ-1962-B full route={route:?} seed_branch={} matched_branch={} nodes={} Y_A(L)={:.8} T(L)={:.6}K epsilon_inlet={epsilon_inlet:.6e} profile_err[T,Y,J,rate]=[{dt:.3e},{dy:.3e},{dj:.3e},{drate:.3e}] rate_peak_eta_error={rate_peak_eta_error:.3e} first_integral={invariant_error:.3e} energy_invariant={energy_invariant_relative_drift:.3e} mesh_eta_max={:.1}",
            seed_branch.name,
            branch.name,
            profile.nrows(),
            profile[(outlet, 2)],
            profile[(outlet, 0)],
            mesh[mesh.len() - 1] / (MASS_FLUX / RHO * TAU),
        );
        let (dt_tolerance, dy_tolerance, dj_tolerance) = match route {
            Route::Sci => (2e-3, 2e-7, 2e-8),
            Route::Damped => (8e-2, 2e-5, 2e-6),
            Route::Frozen => (4e-2, 1e-5, 1e-6),
        };
        assert!(dt < dt_tolerance, "temperature error={dt:e} K");
        assert!(dy < dy_tolerance, "mass-fraction error={dy:e}");
        assert!(dj < dj_tolerance, "species-flux error={dj:e}");
        assert!(drate < 2e-5, "dimensionless reaction-rate error={drate:e}");
        assert!(rate_peak_eta_error < 0.2);
        matched.push((branch.name, profile[(outlet, 2)]));
    }
    assert!(
        !matched.is_empty(),
        "at least one seeded production BVP should solve"
    );
    let distinct: std::collections::BTreeSet<_> = matched.iter().map(|(name, _)| *name).collect();
    println!("ZZ-1962-B discovered KiThe branches: {distinct:?}");
}

#[test]
fn zz1962b_correction_newton_matches_low_branch_and_retains_telemetry() {
    let branch = branches_eta10()[0];
    let coarse = solve_correction_newton(0.05, branch, 129)
        .unwrap_or_else(|error| panic!("CorrectionNewton coarse low-branch solve: {error}"));
    let coarse_errors = compare_to_branch(&coarse, branch);
    let fine = solve_correction_newton(0.05, branch, 257)
        .unwrap_or_else(|error| panic!("CorrectionNewton fine low-branch solve: {error}"));
    let fine_errors = compare_to_branch(&fine, branch);
    let fine_heat_error = max_heat_flux_error_to_branch(&fine, branch);
    println!(
        "ZZ-1962-B CorrectionNewton fixed mesh low branch coarse=129 [dT,dY,dJ,drate]=[{:.3e},{:.3e},{:.3e},{:.3e}] fine=257 [{:.3e},{:.3e},{:.3e},{:.3e}] dq={fine_heat_error:.3e}",
        coarse_errors.0,
        coarse_errors.1,
        coarse_errors.2,
        coarse_errors.3,
        fine_errors.0,
        fine_errors.1,
        fine_errors.2,
        fine_errors.3,
    );
    assert!(
        fine_errors.0 < coarse_errors.0
            && fine_errors.1 < coarse_errors.1
            && fine_errors.2 < coarse_errors.2
            && fine_errors.3 < coarse_errors.3,
        "CorrectionNewton mesh refinement should improve all low-branch profile errors: coarse={coarse_errors:?}, fine={fine_errors:?}"
    );
    assert!(fine_errors.0 < 8e-2);
    assert!(fine_errors.1 < 2e-5);
    assert!(fine_errors.2 < 2e-6);
    assert!(fine_errors.3 < 2e-5);
    assert!(fine_errors.4 < 0.2);
    assert!(fine_heat_error < 0.05);

    let profile = fine.solver.get_solution().expect("published CN profile");
    let mut first_integral_error = 0.0_f64;
    let mut mass_sum_error = 0.0_f64;
    let mut flux_sum_error = 0.0_f64;
    let mut energy_min = f64::INFINITY;
    let mut energy_max = f64::NEG_INFINITY;
    let mut energy_scale = 0.0_f64;
    for node in 0..profile.nrows() {
        let temperature = profile[(node, 0)];
        let ya = profile[(node, 2)];
        let yb = profile[(node, 4)];
        assert!(temperature.is_finite() && temperature > 0.0);
        assert!(ya.is_finite() && (-1e-7..=1.0 + 1e-7).contains(&ya));
        assert!(yb.is_finite() && (-1e-7..=1.0 + 1e-7).contains(&yb));
        first_integral_error =
            first_integral_error.max((temperature - T0_K - DELTA_T_SCALE * yb).abs());
        mass_sum_error = mass_sum_error.max((ya + yb - 1.0).abs());
        flux_sum_error = flux_sum_error.max((profile[(node, 3)] + profile[(node, 5)]).abs());
        let energy = profile[(node, 1)] - MASS_FLUX * CP * temperature
            + Q_MASS * (profile[(node, 3)] - MASS_FLUX * ya);
        energy_min = energy_min.min(energy);
        energy_max = energy_max.max(energy);
        energy_scale = energy_scale.max(
            profile[(node, 1)].abs()
                + (MASS_FLUX * CP * temperature).abs()
                + (Q_MASS * (profile[(node, 3)] - MASS_FLUX * ya)).abs(),
        );
    }
    let energy_relative_drift = (energy_max - energy_min) / energy_scale.max(1e-12);
    assert!(first_integral_error < 2e-3);
    assert!(mass_sum_error < 2e-5);
    assert!(flux_sum_error < 2e-5);
    assert!(energy_relative_drift < 2e-7);

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
fn zz1962b_damped_grid_parameter_sweep_measures_branch_accuracy() {
    let branch = branches_eta10()[0];
    let length_m = 0.05;
    let seed_nodes = 65;
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
    let seed_guess = reference_initial_guess(length_m, branch, seed_nodes);
    let baseline = solve_full_with_grid(
        length_m,
        seed_guess.clone(),
        Route::Damped,
        branch.name,
        20_000,
        1e-9,
        Some((GridRefinementMethod::DoublePoints, 0)),
    )
    .expect("unrefined Damped branch seed should solve");
    let baseline_error = compare_to_branch(&baseline, branch);
    let mut measured = Vec::new();

    for (name, method) in candidates {
        let reactor = solve_full_with_grid(
            length_m,
            seed_guess.clone(),
            Route::Damped,
            branch.name,
            20_000,
            1e-9,
            Some((method, 3)),
        )
        .unwrap_or_else(|error| panic!("ZZ-1962-B {name} failed: {error}"));
        let errors = compare_to_branch(&reactor, branch);
        let nodes = reactor.solver.x_mesh.as_ref().unwrap().len();
        println!(
            "ZZ-1962-B tuned grid={name} nodes={nodes} dT={:.3e} dY={:.3e} dJ={:.3e} drate={:.3e} peak_eta={:.3e}",
            errors.0, errors.1, errors.2, errors.3, errors.4,
        );
        measured.push((name, nodes, errors));
    }
    println!(
        "ZZ-1962-B fixed seed nodes={seed_nodes} dT={:.3e} dY={:.3e} dJ={:.3e} drate={:.3e} peak_eta={:.3e}",
        baseline_error.0, baseline_error.1, baseline_error.2, baseline_error.3, baseline_error.4,
    );
    let grcar = measured
        .iter()
        .find(|row| row.0 == "Grcar-1e-2-5e-2")
        .expect("the Grcar branch candidate must be measured");
    assert_eq!(grcar.1, 207);
    assert!(
        grcar.2.0 < baseline_error.0
            && grcar.2.1 < baseline_error.1
            && grcar.2.2 < baseline_error.2
            && grcar.2.3 < baseline_error.3,
        "selected Grcar settings must improve all branch-profile fields: baseline={baseline_error:?}, candidate={grcar:?}"
    );
    assert!(
        measured.iter().any(|(_, _, errors)| {
            errors.0 < baseline_error.0
                && errors.1 < baseline_error.1
                && errors.2 < baseline_error.2
                && errors.3 < baseline_error.3
        }),
        "no adaptive strategy improved all branch-profile errors: baseline={baseline_error:?}, candidates={measured:?}"
    );
}

#[test]
#[ignore = "release-only high-branch validation; expensive BVP_sci solve"]
fn zz1962b_full_bvp_sci_reproduces_high_temperature_branch_eta10() {
    let high = branches_eta10()[2];
    let guess = reference_initial_guess(0.05, high, 301);
    let reactor = solve_full_with_limits(0.05, guess, Route::Sci, high.name, 1200, 1e-9)
        .unwrap_or_else(|error| panic!("high branch BVP_sci failed: {error}"));
    assert_high_branch_solution(&reactor, high);

    let profile = reactor.solver.get_solution().unwrap();
    let last = profile.nrows() - 1;
    let epsilon_inlet = profile[(0, 3)].abs() / MASS_FLUX;
    let zone_width = reaction_zone_width_10_percent(&reactor);
    assert!((profile[(last, 2)] - 0.1479815712).abs() < 2e-7);
    assert!((profile[(last, 0)] - 2856.055286).abs() < 2e-3);
    assert!((epsilon_inlet - 0.1037427).abs() < 2e-5);
    assert!(zone_width.is_finite() && zone_width > 0.0);
    println!(
        "ZZ-1962-B high branch BVP_sci etaMax=10 nodes={} YAout={:.9} Tout={:.6}K epsilon_inlet={:.7e} reaction_zone_10pct_width_x={:.6e}m",
        profile.nrows(),
        profile[(last, 2)],
        profile[(last, 0)],
        epsilon_inlet,
        zone_width,
    );
}

#[test]
#[ignore = "release-only branch-selection diagnostic; adaptive Damped solve"]
fn zz1962b_full_bvp_damped_reports_branch_selection_eta10() {
    let high = branches_eta10()[2];
    let reactor = solve_full(0.05, high, Route::Damped)
        .unwrap_or_else(|error| panic!("Damped high branch failed: {error}"));
    report_coarse_route_branch_selection(&reactor, "Damped");
}

#[test]
#[ignore = "release-only branch-selection diagnostic; dense Frozen solve"]
fn zz1962b_full_bvp_frozen_reports_branch_selection_eta10() {
    let high = branches_eta10()[2];
    let reactor = solve_full(0.05, high, Route::Frozen)
        .unwrap_or_else(|error| panic!("Frozen high branch failed: {error}"));
    report_coarse_route_branch_selection(&reactor, "Frozen");
}

#[test]
#[ignore = "release-only branch-selection diagnostic; smooth non-reference Damped seed"]
fn zz1962b_full_bvp_damped_reports_selection_from_smooth_seed_eta10() {
    let guess = smooth_high_branch_guess(10.0, 256, 0.85, 0.2);
    let reactor = solve_full_with_guess(0.05, guess, Route::Damped, "smooth theta_R=0.85 a=0.2")
        .unwrap_or_else(|error| panic!("Damped smooth-seed solve failed: {error}"));
    let selected = report_coarse_route_branch_selection(&reactor, "Damped/smooth-seed");
    assert_eq!(selected, "high", "Damped branch selection changed");
}

#[test]
#[ignore = "release-only branch-selection diagnostic; smooth non-reference Frozen seed"]
fn zz1962b_full_bvp_frozen_reports_selection_from_smooth_seed_eta10() {
    let guess = smooth_high_branch_guess(10.0, 129, 0.85, 0.2);
    let reactor = solve_full_with_guess(0.05, guess, Route::Frozen, "smooth theta_R=0.85 a=0.2")
        .unwrap_or_else(|error| panic!("Frozen smooth-seed solve failed: {error}"));
    let selected = report_coarse_route_branch_selection(&reactor, "Frozen/smooth-seed");
    assert_eq!(selected, "high", "Frozen branch selection changed");
}

#[test]
#[ignore = "release-only high-branch validation; expensive BVP_sci solve"]
fn zz1962b_full_bvp_sci_reaches_high_branch_from_smooth_nonreference_seed() {
    let high = branches_eta10()[2];
    let guess = smooth_high_branch_guess(10.0, 301, 0.85, 0.2);
    let reactor = solve_full_with_limits(
        0.05,
        guess,
        Route::Sci,
        "smooth theta_R=0.85 a=0.2",
        1200,
        1e-9,
    )
    .unwrap_or_else(|error| panic!("smooth-seed high-branch BVP_sci failed: {error}"));
    assert_high_branch_solution(&reactor, high);
    println!(
        "ZZ-1962-B nonreference seed converged to high branch; zone_width_10pct_x={:.6e}m",
        reaction_zone_width_10_percent(&reactor),
    );
}

#[test]
#[ignore = "release-only domain-continuation validation; multiple expensive BVP_sci solves"]
fn zz1962b_high_branch_bvp_sci_continues_domain_10_to_15_to_20() {
    let high10 = branches_eta10()[2];
    let high15 = ReferenceBranch {
        name: "high-eta15",
        csv: include_str!("fixtures/zz1962b_eta15_branch2.csv"),
    };
    let high20 = ReferenceBranch {
        name: "high-eta20",
        csv: include_str!("fixtures/zz1962b_eta20_branch1.csv"),
    };

    let guess10 = reference_initial_guess(0.05, high10, 301);
    let solution10 = solve_full_with_limits(0.05, guess10, Route::Sci, high10.name, 1200, 1e-9)
        .unwrap_or_else(|error| panic!("etaMax=10 high branch failed: {error}"));
    assert_high_branch_solution(&solution10, high10);
    let guess15 = prolong_published_solution(&solution10, 10.0, 15.0, 151);
    let solution15 = solve_full_with_limits(
        0.075,
        guess15,
        Route::Sci,
        "etaMax=10 continuation",
        3000,
        1e-9,
    )
    .unwrap_or_else(|error| panic!("etaMax=15 high-branch continuation failed: {error}"));
    assert_high_branch_solution_with_fraction_tolerance(&solution15, high15, 1e-6);

    let guess20 = prolong_published_solution(&solution15, 15.0, 20.0, 201);
    let solution20 = solve_full_with_limits(
        0.1,
        guess20,
        Route::Sci,
        "etaMax=15 continuation",
        3000,
        1e-9,
    )
    .unwrap_or_else(|error| panic!("etaMax=20 high-branch continuation failed: {error}"));
    assert_high_branch_solution_with_fraction_tolerance(&solution20, high20, 1e-6);

    for (eta_max, reactor, reference) in [
        (10.0, &solution10, high10),
        (15.0, &solution15, high15),
        (20.0, &solution20, high20),
    ] {
        let profile = reactor.solver.get_solution().unwrap();
        let last = profile.nrows() - 1;
        let epsilon_inlet = profile[(0, 3)].abs() / MASS_FLUX;
        println!(
            "ZZ-1962-B domain continuation etaMax={eta_max:.0} L={:.3}m nodes={} YAout={:.9} Tout={:.6}K epsilon_inlet={epsilon_inlet:.7e} high_reference_distance={:.3e} reaction_zone_10pct_width_x={:.6e}m",
            reactor.L,
            profile.nrows(),
            profile[(last, 2)],
            profile[(last, 0)],
            max_profile_distance_to_high_branch(reactor, reference),
            reaction_zone_width_10_percent(reactor),
        );
    }
}
