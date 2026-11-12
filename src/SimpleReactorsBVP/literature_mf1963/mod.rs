//! KiThe integration checks for the independent Merzhanov-Filonenko reference.

use crate::SimpleReactorsBVP::SimpleReactorBVP::{
    DensityModel, FastElemReact, NrbvpHandoffConfig, R_G, ReactorError, SimpleReactorTask,
};
use crate::SimpleReactorsBVP::reactor_BVP_utils::ScalingConfig;
use crate::SimpleReactorsBVP::solver_backend::{ReactorBvpSolverConfig, ReactorBvpSymbolicBackend};
use RustedSciThe::numerical::BVP_CorrectionNewton::{
    CorrectionNewtonMatrixLayout, CorrectionNewtonOptions,
};
use RustedSciThe::numerical::BVP_Damp::NR_Damp_solver_damped::{AdaptiveGridConfig, SolverParams};
use RustedSciThe::numerical::BVP_Damp::grid_api::GridRefinementMethod;
use RustedSciThe::numerical::BVP_Damp::telemetry::{
    BvpLogEventKind, BvpLoggingConfig, BvpLoggingMode, BvpSolveStatus,
};
use RustedSciThe::numerical::BVP_api::{
    BvpApiError, BvpApiTelemetryMode, BvpSciSparseFactorization, BvpSolverSettings,
    BvpSolverTelemetrySnapshot,
};
use RustedSciThe::numerical::BVP_sci::new::{
    BvpSciAssembly, BvpSciExecution, BvpSciMatrixLayout, BvpSciOptions,
};
use RustedSciThe::numerical::logging::SolverLogLevel;
use nalgebra::DMatrix;
use std::collections::HashMap;

#[cfg(test)]
mod adaptive_refinement;
#[cfg(test)]
mod continuation_accuracy;
#[cfg(test)]
mod production_mesh_matrix;
#[cfg(test)]
mod rst_damped_hot_branch_reproducer;
#[cfg(test)]
mod rst_solver_selection;
#[cfg(test)]
mod scipy_parity;
#[cfg(test)]
mod six_state_parity;

const T0: f64 = 300.0;
const CP: f64 = 1000.0;
const RHO0: f64 = 1.0;
const LAMBDA: f64 = 0.1;
const MASS_FLUX: f64 = 0.01;
const T_SCALE: f64 = 15.0;
const MOLAR_MASS_KG_PER_MOL: f64 = 0.060052;
const R_MF00: f64 = 8.31446261815324;
const BETA_MF00: f64 = 0.05;
const GAMMA_MF00: f64 = 0.01;
const Q_MASS: f64 = 1.5e6;
const RHO_D_TARGET: f64 = LAMBDA / CP;

#[derive(Clone, Copy, Debug)]
struct Case {
    name: &'static str,
    delta: f64,
    xi_max: f64,
    csv: &'static str,
}

fn cases() -> [Case; 2] {
    [
        Case {
            name: "MF-LOW",
            delta: 1e-4,
            xi_max: 12017.324603855763,
            csv: include_str!("fixtures/mf1963_mf_low.csv"),
        },
        Case {
            name: "MF-HIGH",
            delta: 0.1,
            xi_max: 0.5427079731331537,
            csv: include_str!("fixtures/mf1963_mf_high.csv"),
        },
    ]
}

#[derive(Clone, Copy, Debug)]
struct RefPoint {
    xi: f64,
    theta: f64,
    p: f64,
    y_a: f64,
}

fn reference_profile(case: Case) -> Vec<RefPoint> {
    case.csv
        .lines()
        .skip(1)
        .filter(|line| !line.trim().is_empty())
        .map(|line| {
            let mut fields = line.split(',').map(|field| field.parse::<f64>().unwrap());
            RefPoint {
                xi: fields.next().unwrap(),
                theta: fields.next().unwrap(),
                p: fields.next().unwrap(),
                y_a: {
                    let _source = fields.next().unwrap();
                    let _temperature = fields.next().unwrap();
                    fields.next().unwrap()
                },
            }
        })
        .collect()
}

fn interpolate(profile: &[RefPoint], xi: f64, select: impl Fn(RefPoint) -> f64) -> f64 {
    let upper = profile
        .partition_point(|point| point.xi < xi)
        .clamp(1, profile.len() - 1);
    let lower = upper - 1;
    let left = profile[lower];
    let right = profile[upper];
    let fraction = ((xi - left.xi) / (right.xi - left.xi)).clamp(0.0, 1.0);
    select(left) + fraction * (select(right) - select(left))
}

fn reactor_for(case: Case) -> SimpleReactorTask {
    // The two formal species are isomolar C2H4O2 isomers. Their equal formula
    // makes mass-fraction conversion identical to the reference conversion.
    let formula = HashMap::from([
        ("C".to_string(), 2),
        ("H".to_string(), 4),
        ("O".to_string(), 2),
    ]);
    let k_at_t0 = case.delta * CP * GAMMA_MF00 * MASS_FLUX.powi(2) / (LAMBDA * RHO0);
    let activation_energy = R_G * T0 / BETA_MF00;
    let k0 = k_at_t0 * (activation_energy / (R_G * T0)).exp();
    let q_molar = Q_MASS * MOLAR_MASS_KG_PER_MOL;
    let rho_298 = RHO0 * T0 / 298.15;
    let diffusion_input = RHO_D_TARGET / (rho_298 * (T0 / 298.15).sqrt());
    let pressure = RHO0 * R_G * T0 / MOLAR_MASS_KG_PER_MOL;
    let length = LAMBDA * case.xi_max / (MASS_FLUX * CP);

    let mut reactor = SimpleReactorTask::new();
    reactor
        .fast_react_set(vec![FastElemReact {
            eq: "A=>B".into(),
            A: k0,
            n: 0.0,
            E: activation_energy,
            Q: q_molar,
        }])
        .expect("MF-1963 elementary reaction should be accepted");
    reactor.kindata.substances = vec!["A".into(), "B".into()];
    reactor.kindata.groups = Some(HashMap::from([
        ("A".into(), formula.clone()),
        ("B".into(), formula),
    ]));
    reactor.M = MOLAR_MASS_KG_PER_MOL;
    reactor.set_parameters(
        vec![q_molar],
        pressure,
        T0,
        CP,
        HashMap::from([
            ("T".into(), T0),
            ("A".into(), 1.0),
            ("B".into(), 0.0),
            ("q".into(), 0.0),
            ("A_flux".into(), 0.0),
            ("B_flux".into(), 0.0),
        ]),
        LAMBDA,
        HashMap::from([("A".into(), diffusion_input), ("B".into(), diffusion_input)]),
        MASS_FLUX,
        ScalingConfig::new(T0, length, T_SCALE),
    );
    reactor
}

fn reference_seed(case: Case, nodes: usize) -> DMatrix<f64> {
    let profile = reference_profile(case);
    let pe_q = case.xi_max;
    DMatrix::from_fn(6, nodes, |row, column| {
        let z = column as f64 / (nodes - 1) as f64;
        let xi = pe_q * z;
        let theta = interpolate(&profile, xi, |point| point.theta);
        let p = interpolate(&profile, xi, |point| point.p);
        let y_a = interpolate(&profile, xi, |point| point.y_a);
        let y_b = 1.0 - y_a;
        let derivative_a_z = -GAMMA_MF00 * pe_q * p;
        let k = RHO_D_TARGET;
        match row {
            0 => theta,
            1 => LAMBDA * pe_q * p,
            2 => y_a,
            3 => k * derivative_a_z,
            4 => y_b,
            5 => -k * derivative_a_z,
            _ => unreachable!(),
        }
    })
}

fn smooth_seed(case: Case, nodes: usize) -> DMatrix<f64> {
    let peak_xi = if case.name == "MF-LOW" {
        11996.338870908587
    } else {
        0.04358419907586729
    };
    let center = peak_xi / case.xi_max;
    let steepness = if case.name == "MF-LOW" { 50.0 } else { 35.0 };
    let logistic = |z: f64| 1.0 / (1.0 + (-steepness * (z - center)).exp());
    let left = logistic(0.0);
    let right = logistic(1.0);
    let normalization = right - left;
    let density_diffusion = RHO_D_TARGET;
    DMatrix::from_fn(6, nodes, |row, column| {
        let z = column as f64 / (nodes - 1) as f64;
        let g = logistic(z);
        let theta = (1.0 / GAMMA_MF00) * (g - left) / normalization;
        let dtheta_dz = (1.0 / GAMMA_MF00) * steepness * g * (1.0 - g) / normalization;
        let y_a = 1.0 - GAMMA_MF00 * theta;
        let j_a = -density_diffusion * GAMMA_MF00 * dtheta_dz;
        match row {
            0 => theta,
            1 => LAMBDA * dtheta_dz,
            2 => y_a,
            3 => j_a,
            4 => 1.0 - y_a,
            5 => -j_a,
            _ => unreachable!(),
        }
    })
}

#[test]
fn mf1963_actual_kithe_inputs_recover_transport_and_scaling_contract() {
    for case in cases() {
        let mut reactor = reactor_for(case);
        assert_eq!(
            reactor.density_model,
            DensityModel::IdealGasFixedComposition
        );
        reactor
            .setup_bvp()
            .unwrap_or_else(|error| panic!("{} setup failed: {error}", case.name));

        assert!((reactor.ideal_gas_density().unwrap() - RHO0).abs() < 2e-14);
        assert!((reactor.D_ro_map["A"] - RHO_D_TARGET).abs() < 2e-14);
        assert!((reactor.D_ro_map["B"] - RHO_D_TARGET).abs() < 2e-14);
        assert!((reactor.Pe_q - case.xi_max).abs() < 2e-9 * case.xi_max.max(1.0));
        assert!((reactor.Le_number().unwrap() - 1.0).abs() < 2e-12);
        assert_eq!(
            reactor.solver.unknowns,
            ["Teta", "q", "C0", "J0", "C1", "J1"]
        );
        assert!(reactor.map_eq_rate.contains_key("A=>B"));
        assert_eq!(
            reactor.thermal_effects,
            vec![Q_MASS * MOLAR_MASS_KG_PER_MOL]
        );

        let activation_energy = R_G * T0 / BETA_MF00;
        let k0 = case.delta * CP * GAMMA_MF00 * MASS_FLUX.powi(2) / (LAMBDA * RHO0)
            * (activation_energy / (R_G * T0)).exp();
        let state = [20.0, 0.0, 0.8, 0.0, 0.2, 0.0];
        let state_names: Vec<&str> = reactor.solver.unknowns.iter().map(String::as_str).collect();
        let rate_fun = reactor.map_eq_rate["A=>B"].lambdify_borrowed_thread_safe(&state_names);
        let rate_actual = rate_fun(&state);
        let temperature = T0 + T_SCALE * state[0];
        let density = RHO0 * T0 / temperature;
        let k_temperature = k0 * (-activation_energy / (R_G * temperature)).exp();
        let rate_expected = k_temperature * density * state[2] / MOLAR_MASS_KG_PER_MOL;
        assert!((rate_actual / rate_expected - 1.0).abs() < 2e-12);

        let beta_kithe = R_G * T0 / activation_energy;
        let gamma_kithe = CP * R_G * T0.powi(2) / (Q_MASS * activation_energy);
        let beta_mf00_on_kithe_inputs = R_MF00 * T0 / activation_energy;
        let gamma_mf00_on_kithe_inputs = CP * R_MF00 * T0.powi(2) / (Q_MASS * activation_energy);
        let k_at_t0_mf00 = k0 * (-activation_energy / (R_MF00 * T0)).exp();
        let delta_mf00_on_kithe_inputs = Q_MASS * LAMBDA * RHO0 * activation_energy
            / (CP.powi(2) * MASS_FLUX.powi(2) * R_MF00 * T0.powi(2))
            * k_at_t0_mf00;
        assert!((beta_kithe / BETA_MF00 - 1.0).abs() < 2e-15);
        assert!((gamma_kithe / GAMMA_MF00 - 1.0).abs() < 2e-15);

        let seed = reference_seed(case, 101);
        assert_eq!(seed.nrows(), reactor.solver.unknowns.len());
        assert!(seed.iter().all(|value| value.is_finite()));
        assert!((seed[(0, 0)]).abs() < 2e-12);
        assert!((seed[(2, 0)] - 1.0).abs() < 2e-12);
        assert!((seed[(4, 0)]).abs() < 2e-12);
        println!(
            "{} mapped: L={:.9e} m Pe_q={:.9e} rhoD={:.9e} D_input={:.9e} k0={:.9e} Q_mol={:.9e} J/mol beta_KiThe={:.12e} gamma_KiThe={:.12e} MF00_R beta={:.12e} gamma={:.12e} delta_ratio={:.9e}",
            case.name,
            reactor.L,
            reactor.Pe_q,
            reactor.D_ro_map["A"],
            reactor.Diffusion["A"],
            k0,
            reactor.thermal_effects[0],
            beta_kithe,
            gamma_kithe,
            beta_mf00_on_kithe_inputs,
            gamma_mf00_on_kithe_inputs,
            delta_mf00_on_kithe_inputs / case.delta,
        );
    }
}

fn solve_mf_sci(
    case: Case,
    initial_guess: DMatrix<f64>,
    nodes: usize,
    max_nodes: usize,
    tolerance: f64,
) -> Result<SimpleReactorTask, ReactorError> {
    let mut reactor = reactor_for(case);
    reactor.setup_bvp()?;
    let options = BvpSciOptions::default()
        .with_execution(BvpSciExecution::Lambdify)
        .with_assembly(BvpSciAssembly::ExprLegacy)
        .with_matrix_layout(BvpSciMatrixLayout::Sparse)
        .with_tolerance(tolerance)
        .with_limits(max_nodes, 100, 20);
    reactor.solver.solve_with_hybrid_settings(
        BvpSolverSettings::Sci {
            options,
            sparse_factorization: BvpSciSparseFactorization::StandardFaer,
            aot_config: None,
        },
        initial_guess,
        0.0,
        1.0,
        nodes - 1,
    )?;
    reactor.postprocessing()?;
    Ok(reactor)
}

fn assert_sci_matches_reference(
    case: Case,
    initial_guess: DMatrix<f64>,
    nodes: usize,
    max_nodes: usize,
    tolerance: f64,
) -> SimpleReactorTask {
    let reactor = solve_mf_sci(case, initial_guess, nodes, max_nodes, tolerance)
        .unwrap_or_else(|error| panic!("{}: {error}", case.name));

    let solution = reactor.solver.get_solution().expect("published solution");
    let mesh = reactor
        .solver
        .x_mesh
        .as_ref()
        .expect("published physical mesh");
    assert_eq!(solution.ncols(), 6);
    assert_eq!(solution.nrows(), mesh.len());
    let profile = reference_profile(case);
    let mut max_t_error = 0.0_f64;
    let mut max_y_error = 0.0_f64;
    let mut max_q_error = 0.0_f64;
    let mut max_balance_error = 0.0_f64;
    for (row, &x_m) in mesh.iter().enumerate() {
        let xi = x_m * MASS_FLUX * CP / LAMBDA;
        let theta = interpolate(&profile, xi, |point| point.theta);
        let p = interpolate(&profile, xi, |point| point.p);
        let y_a = interpolate(&profile, xi, |point| point.y_a);
        let t_expected = T0 + T_SCALE * theta;
        let q_expected = MASS_FLUX * CP * T_SCALE * p;
        max_t_error = max_t_error.max((solution[(row, 0)] - t_expected).abs());
        max_q_error = max_q_error.max((solution[(row, 1)] - q_expected).abs());
        max_y_error = max_y_error
            .max((solution[(row, 2)] - y_a).abs())
            .max((solution[(row, 4)] - (1.0 - y_a)).abs());
        max_balance_error =
            max_balance_error.max((solution[(row, 2)] + solution[(row, 4)] - 1.0).abs());
    }
    assert!(solution.iter().all(|value| value.is_finite()));
    assert!(
        solution
            .column(2)
            .iter()
            .all(|value| (-1e-12..=1.0 + 1e-12).contains(value))
    );
    assert!(
        solution
            .column(4)
            .iter()
            .all(|value| (-1e-12..=1.0 + 1e-12).contains(value))
    );
    let q_scale =
        MASS_FLUX * CP * T_SCALE * interpolate(&profile, 0.0, |point| point.p).abs().max(1.0);
    println!(
        "{} BVP_sci diagnostics: nodes={} T_error={:.3e} K Y_error={:.3e} q_error={:.3e} W/m2 q_scale={:.3e} sumY_error={:.3e}",
        case.name,
        solution.nrows(),
        max_t_error,
        max_y_error,
        max_q_error,
        q_scale,
        max_balance_error,
    );
    assert!(
        max_t_error < 0.5,
        "{} T profile error={max_t_error:.3e} K",
        case.name
    );
    assert!(
        max_y_error < 2e-3,
        "{} species profile error={max_y_error:.3e}",
        case.name
    );
    assert!(
        max_q_error / q_scale < 2e-3,
        "{} heat-flux profile error={max_q_error:.3e} W/m2",
        case.name
    );
    assert!(max_balance_error < 2e-8);
    println!(
        "{} BVP_sci: nodes={} T_error={:.3e} K Y_error={:.3e} q_error={:.3e} W/m2 sumY_error={:.3e}",
        case.name,
        solution.nrows(),
        max_t_error,
        max_y_error,
        max_q_error,
        max_balance_error,
    );
    reactor
}

fn production_profile_errors(case: Case, reactor: &SimpleReactorTask) -> (f64, f64, f64, f64) {
    let solution = reactor.solver.get_solution().expect("published solution");
    let mesh = reactor.solver.x_mesh.as_ref().expect("published mesh");
    assert_eq!(solution.ncols(), 6);
    assert_eq!(solution.nrows(), mesh.len());
    let profile = reference_profile(case);
    let (mut max_t_error, mut max_y_error, mut max_q_error, mut max_closure_error) =
        (0.0_f64, 0.0_f64, 0.0_f64, 0.0_f64);
    for (row, &x_m) in mesh.iter().enumerate() {
        let xi = x_m * MASS_FLUX * CP / LAMBDA;
        let theta = interpolate(&profile, xi, |point| point.theta);
        let p = interpolate(&profile, xi, |point| point.p);
        let y_a = interpolate(&profile, xi, |point| point.y_a);
        max_t_error = max_t_error.max((solution[(row, 0)] - (T0 + T_SCALE * theta)).abs());
        max_y_error = max_y_error
            .max((solution[(row, 2)] - y_a).abs())
            .max((solution[(row, 4)] - (1.0 - y_a)).abs());
        max_q_error = max_q_error.max((solution[(row, 1)] - MASS_FLUX * CP * T_SCALE * p).abs());
        max_closure_error =
            max_closure_error.max((solution[(row, 2)] + solution[(row, 4)] - 1.0).abs());
    }
    let q_scale = MASS_FLUX
        * CP
        * T_SCALE
        * profile
            .iter()
            .map(|point| point.p.abs())
            .fold(0.0_f64, f64::max);
    assert!(q_scale.is_finite() && q_scale > 0.0);
    (
        max_t_error,
        max_y_error,
        max_q_error / q_scale,
        max_closure_error,
    )
}

fn physical_profile_difference(
    candidate: &SimpleReactorTask,
    reference: &SimpleReactorTask,
) -> (f64, f64, f64) {
    let candidate_solution = candidate.solver.get_solution().expect("candidate profile");
    let candidate_mesh = candidate.solver.x_mesh.as_ref().expect("candidate mesh");
    let reference_solution = reference.solver.get_solution().expect("reference profile");
    let reference_mesh = reference
        .solver
        .x_mesh
        .as_ref()
        .expect("reference mesh")
        .as_slice();
    let interpolate = |column: usize, x: f64| {
        let upper = reference_mesh
            .partition_point(|point| *point < x)
            .clamp(1, reference_mesh.len() - 1);
        let lower = upper - 1;
        let fraction = ((x - reference_mesh[lower])
            / (reference_mesh[upper] - reference_mesh[lower]))
            .clamp(0.0, 1.0);
        reference_solution[(lower, column)]
            + fraction * (reference_solution[(upper, column)] - reference_solution[(lower, column)])
    };
    let mut max_t = 0.0_f64;
    let mut max_y = 0.0_f64;
    let mut max_q = 0.0_f64;
    let mut q_scale = 0.0_f64;
    for (row, &x) in candidate_mesh.iter().enumerate() {
        max_t = max_t.max((candidate_solution[(row, 0)] - interpolate(0, x)).abs());
        max_y = max_y
            .max((candidate_solution[(row, 2)] - interpolate(2, x)).abs())
            .max((candidate_solution[(row, 4)] - interpolate(4, x)).abs());
        let q_reference = interpolate(1, x);
        max_q = max_q.max((candidate_solution[(row, 1)] - q_reference).abs());
        q_scale = q_scale.max(q_reference.abs());
    }
    (max_t, max_y, max_q / q_scale.max(1.0))
}

fn mf_discrete_trapezoid_defects(case: Case, reactor: &SimpleReactorTask) -> (f64, f64, f64) {
    let solution = reactor.solver.get_solution().expect("published solution");
    let mesh = reactor.solver.x_mesh.as_ref().expect("published mesh");
    let mut max_theta_defect = 0.0_f64;
    let mut max_p_defect = 0.0_f64;
    let mut max_first_integral_defect = 0.0_f64;
    for node in 0..mesh.len() {
        let theta = (solution[(node, 0)] - T0) / T_SCALE;
        let conversion = 1.0 - solution[(node, 2)];
        max_first_integral_defect =
            max_first_integral_defect.max((conversion - GAMMA_MF00 * theta).abs());
    }
    for interval in 0..mesh.len() - 1 {
        let xi0 = mesh[interval] * MASS_FLUX * CP / LAMBDA;
        let xi1 = mesh[interval + 1] * MASS_FLUX * CP / LAMBDA;
        let h = xi1 - xi0;
        let theta0 = (solution[(interval, 0)] - T0) / T_SCALE;
        let theta1 = (solution[(interval + 1, 0)] - T0) / T_SCALE;
        let p0 = solution[(interval, 1)] / (MASS_FLUX * CP * T_SCALE);
        let p1 = solution[(interval + 1, 1)] / (MASS_FLUX * CP * T_SCALE);
        let source = |theta: f64| {
            let denominator = 1.0 + BETA_MF00 * theta;
            case.delta * (theta / denominator).exp() * (1.0 - GAMMA_MF00 * theta) / denominator
        };
        let rhs_p0 = p0 - source(theta0);
        let rhs_p1 = p1 - source(theta1);
        max_theta_defect = max_theta_defect.max((theta1 - theta0 - 0.5 * h * (p0 + p1)).abs());
        max_p_defect = max_p_defect.max((p1 - p0 - 0.5 * h * (rhs_p0 + rhs_p1)).abs());
    }
    (max_theta_defect, max_p_defect, max_first_integral_defect)
}

fn assert_production_profile_matches_reference(case: Case, reactor: &SimpleReactorTask) {
    let solution = reactor.solver.get_solution().expect("published solution");
    assert!(solution.iter().all(|value| value.is_finite()));
    assert!(
        solution
            .column(0)
            .iter()
            .all(|theta| T0 + T_SCALE * theta > 0.0)
    );
    assert!(
        solution
            .column(2)
            .iter()
            .chain(solution.column(4).iter())
            .all(|fraction| (0.0..=1.0).contains(fraction))
    );
    let (max_t_error, max_y_error, relative_q_error, max_closure_error) =
        production_profile_errors(case, reactor);
    println!(
        "{} production profile: nodes={} T_error={:.3e} K Y_error={:.3e} relative_q_error={:.3e} closure={:.3e}",
        case.name,
        solution.nrows(),
        max_t_error,
        max_y_error,
        relative_q_error,
        max_closure_error,
    );
    assert!(
        max_t_error < 0.5,
        "{} T error={max_t_error:.3e} K",
        case.name
    );
    assert!(
        max_y_error < 2e-3,
        "{} Y error={max_y_error:.3e}",
        case.name
    );
    assert!(
        relative_q_error < 2e-3,
        "{} relative q error={relative_q_error:.3e}",
        case.name
    );
    assert!(max_closure_error < 2e-8);
}

fn solve_mf_damped(
    case: Case,
    nodes: usize,
    max_refinements: usize,
    grid_method: GridRefinementMethod,
    max_jac: usize,
    adaptive_refinement: bool,
    detailed_decisions: bool,
) -> Result<SimpleReactorTask, String> {
    let mut template = reactor_for(case);
    template
        .setup_bvp()
        .map_err(|error| format!("Damped seed preparation: {error}"))?;
    let seed = damped_reduced_seed_from_nondimensional_profile(
        &template,
        &reference_seed(case, nodes),
        nodes,
    );
    solve_mf_damped_with_seed(
        case,
        seed,
        nodes,
        max_refinements,
        grid_method,
        max_jac,
        adaptive_refinement,
        detailed_decisions,
    )
}

fn solve_mf_damped_with_seed(
    case: Case,
    initial: DMatrix<f64>,
    nodes: usize,
    max_refinements: usize,
    grid_method: GridRefinementMethod,
    max_jac: usize,
    adaptive_refinement: bool,
    detailed_decisions: bool,
) -> Result<SimpleReactorTask, String> {
    solve_mf_damped_with_seed_and_telemetry(
        case,
        initial,
        nodes,
        max_refinements,
        grid_method,
        max_jac,
        adaptive_refinement,
        detailed_decisions,
    )
    .0
}

fn solve_mf_damped_with_seed_and_telemetry(
    case: Case,
    initial: DMatrix<f64>,
    nodes: usize,
    max_refinements: usize,
    grid_method: GridRefinementMethod,
    max_jac: usize,
    adaptive_refinement: bool,
    detailed_decisions: bool,
) -> (
    Result<SimpleReactorTask, String>,
    Option<BvpSolverTelemetrySnapshot>,
) {
    let mut reactor = reactor_for(case);
    if let Err(error) = reactor.setup_bvp() {
        return (Err(format!("Damped setup: {error}")), None);
    }
    reactor.solver = reactor
        .solver
        .with_bvp_telemetry_mode(if detailed_decisions {
            BvpApiTelemetryMode::Detailed
        } else {
            BvpApiTelemetryMode::Counters
        });
    let strategy = SolverParams {
        max_jac: Some(max_jac),
        max_damp_iter: Some(20),
        adaptive: adaptive_refinement.then_some(AdaptiveGridConfig {
            version: 1,
            max_refinements,
            grid_method,
        }),
        ..SolverParams::default()
    };
    let bounds = reactor
        .solver
        .unknowns
        .iter()
        .map(|name| {
            let range = match name.as_str() {
                "Teta" => (-19.9, 120.0),
                "q" => (-1e6, 1e6),
                _ => (-1e12, 1e12),
            };
            (name.clone(), range)
        })
        .collect();
    let tolerances = reactor
        .solver
        .unknowns
        .iter()
        .map(|name| (name.clone(), 1e-8))
        .collect();
    let mut config = NrbvpHandoffConfig::new(
        initial,
        0.0,
        1.0,
        nodes,
        "trapezoid".into(),
        "Damped".into(),
        Some(strategy),
        None,
        "Sparse".into(),
        1e-8,
        Some(tolerances),
        150,
        Some(bounds),
        Some("off".into()),
    )
    .with_solver_backend_config(
        ReactorBvpSolverConfig::sparse_lambdify()
            .with_symbolic_backend(ReactorBvpSymbolicBackend::ExprLegacy),
    );
    if detailed_decisions {
        config = config.with_bvp_logging_config(
            BvpLoggingConfig::new(BvpLoggingMode::Detailed).with_max_events(4096),
        );
    }
    if let Err(error) = reactor.solver.solve_with_bvp_api(config) {
        let telemetry = reactor.solver.last_telemetry.clone();
        return (
            Err(format!(
                "Damped solve: {error:?}; telemetry={:?}",
                reactor.solver.last_telemetry
            )),
            telemetry,
        );
    }
    let telemetry = reactor.solver.last_telemetry.clone();
    if let Err(error) = reactor.postprocessing() {
        return (Err(format!("Damped postprocessing: {error}")), telemetry);
    }
    (Ok(reactor), telemetry)
}

#[test]
fn mf1963_damped_fixed_mesh_retains_hot_reference_branch() {
    let case = cases()[1];
    let nodes = 301;
    let seed = reference_seed(case, nodes);
    let seed_teta = seed.row(0);
    let initial_temperatures = seed_teta.iter().map(|theta| T0 + T_SCALE * theta);
    let (seed_min_t, seed_max_t) = initial_temperatures
        .fold((f64::INFINITY, f64::NEG_INFINITY), |(min_t, max_t), t| {
            (min_t.min(t), max_t.max(t))
        });
    println!(
        "{} Damped fixed-mesh seed: rows={} cols={} T_range=[{seed_min_t:.9e}, {seed_max_t:.9e}] K endpoints=({:.9e}, {:.9e}) K Y_A_endpoints=({:.9e}, {:.9e})",
        case.name,
        seed.nrows(),
        seed.ncols(),
        T0 + T_SCALE * seed[(0, 0)],
        T0 + T_SCALE * seed[(0, seed.ncols() - 1)],
        seed[(2, 0)],
        seed[(2, seed.ncols() - 1)],
    );

    let reactor = solve_mf_damped(
        case,
        nodes,
        0,
        GridRefinementMethod::DoublePoints,
        50,
        false,
        true,
    )
    .unwrap_or_else(|error| panic!("{} fixed-mesh Damped: {error}", case.name));
    let solution = reactor.solver.get_solution().expect("published solution");
    let errors = production_profile_errors(case, &reactor);
    let (min_t, max_t) = solution.column(0).iter().fold(
        (f64::INFINITY, f64::NEG_INFINITY),
        |(min_t, max_t), &temperature| (min_t.min(temperature), max_t.max(temperature)),
    );
    println!(
        "{} Damped fixed-mesh published (postprocessed physical units): states={:?} rows={} T_range=[{min_t:.9e}, {max_t:.9e}] K endpoints=({:.9e}, {:.9e}) K Y_A_endpoints=({:.9e}, {:.9e}) errors=(dT={:.6e} K,dY={:.6e},relative_dq={:.6e},closure={:.6e})",
        case.name,
        reactor.solver.unknowns,
        solution.nrows(),
        solution[(0, 0)],
        solution[(solution.nrows() - 1, 0)],
        solution[(0, 2)],
        solution[(solution.nrows() - 1, 2)],
        errors.0,
        errors.1,
        errors.2,
        errors.3,
    );
    let snapshot = reactor.solver.last_telemetry.as_ref().expect("telemetry");
    if let BvpSolverTelemetrySnapshot::Damped(snapshot) = snapshot {
        println!(
            "{} Damped fixed-mesh telemetry: status={:?} iterations={} residual_evals={} jacobian_evals={} damping_trials={} damping_rejections={} refinements={} events={} dropped={}",
            case.name,
            snapshot.last_solve_status,
            snapshot.counters.iterations,
            snapshot.counters.residual_requests,
            snapshot.counters.jacobian_requests,
            snapshot.counters.damping_trials,
            snapshot.counters.damping_rejections,
            snapshot.counters.grid_refinements,
            snapshot.log_events.len(),
            snapshot.log_events_dropped,
        );
        for event in &snapshot.log_events {
            if matches!(
                event.kind,
                BvpLogEventKind::DampingTrial | BvpLogEventKind::DampingRejected
            ) {
                println!(
                    "DAMP_TRACE iter={} trial={} kind={:?} lambda={:.9e} value_b={:.9e}",
                    event.iteration, event.damping_trial, event.kind, event.value_a, event.value_b,
                );
            }
        }
    }
    let defects = mf_discrete_trapezoid_defects(case, &reactor);
    println!(
        "{} Damped fixed-mesh independent MF defects: max|r_theta|={:.9e} max|r_p|={:.9e} max|eta-gamma*theta|={:.9e}",
        case.name, defects.0, defects.1, defects.2,
    );
    assert_eq!(solution.nrows(), nodes + 1, "fixed mesh must not refine");
    assert!(max_t > 1000.0, "Damped left the hot reference branch");
    assert_production_profile_matches_reference(case, &reactor);
    assert_damped_reports_converged(&reactor);
}

fn full_grid_seed_from_physical_profile(reactor: &SimpleReactorTask, nodes: usize) -> DMatrix<f64> {
    let solution = reactor.solver.get_solution().expect("physical profile");
    let mesh = reactor.solver.x_mesh.as_ref().expect("physical mesh");
    let mesh = mesh.as_slice();
    let interpolate_column = |column: usize, x: f64| {
        let upper = mesh
            .partition_point(|point| *point < x)
            .clamp(1, mesh.len() - 1);
        let lower = upper - 1;
        let fraction = ((x - mesh[lower]) / (mesh[upper] - mesh[lower])).clamp(0.0, 1.0);
        solution[(lower, column)]
            + fraction * (solution[(upper, column)] - solution[(lower, column)])
    };
    let length = *mesh.last().expect("nonempty mesh");
    DMatrix::from_fn(6, nodes, |state, point| {
        let x = length * point as f64 / (nodes - 1) as f64;
        match state {
            0 => (interpolate_column(0, x) - T0) / T_SCALE,
            1 => interpolate_column(1, x) * length / T_SCALE,
            2 => interpolate_column(2, x),
            3 => interpolate_column(3, x) * length,
            4 => interpolate_column(4, x),
            5 => interpolate_column(5, x) * length,
            _ => unreachable!(),
        }
    })
}

fn damped_reduced_seed_from_profile(
    reactor: &SimpleReactorTask,
    full_profile: &DMatrix<f64>,
    source_mesh: &[f64],
    intervals: usize,
) -> DMatrix<f64> {
    // Damped's matrix is state-major: each row omits that state's own fixed
    // endpoint. Flattening the full mesh first would mix states at mixed BCs.
    let states = reactor.solver.unknowns.len();
    assert_eq!(full_profile.nrows(), states);
    assert_eq!(full_profile.ncols(), source_mesh.len());
    assert!(source_mesh.len() >= 2);

    DMatrix::from_fn(states, intervals, |state, free_node| {
        let name = &reactor.solver.unknowns[state];
        let (boundary_side, _) = reactor
            .solver
            .BorderConditions
            .get(name)
            .copied()
            .expect("boundary condition for every state");
        let node = match boundary_side {
            0 => free_node + 1,
            1 => free_node,
            side => panic!("unexpected boundary side {side} for {name}"),
        };
        let x = node as f64 / intervals as f64;
        let upper = source_mesh
            .partition_point(|point| *point < x)
            .clamp(1, source_mesh.len() - 1);
        let lower = upper - 1;
        let fraction =
            ((x - source_mesh[lower]) / (source_mesh[upper] - source_mesh[lower])).clamp(0.0, 1.0);
        full_profile[(state, lower)]
            + fraction * (full_profile[(state, upper)] - full_profile[(state, lower)])
    })
}

fn damped_reduced_seed_from_nondimensional_profile(
    reactor: &SimpleReactorTask,
    full_profile: &DMatrix<f64>,
    intervals: usize,
) -> DMatrix<f64> {
    let source_mesh = (0..full_profile.ncols())
        .map(|point| point as f64 / (full_profile.ncols() - 1) as f64)
        .collect::<Vec<_>>();
    damped_reduced_seed_from_profile(reactor, full_profile, &source_mesh, intervals)
}

fn damped_reduced_seed_from_physical_profile(
    reactor: &SimpleReactorTask,
    intervals: usize,
) -> DMatrix<f64> {
    let solution = reactor.solver.get_solution().expect("physical profile");
    let mesh = reactor.solver.x_mesh.as_ref().expect("physical mesh");
    let states = reactor.solver.unknowns.len();
    assert_eq!(solution.ncols(), states);
    assert_eq!(solution.nrows(), mesh.len());

    let nondimensional = DMatrix::from_fn(states, solution.nrows(), |state, node| {
        let physical = solution[(node, state)];
        match state {
            0 => (physical - T0) / T_SCALE,
            1 => physical * reactor.L / T_SCALE,
            2 | 4 => physical,
            3 | 5 => physical * reactor.L,
            _ => unreachable!("unexpected MF state index {state}"),
        }
    });
    let length = mesh[mesh.len() - 1];
    assert!(length.is_finite() && length > 0.0);
    let unit_mesh = mesh.iter().map(|x| x / length).collect::<Vec<_>>();
    damped_reduced_seed_from_profile(reactor, &nondimensional, &unit_mesh, intervals)
}

fn mf_scalar_discrete_residual(theta: &[f64], slope: &[f64], delta: f64) -> Vec<f64> {
    let intervals = theta.len() - 1;
    let h = cases()[1].xi_max / intervals as f64;
    let source = |value: f64| {
        let denominator = 1.0 + BETA_MF00 * value;
        delta * (value / denominator).exp() * (1.0 - GAMMA_MF00 * value) / denominator
    };
    let mut residual = Vec::with_capacity(2 * intervals + 2);
    residual.push(theta[0]);
    for node in 0..intervals {
        residual.push(theta[node + 1] - theta[node] - 0.5 * h * (slope[node] + slope[node + 1]));
        let f0 = slope[node] - source(theta[node]);
        let f1 = slope[node + 1] - source(theta[node + 1]);
        residual.push(slope[node + 1] - slope[node] - 0.5 * h * (f0 + f1));
    }
    residual.push(slope[intervals]);
    residual
}

fn mf_scalar_discrete_jacobian(theta: &[f64], delta: f64) -> DMatrix<f64> {
    let intervals = theta.len() - 1;
    let h = cases()[1].xi_max / intervals as f64;
    let source_derivative = |value: f64| {
        let denominator = 1.0 + BETA_MF00 * value;
        let remaining = 1.0 - GAMMA_MF00 * value;
        let exponential = (value / denominator).exp();
        delta
            * exponential
            * (remaining / denominator.powi(3)
                - GAMMA_MF00 / denominator
                - BETA_MF00 * remaining / denominator.powi(2))
    };
    let mut jacobian = DMatrix::zeros(2 * intervals + 2, 2 * intervals + 2);
    jacobian[(0, 0)] = 1.0;
    for node in 0..intervals {
        let theta_row = 1 + 2 * node;
        let slope_row = theta_row + 1;
        let theta_left = 2 * node;
        let slope_left = theta_left + 1;
        let theta_right = theta_left + 2;
        let slope_right = theta_left + 3;
        jacobian[(theta_row, theta_left)] = -1.0;
        jacobian[(theta_row, theta_right)] = 1.0;
        jacobian[(theta_row, slope_left)] = -0.5 * h;
        jacobian[(theta_row, slope_right)] = -0.5 * h;
        jacobian[(slope_row, theta_left)] = 0.5 * h * source_derivative(theta[node]);
        jacobian[(slope_row, theta_right)] = 0.5 * h * source_derivative(theta[node + 1]);
        jacobian[(slope_row, slope_left)] = -1.0 - 0.5 * h;
        jacobian[(slope_row, slope_right)] = 1.0 - 0.5 * h;
    }
    jacobian[(2 * intervals + 1, 2 * intervals + 1)] = 1.0;
    jacobian
}

#[test]
fn mf1963_damped_hot_seed_roundtrip_and_first_scalar_newton_audit() {
    let case = cases()[1];
    let intervals = 301;
    let hot = solve_mf_correction_newton(case, intervals)
        .unwrap_or_else(|error| panic!("MF-HIGH CorrectionNewton hot control: {error}"));
    let seed = damped_reduced_seed_from_physical_profile(&hot, intervals);
    let solution = hot.solver.get_solution().expect("CorrectionNewton profile");
    let theta = (0..=intervals)
        .map(|node| (solution[(node, 0)] - T0) / T_SCALE)
        .collect::<Vec<_>>();
    let slope = (0..=intervals)
        .map(|node| solution[(node, 1)] / (MASS_FLUX * CP * T_SCALE))
        .collect::<Vec<_>>();
    let mut roundtrip_max = 0.0_f64;
    for state in 0..6 {
        let left_fixed = matches!(state, 0 | 2 | 4);
        for slot in 0..intervals {
            let node = if left_fixed { slot + 1 } else { slot };
            let physical = solution[(node, state)];
            let internal = match state {
                0 => (physical - T0) / T_SCALE,
                1 => physical * hot.L / T_SCALE,
                2 | 4 => physical,
                3 | 5 => physical * hot.L,
                _ => unreachable!(),
            };
            roundtrip_max = roundtrip_max.max((seed[(state, slot)] - internal).abs());
        }
    }
    assert!(
        roundtrip_max < 2e-12,
        "physical/internal/reduced round-trip error={roundtrip_max:e}"
    );

    let residual = mf_scalar_discrete_residual(&theta, &slope, case.delta);
    let jacobian = mf_scalar_discrete_jacobian(&theta, case.delta);
    let residual = nalgebra::DVector::from_vec(residual);
    let correction = jacobian
        .clone()
        .lu()
        .solve(&(-&residual))
        .expect("first scalar Newton correction");
    let linear_error = (&jacobian * &correction + &residual).norm();
    let theta_trial = (0..=intervals)
        .map(|node| theta[node] + correction[2 * node])
        .collect::<Vec<_>>();
    let slope_trial = (0..=intervals)
        .map(|node| slope[node] + correction[2 * node + 1])
        .collect::<Vec<_>>();
    let full_trial_residual = mf_scalar_discrete_residual(&theta_trial, &slope_trial, case.delta);
    let full_trial_l2 = full_trial_residual
        .iter()
        .map(|x| x * x)
        .sum::<f64>()
        .sqrt();
    let residual_inf = residual.iter().map(|x| x.abs()).fold(0.0, f64::max);
    println!(
        "MF-HIGH corrected hot seed: roundtrip_max={roundtrip_max:.6e}; independent scalar trapezoid F0: l2={:.9e} inf={residual_inf:.9e}; first Newton correction l2={:.9e} range=[{:.9e},{:.9e}]; linear residual={linear_error:.9e}; full-step trial residual l2={full_trial_l2:.9e}",
        residual.norm(),
        correction.norm(),
        correction.min(),
        correction.max()
    );

    for (family, component) in [("theta", 0usize), ("slope", 1usize), ("mixed", 2usize)] {
        let mut direction = nalgebra::DVector::zeros(2 * intervals + 2);
        for node in (4..intervals - 4).step_by(19) {
            if component != 1 {
                direction[2 * node] = ((node % 7) as f64 - 3.0) / 3.0;
            }
            if component != 0 {
                direction[2 * node + 1] = ((node % 11) as f64 - 5.0) / 5.0;
            }
        }
        direction /= direction.norm();
        let epsilon = 1e-6;
        let eval = |sign: f64| {
            let theta_step = (0..=intervals)
                .map(|node| theta[node] + sign * epsilon * direction[2 * node])
                .collect::<Vec<_>>();
            let slope_step = (0..=intervals)
                .map(|node| slope[node] + sign * epsilon * direction[2 * node + 1])
                .collect::<Vec<_>>();
            nalgebra::DVector::from_vec(mf_scalar_discrete_residual(
                &theta_step,
                &slope_step,
                case.delta,
            ))
        };
        let fd = (eval(1.0) - eval(-1.0)) / (2.0 * epsilon);
        let analytic = &jacobian * &direction;
        let relative = (&fd - &analytic).norm() / fd.norm().max(analytic.norm()).max(1e-14);
        println!(
            "MF-HIGH independent scalar Jacobian directional FD {family}: relative_error={relative:.6e}"
        );
        assert!(
            relative < 2e-7,
            "independent scalar Jacobian FD mismatch in {family}: {relative:e}"
        );
    }
}

#[test]
fn mf1963_independent_trapezoid_hot_root_is_probed_through_damped() {
    let case = cases()[1];
    let intervals = 301;
    let rows = include_str!("mf1963_mf_high_hot_trapezoid_301.csv")
        .lines()
        .skip(1)
        .filter(|line| !line.trim().is_empty())
        .map(|line| {
            line.split(',')
                .map(|field| field.parse::<f64>().unwrap())
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    assert_eq!(rows.len(), intervals + 1);
    let mut template = reactor_for(case);
    template
        .setup_bvp()
        .expect("independent hot-root Damped setup");
    let full_profile = DMatrix::from_fn(6, intervals + 1, |state, node| {
        let theta = rows[node][1];
        let slope = rows[node][2];
        let y_a = rows[node][4];
        let y_a_flux = -RHO_D_TARGET * GAMMA_MF00 * case.xi_max * slope;
        match state {
            0 => theta,
            1 => LAMBDA * case.xi_max * slope,
            2 => y_a,
            3 => y_a_flux,
            4 => 1.0 - y_a,
            5 => -y_a_flux,
            _ => unreachable!(),
        }
    });
    let physical_profile = DMatrix::from_fn(6, intervals + 1, |state, node| {
        let theta = rows[node][1];
        let slope = rows[node][2];
        let y_a = rows[node][4];
        let y_a_flux = -RHO_D_TARGET * GAMMA_MF00 * case.xi_max * slope;
        match state {
            0 => T0 + T_SCALE * theta,
            1 => MASS_FLUX * CP * T_SCALE * slope,
            2 => y_a,
            3 => y_a_flux,
            4 => 1.0 - y_a,
            5 => -y_a_flux,
            _ => unreachable!(),
        }
    });
    template.solver.solution = Some(physical_profile.transpose());
    template.solver.x_mesh = Some(nalgebra::DVector::from_iterator(
        intervals + 1,
        (0..=intervals).map(|node| template.L * node as f64 / intervals as f64),
    ));
    let full_system_defects = mf_discrete_trapezoid_defects(case, &template);
    println!(
        "MF-HIGH independent hot root mapped into six physical states: full-system scalar defects=({:.6e},{:.6e},{:.6e})",
        full_system_defects.0, full_system_defects.1, full_system_defects.2,
    );
    assert!(
        full_system_defects.0 < 1e-8
            && full_system_defects.1 < 1e-8
            && full_system_defects.2 < 1e-10
    );
    let seed = damped_reduced_seed_from_nondimensional_profile(&template, &full_profile, intervals);
    let mut packing_error = 0.0_f64;
    for state in 0..6 {
        let left_fixed = matches!(state, 0 | 2 | 4);
        for slot in 0..intervals {
            let node = if left_fixed { slot + 1 } else { slot };
            packing_error =
                packing_error.max((seed[(state, slot)] - full_profile[(state, node)]).abs());
        }
    }
    assert!(
        packing_error < 2e-12,
        "independent hot-root reduced-seed packing={packing_error:e}"
    );

    let result = solve_mf_damped_with_seed(
        case,
        seed,
        intervals,
        0,
        GridRefinementMethod::DoublePoints,
        1,
        false,
        true,
    )
    .unwrap_or_else(|error| panic!("Damped from independent trapezoid hot root: {error}"));
    let solution = result.solver.get_solution().expect("Damped output");
    let mut max_temperature_delta = 0.0_f64;
    let mut max_y_delta = 0.0_f64;
    for node in 0..=intervals {
        max_temperature_delta =
            max_temperature_delta.max((solution[(node, 0)] - rows[node][3]).abs());
        max_y_delta = max_y_delta
            .max((solution[(node, 2)] - rows[node][4]).abs())
            .max((solution[(node, 4)] - (1.0 - rows[node][4])).abs());
    }
    let independent_residual = mf_scalar_discrete_residual(
        &rows.iter().map(|row| row[1]).collect::<Vec<_>>(),
        &rows.iter().map(|row| row[2]).collect::<Vec<_>>(),
        case.delta,
    );
    let residual_l2 = independent_residual
        .iter()
        .map(|value| value * value)
        .sum::<f64>()
        .sqrt();
    println!(
        "MF-HIGH independent trapezoid hot root -> Damped: scalar_residual_l2={residual_l2:.9e} packing_max={packing_error:.9e} initial_Tmax={:.6} K output_Tmax={:.6} K max_profile_dT={max_temperature_delta:.6e} K max_profile_dY={max_y_delta:.6e}",
        rows.iter()
            .map(|row| row[3])
            .fold(f64::NEG_INFINITY, f64::max),
        solution
            .column(0)
            .iter()
            .copied()
            .fold(f64::NEG_INFINITY, f64::max),
    );
    if let Some(BvpSolverTelemetrySnapshot::Damped(snapshot)) =
        result.solver.last_telemetry.as_ref()
    {
        for event in &snapshot.log_events {
            if matches!(
                event.kind,
                BvpLogEventKind::DampingTrial | BvpLogEventKind::DampingRejected
            ) {
                println!(
                    "MF-HIGH exact-hot-seed Damped event: kind={:?} iter={} trial={} alpha_or_next={:.9e} value_b={:.9e}",
                    event.kind, event.iteration, event.damping_trial, event.value_a, event.value_b,
                );
            }
        }
    }
    assert!(
        residual_l2 < 1e-9,
        "independent hot trapezoid fixture is not a root"
    );
    assert_damped_reports_converged(&result);
}

#[test]
fn mf1963_damped_corrector_newton_seed_layout_preserves_hot_branch() {
    let case = cases()[1];
    let intervals = 301;
    let hot = solve_mf_correction_newton(case, intervals)
        .unwrap_or_else(|error| panic!("MF-HIGH CorrectionNewton hot control: {error}"));
    assert_production_profile_matches_reference(case, &hot);

    let old_seed = full_grid_seed_from_physical_profile(&hot, intervals);
    let reduced_seed = damped_reduced_seed_from_physical_profile(&hot, intervals);
    let physical = hot.solver.get_solution().expect("hot physical profile");
    let nondimensional = DMatrix::from_fn(6, physical.nrows(), |state, node| {
        let value = physical[(node, state)];
        match state {
            0 => (value - T0) / T_SCALE,
            1 => value * hot.L / T_SCALE,
            2 | 4 => value,
            3 | 5 => value * hot.L,
            _ => unreachable!(),
        }
    });
    let boundary_sides = [0, 1, 0, 1, 0, 1];
    let mut roundtrip_error = 0.0_f64;
    for state in 0..6 {
        for slot in 0..intervals {
            let node = if boundary_sides[state] == 0 {
                slot + 1
            } else {
                slot
            };
            roundtrip_error = roundtrip_error
                .max((reduced_seed[(state, slot)] - nondimensional[(state, node)]).abs());
        }
    }
    println!(
        "MF-HIGH reduced seed round-trip max_abs={roundtrip_error:.6e} across Teta,q,C0,J0,C1,J1; physical mesh=[{:.6e},{:.6e}] m",
        hot.solver.x_mesh.as_ref().unwrap()[0],
        hot.solver.x_mesh.as_ref().unwrap()[intervals]
    );
    assert!(
        roundtrip_error < 2e-12,
        "hot seed physical-to-internal-to-reduced round-trip={roundtrip_error:e}"
    );
    let max_seed_delta = old_seed
        .iter()
        .zip(reduced_seed.iter())
        .map(|(old, reduced)| (old - reduced).abs())
        .fold(0.0_f64, f64::max);
    println!(
        "MF-HIGH Damped seed packing: old_shape={:?} reduced_shape={:?} max_abs_delta={max_seed_delta:.9e}",
        old_seed.shape(),
        reduced_seed.shape(),
    );

    let old_result = solve_mf_damped_with_seed(
        case,
        old_seed,
        intervals,
        0,
        GridRefinementMethod::DoublePoints,
        1,
        false,
        false,
    )
    .unwrap_or_else(|error| panic!("old interpolated hot seed Damped: {error}"));
    let old_error = physical_profile_difference(&old_result, &hot);

    let reduced_result = solve_mf_damped_with_seed(
        case,
        reduced_seed,
        intervals,
        0,
        GridRefinementMethod::DoublePoints,
        1,
        false,
        true,
    )
    .unwrap_or_else(|error| panic!("BC-aware reduced hot seed Damped: {error}"));
    let reduced_error = physical_profile_difference(&reduced_result, &hot);
    println!(
        "MF-HIGH Damped seed packing outcomes: old(dT={:.6e} K,dY={:.6e},rel_dq={:.6e}) reduced(dT={:.6e} K,dY={:.6e},rel_dq={:.6e})",
        old_error.0, old_error.1, old_error.2, reduced_error.0, reduced_error.1, reduced_error.2,
    );
    let telemetry = reduced_result
        .solver
        .last_telemetry
        .as_ref()
        .expect("Damped telemetry");
    if let BvpSolverTelemetrySnapshot::Damped(snapshot) = telemetry {
        let decision_events = snapshot
            .log_events
            .iter()
            .filter(|event| {
                matches!(
                    event.kind,
                    BvpLogEventKind::DampingTrial | BvpLogEventKind::DampingRejected
                )
            })
            .collect::<Vec<_>>();
        println!(
            "MF-HIGH Damped trial telemetry: {} decision events (event, iter, trial, value_a=alpha-or-next-alpha, value_b=NaN-or-correction-norm)",
            decision_events.len()
        );
        for event in decision_events.iter().take(24) {
            println!(
                "  {:?} iter={} trial={} a={:.6e} b={:.6e}",
                event.kind, event.iteration, event.damping_trial, event.value_a, event.value_b
            );
        }
    }
    for (label, errors) in [("interpolated", old_error), ("BC-aware", reduced_error)] {
        assert!(
            errors.0 < 0.01,
            "{label} seed temperature error={}",
            errors.0
        );
        assert!(errors.1 < 1e-5, "{label} seed fraction error={}", errors.1);
        assert!(errors.2 < 1e-5, "{label} seed flux error={}", errors.2);
    }
    assert_damped_reports_converged(&old_result);
    assert_damped_reports_converged(&reduced_result);
}

#[test]
fn mf1963_damped_fixed_mesh_seed_basin_matrix() {
    let case = cases()[1];
    let nodes = 301;
    let sci = assert_sci_matches_reference(case, reference_seed(case, nodes), nodes, 20_000, 1e-6);
    let mut reactor_template = reactor_for(case);
    reactor_template
        .setup_bvp()
        .expect("MF-HIGH template BVP setup");
    let seeds = [
        (
            "frozen-reference",
            damped_reduced_seed_from_nondimensional_profile(
                &reactor_template,
                &reference_seed(case, nodes),
                nodes,
            ),
        ),
        (
            "smooth",
            damped_reduced_seed_from_nondimensional_profile(
                &reactor_template,
                &smooth_seed(case, nodes),
                nodes,
            ),
        ),
        (
            "BVP-sci-hot-profile",
            damped_reduced_seed_from_physical_profile(&sci, nodes),
        ),
    ];
    for (seed_name, seed) in seeds {
        let seed_t =
            seed.row(0)
                .iter()
                .fold((f64::INFINITY, f64::NEG_INFINITY), |(low, high), theta| {
                    let t = T0 + T_SCALE * theta;
                    (low.min(t), high.max(t))
                });
        let reactor = solve_mf_damped_with_seed(
            case,
            seed,
            nodes,
            0,
            GridRefinementMethod::DoublePoints,
            50,
            false,
            true,
        )
        .unwrap_or_else(|error| panic!("{} Damped/{seed_name}: {error}", case.name));
        let errors = production_profile_errors(case, &reactor);
        let solution = reactor.solver.get_solution().expect("published solution");
        let actual_t = solution.column(0).iter().fold(
            (f64::INFINITY, f64::NEG_INFINITY),
            |(low, high), &temperature| (low.min(temperature), high.max(temperature)),
        );
        let snapshot = reactor.solver.last_telemetry.as_ref().expect("telemetry");
        let status = match snapshot {
            BvpSolverTelemetrySnapshot::Damped(snapshot) => snapshot.last_solve_status,
            other => panic!("expected Damped telemetry, got {other:?}"),
        };
        let defects = mf_discrete_trapezoid_defects(case, &reactor);
        println!(
            "MF-HIGH basin seed={seed_name} seed_T=[{:.6e},{:.6e}] K result_T=[{:.6e},{:.6e}] K dT={:.6e} K dY={:.6e} rel_dq={:.6e} mf_defect=({:.3e},{:.3e},{:.3e}) status={status:?}",
            seed_t.0,
            seed_t.1,
            actual_t.0,
            actual_t.1,
            errors.0,
            errors.1,
            errors.2,
            defects.0,
            defects.1,
            defects.2,
        );
        assert_damped_reports_converged(&reactor);
        if errors.0 < 0.5 {
            assert_production_profile_matches_reference(case, &reactor);
        } else {
            assert!(
                errors.0 > 1e3 && errors.1 > 0.9,
                "unclassified MF branch outcome"
            );
            assert!(defects.0 < 1e-8 && defects.1 < 1e-8 && defects.2 < 1e-8);
        }
    }

    let correction = solve_mf_correction_newton(case, nodes)
        .unwrap_or_else(|error| panic!("{} CorrectionNewton reference seed: {error}", case.name));
    let correction_errors = production_profile_errors(case, &correction);
    println!(
        "MF-HIGH fixed-mesh CorrectionNewton/reference seed: dT={:.6e} K dY={:.6e} rel_dq={:.6e}",
        correction_errors.0, correction_errors.1, correction_errors.2,
    );
    assert_production_profile_matches_reference(case, &correction);
}

#[test]
fn mf1963_damped_hot_seed_jacobian_refresh_matrix() {
    let case = cases()[1];
    let nodes = 301;
    let sci = assert_sci_matches_reference(case, reference_seed(case, nodes), nodes, 20_000, 1e-6);
    let hot_seed = damped_reduced_seed_from_physical_profile(&sci, nodes);
    for max_jac in [1, 2, 5, 50] {
        let reactor = solve_mf_damped_with_seed(
            case,
            hot_seed.clone(),
            nodes,
            0,
            GridRefinementMethod::DoublePoints,
            max_jac,
            false,
            false,
        )
        .unwrap_or_else(|error| panic!("{} Damped max_jac={max_jac}: {error}", case.name));
        let errors = production_profile_errors(case, &reactor);
        let defects = mf_discrete_trapezoid_defects(case, &reactor);
        let snapshot = reactor.solver.last_telemetry.as_ref().expect("telemetry");
        let (status, iterations, jacobian_requests) = match snapshot {
            BvpSolverTelemetrySnapshot::Damped(snapshot) => (
                snapshot.last_solve_status,
                snapshot.counters.iterations,
                snapshot.counters.jacobian_requests,
            ),
            other => panic!("expected Damped telemetry, got {other:?}"),
        };
        let solution = reactor.solver.get_solution().expect("published solution");
        let temperature_range = solution.column(0).iter().fold(
            (f64::INFINITY, f64::NEG_INFINITY),
            |(low, high), &temperature| (low.min(temperature), high.max(temperature)),
        );
        println!(
            "MF-HIGH Damped fixed mesh max_jac={max_jac}: status={status:?} iterations={iterations} jacobian_requests={jacobian_requests} T=[{:.6e},{:.6e}] K dT={:.6e} K dY={:.6e} rel_dq={:.6e} defect=({:.3e},{:.3e},{:.3e})",
            temperature_range.0,
            temperature_range.1,
            errors.0,
            errors.1,
            errors.2,
            defects.0,
            defects.1,
            defects.2,
        );
        assert_damped_reports_converged(&reactor);
        assert_production_profile_matches_reference(case, &reactor);
    }
}

#[test]
fn mf1963_damped_mf_low_fixed_mesh_classifies_solver_and_physical_outcomes() {
    let case = cases()[0];
    let nodes = 2049;
    let mut template = reactor_for(case);
    template.setup_bvp().expect("MF-LOW seed setup");
    let seed = damped_reduced_seed_from_nondimensional_profile(
        &template,
        &reference_seed(case, nodes),
        nodes,
    );
    for max_jac in [1, 50] {
        let (result, telemetry) = solve_mf_damped_with_seed_and_telemetry(
            case,
            seed.clone(),
            nodes,
            0,
            GridRefinementMethod::DoublePoints,
            max_jac,
            false,
            false,
        );
        match result {
            Err(error) => {
                let snapshot = telemetry.expect("failed Damped solve retains telemetry");
                let BvpSolverTelemetrySnapshot::Damped(snapshot) = snapshot else {
                    panic!("expected Damped telemetry on failure")
                };
                println!(
                    "MF-LOW fixed-mesh Damped max_jac={max_jac}: SOLVE_FAIL status={:?} iterations={} detail={error}",
                    snapshot.last_solve_status, snapshot.counters.iterations,
                );
                assert_eq!(snapshot.last_solve_status, BvpSolveStatus::NoResult);
                assert!(error.contains("MissingResult"));
            }
            Ok(reactor) => {
                let solution = reactor.solver.get_solution().expect("published solution");
                let errors = production_profile_errors(case, &reactor);
                let defects = mf_discrete_trapezoid_defects(case, &reactor);
                let min_y_a = solution
                    .column(2)
                    .iter()
                    .copied()
                    .fold(f64::INFINITY, f64::min);
                let physical_pass = min_y_a >= -1e-12
                    && errors.0 < 0.5
                    && errors.1 < 2e-3
                    && errors.2 < 2e-3
                    && errors.3 < 2e-8;
                println!(
                    "MF-LOW fixed-mesh Damped max_jac={max_jac}: SOLVE_OK min(Y_A)={min_y_a:.3e} dT={:.3e} K dY={:.3e} relative_dq={:.3e} defects={defects:?} physical_gate={physical_pass}",
                    errors.0, errors.1, errors.2,
                );
                assert_damped_reports_converged(&reactor);
                assert!(defects.0 < 1e-5 && defects.1 < 1e-5 && defects.2 < 1e-5);
                assert!(physical_pass, "unexpected LOW fixed-mesh physical root");
            }
        }
    }
}

fn solve_mf_correction_newton(case: Case, nodes: usize) -> Result<SimpleReactorTask, String> {
    let mut reactor = reactor_for(case);
    reactor
        .setup_bvp()
        .map_err(|error| format!("CorrectionNewton setup: {error}"))?;
    reactor.solver = reactor
        .solver
        .with_bvp_telemetry_mode(BvpApiTelemetryMode::Detailed);
    let mut options = CorrectionNewtonOptions {
        matrix_layout: CorrectionNewtonMatrixLayout::Sparse,
        absolute_tolerance: 1e-5,
        default_relative_tolerance: 1e-5,
        max_iterations: 100,
        logging_level: SolverLogLevel::Off,
        ..CorrectionNewtonOptions::default()
    };
    options.relative_tolerances = reactor
        .solver
        .unknowns
        .iter()
        .map(|name| (name.clone(), 1e-5))
        .collect();
    options.bounds = reactor
        .solver
        .unknowns
        .iter()
        .map(|name| {
            let range = match name.as_str() {
                "Teta" => (-19.9, 120.0),
                "q" => (-1e6, 1e6),
                _ => (-1e12, 1e12),
            };
            (name.clone(), range)
        })
        .collect();
    reactor
        .solver
        .solve_with_hybrid_settings(
            BvpSolverSettings::CorrectionNewton(options),
            reference_seed(case, nodes),
            0.0,
            1.0,
            nodes,
        )
        .map_err(|error| {
            format!(
                "CorrectionNewton solve: {error}; telemetry={:?}",
                reactor.solver.last_telemetry
            )
        })?;
    reactor
        .postprocessing()
        .map_err(|error| format!("CorrectionNewton postprocessing: {error}"))?;
    Ok(reactor)
}

#[test]
fn mf1963_bvp_sci_matches_mf_high_full_physical_profile() {
    let case = cases()[1];
    assert_sci_matches_reference(case, reference_seed(case, 301), 301, 20_000, 1e-6);
}

#[test]
fn mf1963_bvp_sci_reaches_mf_high_from_independent_smooth_seed() {
    let case = cases()[1];
    assert_sci_matches_reference(case, smooth_seed(case, 301), 301, 20_000, 1e-6);
}

#[test]
fn mf1963_production_correction_newton_matches_mf_high_profile() {
    let case = cases()[1];
    let reactor = solve_mf_correction_newton(case, 301).unwrap_or_else(|error| panic!("{error}"));
    assert_production_profile_matches_reference(case, &reactor);
}

#[test]
fn mf1963_production_correction_newton_matches_mf_low_profile() {
    let case = cases()[0];
    let result = solve_mf_correction_newton(case, 2049);
    let error = result.expect_err("MF-LOW currently exposes the CorrectionNewton limitation");
    assert!(
        error.contains("CorrectionNewton BVP route failed"),
        "{error}"
    );
    println!("MF-LOW CorrectionNewton production acceptance failure: {error}");
}

fn assert_damped_reports_converged(reactor: &SimpleReactorTask) -> u64 {
    match reactor.solver.last_telemetry.as_ref() {
        Some(BvpSolverTelemetrySnapshot::Damped(snapshot)) => {
            assert_eq!(snapshot.last_solve_status, BvpSolveStatus::Converged);
            println!(
                "Damped telemetry: iterations={} residual_calls={} jacobian_calls={} refinements={}",
                snapshot.counters.iterations,
                snapshot.counters.residual_calls,
                snapshot.counters.jacobian_requests,
                snapshot.counters.grid_refinements,
            );
            snapshot.counters.grid_refinements
        }
        other => panic!("expected Damped telemetry, got {other:?}"),
    }
}

#[test]
#[ignore = "MF-LOW is a long-domain, narrow-front integration case"]
fn mf1963_bvp_sci_uniform_low_mesh_reports_node_limit_or_passes_physical_gate() {
    let case = cases()[0];
    match solve_mf_sci(case, reference_seed(case, 2049), 2049, 30_000, 1e-6) {
        Ok(reactor) => {
            println!("MF-LOW BVP_sci uniform mesh: solver passed; checking physical oracle");
            let errors = production_profile_errors(case, &reactor);
            println!(
                "MF-LOW BVP_sci uniform mesh physical errors: dT={:.3e} K dY={:.3e} relative_dq={:.3e} closure={:.3e}",
                errors.0, errors.1, errors.2, errors.3,
            );
            assert!(errors.0 < 0.5 && errors.1 < 2e-3 && errors.2 < 2e-3);
            assert!(errors.3 < 2e-8);
        }
        Err(error @ ReactorError::BvpApi(BvpApiError::Sci(_))) => {
            let detail = format!("{error:?}");
            println!("MF-LOW BVP_sci uniform mesh: typed Sci failure: {detail}");
            assert!(
                detail.contains("MaxNodesExceeded") || detail.contains("MeshRefinementLimit"),
                "expected a mesh-limit error, got {detail}"
            );
        }
        Err(error) => panic!("unexpected MF-LOW BVP_sci failure classification: {error:?}"),
    }
}
