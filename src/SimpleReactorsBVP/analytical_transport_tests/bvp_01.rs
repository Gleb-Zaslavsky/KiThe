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
use RustedSciThe::numerical::BVP_Damp::telemetry::BvpTelemetryMode;
use RustedSciThe::numerical::BVP_api::{
    BvpApiTelemetryMode, BvpSolverSettings, BvpSolverTelemetrySnapshot,
};
use RustedSciThe::numerical::logging::SolverLogLevel;
use std::collections::HashMap;

const SPECIES: [&str; 2] = ["X0", "X1"];
const MOLAR_MASS_KG_PER_MOL: f64 = 0.060052;
const T_REFERENCE_K: f64 = 298.15;
const T_SCALE_K: f64 = 100.0;
const LENGTH_M: f64 = 0.01;
const CP_J_PER_KG_K: f64 = 1000.0;
const MASS_FLUX_KG_PER_M2_S: f64 = 0.01;
const DENSITY_KG_PER_M3: f64 = 1.0;

#[derive(Clone, Copy)]
struct TransportCase {
    pe_q: f64,
    pe_d: [f64; 2],
    q_right_w_per_m2: f64,
    j_right_kg_per_m2_s: [f64; 2],
}

impl TransportCase {
    fn moderate_equal_transport() -> Self {
        Self {
            pe_q: 1.0,
            pe_d: [1.0, 1.0],
            q_right_w_per_m2: 1000.0,
            j_right_kg_per_m2_s: [0.002, -0.002],
        }
    }

    fn unequal_transport() -> Self {
        Self {
            pe_q: 1.0,
            pe_d: [1.0, 0.5],
            q_right_w_per_m2: 1000.0,
            j_right_kg_per_m2_s: [0.002, -0.002],
        }
    }

    fn stiff_boundary_layer() -> Self {
        Self {
            pe_q: 30.0,
            pe_d: [30.0, 30.0],
            q_right_w_per_m2: 250.0,
            j_right_kg_per_m2_s: [0.0015, -0.0015],
        }
    }
}

fn analytical_transport_fixture(case: TransportCase) -> SimpleReactorTask {
    let lambda = LENGTH_M * MASS_FLUX_KG_PER_M2_S * CP_J_PER_KG_K / case.pe_q;
    let diffusion = case
        .pe_d
        .map(|pe| LENGTH_M * MASS_FLUX_KG_PER_M2_S / (DENSITY_KG_PER_M3 * pe));
    let pressure = DENSITY_KG_PER_M3 * R_G * T_REFERENCE_K / MOLAR_MASS_KG_PER_MOL;

    let mut reactor = SimpleReactorTask::new();
    reactor.set_density_model(DensityModel::Constant);
    reactor
        .fast_react_set(vec![FastElemReact {
            eq: "X0=>X1".to_string(),
            A: 0.0,
            n: 0.0,
            E: 0.0,
            Q: 0.0,
        }])
        .expect("finite zero-rate elementary reaction should be accepted");
    reactor.kindata.substances = SPECIES.iter().map(|name| (*name).to_string()).collect();
    reactor.kindata.groups = Some(HashMap::from([
        (
            "X0".to_string(),
            HashMap::from([
                ("C".to_string(), 2),
                ("H".to_string(), 4),
                ("O".to_string(), 2),
            ]),
        ),
        (
            "X1".to_string(),
            HashMap::from([
                ("C".to_string(), 2),
                ("H".to_string(), 4),
                ("O".to_string(), 2),
            ]),
        ),
    ]));

    let mut boundary_conditions = HashMap::from([
        ("T".to_string(), T_REFERENCE_K),
        ("q".to_string(), case.q_right_w_per_m2),
        ("X0".to_string(), 0.6),
        ("X1".to_string(), 0.4),
        ("X0_flux".to_string(), case.j_right_kg_per_m2_s[0]),
        ("X1_flux".to_string(), case.j_right_kg_per_m2_s[1]),
    ]);
    for species in SPECIES {
        boundary_conditions
            .entry(species.to_string())
            .or_insert(0.0);
    }

    reactor.set_parameters(
        vec![0.0],
        pressure,
        T_REFERENCE_K,
        CP_J_PER_KG_K,
        boundary_conditions,
        lambda,
        HashMap::from([
            (SPECIES[0].to_string(), diffusion[0]),
            (SPECIES[1].to_string(), diffusion[1]),
        ]),
        MASS_FLUX_KG_PER_M2_S,
        ScalingConfig::new(T_REFERENCE_K, LENGTH_M, T_SCALE_K),
    );
    // The production task contract requires the mean molar mass before setup.
    reactor.M = MOLAR_MASS_KG_PER_MOL;
    reactor
}

fn assert_close(actual: f64, expected: f64, tolerance: f64, label: &str) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{label}: actual={actual:.16e}, expected={expected:.16e}, error={:.3e}, tolerance={tolerance:.3e}",
        (actual - expected).abs()
    );
}

fn integral_f(p: f64, z: f64) -> f64 {
    if p.abs() < 1e-12 {
        z
    } else {
        ((-p * (1.0 - z)).exp() - (-p).exp()) / p
    }
}

fn analytic_physical_state(case: TransportCase, z: f64) -> [f64; 6] {
    let q_right_hat = case.q_right_w_per_m2 * LENGTH_M / T_SCALE_K;
    let j_right_hat = case.j_right_kg_per_m2_s.map(|flux| flux * LENGTH_M);
    let q_hat = q_right_hat * (case.pe_q * (z - 1.0)).exp();
    let teta = (T_REFERENCE_K - T_REFERENCE_K) / T_SCALE_K
        + (q_right_hat / (LENGTH_M * MASS_FLUX_KG_PER_M2_S * CP_J_PER_KG_K / case.pe_q))
            * integral_f(case.pe_q, z);
    let j_hat = [
        j_right_hat[0] * (case.pe_d[0] * (z - 1.0)).exp(),
        j_right_hat[1] * (case.pe_d[1] * (z - 1.0)).exp(),
    ];
    let c0 = 0.6
        + j_right_hat[0] / (DENSITY_KG_PER_M3 * LENGTH_M * MASS_FLUX_KG_PER_M2_S / case.pe_d[0])
            * integral_f(case.pe_d[0], z);
    let c1 = 0.4
        + j_right_hat[1] / (DENSITY_KG_PER_M3 * LENGTH_M * MASS_FLUX_KG_PER_M2_S / case.pe_d[1])
            * integral_f(case.pe_d[1], z);

    [
        T_REFERENCE_K + T_SCALE_K * teta,
        q_hat * T_SCALE_K / LENGTH_M,
        c0,
        j_hat[0] / LENGTH_M,
        c1,
        j_hat[1] / LENGTH_M,
    ]
}

fn basic_initial_guess(nodes_or_intervals: usize) -> nalgebra::DMatrix<f64> {
    nalgebra::DMatrix::from_fn(6, nodes_or_intervals, |state, _| match state {
        0 => 0.25,
        1 => 0.04,
        2 => 0.6,
        3 => 0.0,
        4 => 0.4,
        5 => 0.0,
        _ => unreachable!(),
    })
}

fn solve_sci(
    case: TransportCase,
    seed_nodes: usize,
    max_nodes: usize,
    tolerance: f64,
) -> SimpleReactorTask {
    let mut reactor = analytical_transport_fixture(case);
    reactor.setup_bvp().expect("fixture should prepare");
    reactor.solver.arg_name = "z".to_string();
    reactor.solver.x_range = (0.0, 1.0);
    reactor
        .solver
        .solve_BVPsci(
            basic_initial_guess(seed_nodes),
            seed_nodes,
            max_nodes,
            tolerance,
        )
        .expect("BVP_sci should solve source-free transport fixture");
    reactor
}

fn published_profile_max_error_detail(
    case: TransportCase,
    reactor: &mut SimpleReactorTask,
) -> (f64, usize, usize, f64) {
    reactor
        .postprocessing()
        .expect("dimensional publication should succeed");
    let solution = reactor.solver.get_solution().unwrap();
    let mesh = reactor.solver.x_mesh.as_ref().unwrap();
    let mut max_abs_error = 0.0_f64;
    let mut worst = (0, 0, 0.0);
    for (node, x) in mesh.iter().enumerate() {
        let expected = analytic_physical_state(case, x / LENGTH_M);
        for state in 0..6 {
            let error = (solution[(node, state)] - expected[state]).abs();
            if error > max_abs_error {
                max_abs_error = error;
                worst = (node, state, *x / LENGTH_M);
            }
        }
    }
    (max_abs_error, worst.0, worst.1, worst.2)
}

fn published_profile_max_error(case: TransportCase, reactor: &mut SimpleReactorTask) -> f64 {
    published_profile_max_error_detail(case, reactor).0
}

fn assert_published_profile_matches(
    case: TransportCase,
    reactor: &mut SimpleReactorTask,
    tolerance: f64,
) -> f64 {
    let max_abs_error = published_profile_max_error(case, reactor);
    let solution = reactor.solver.get_solution().unwrap();
    let mesh = reactor.solver.x_mesh.as_ref().unwrap();
    let mut max_outlet_layer_error = 0.0_f64;
    let mut worst = (0usize, 0usize, 0.0_f64, 0.0_f64);
    let lambda = LENGTH_M * MASS_FLUX_KG_PER_M2_S * CP_J_PER_KG_K / case.pe_q;
    let mut max_thermal_first_integral_error = 0.0_f64;
    let mut max_species_first_integral_error = 0.0_f64;
    for (node, x) in mesh.iter().enumerate() {
        let expected = analytic_physical_state(case, x / LENGTH_M);
        for state in 0..6 {
            let error = (solution[(node, state)] - expected[state]).abs();
            if error > worst.2 {
                worst = (node, state, error, expected[state]);
            }
            if x / LENGTH_M >= 0.8 {
                max_outlet_layer_error = max_outlet_layer_error.max(error);
            }
        }
        max_thermal_first_integral_error = max_thermal_first_integral_error.max(
            (solution[(node, 1)]
                - solution[(0, 1)]
                - lambda * case.pe_q / LENGTH_M * (solution[(node, 0)] - solution[(0, 0)]))
                .abs(),
        );
        for (fraction_state, flux_state) in [(2usize, 3usize), (4, 5)] {
            max_species_first_integral_error = max_species_first_integral_error.max(
                (solution[(node, flux_state)]
                    - solution[(0, flux_state)]
                    - MASS_FLUX_KG_PER_M2_S
                        * (solution[(node, fraction_state)] - solution[(0, fraction_state)]))
                    .abs(),
            );
        }
    }
    assert!(
        max_abs_error <= tolerance,
        "profile max error {max_abs_error:.3e} exceeds {tolerance:.3e}; outlet-layer max={max_outlet_layer_error:.3e}, nodes={}, worst node/state/value/expected={worst:?}",
        mesh.len(),
    );
    assert_close(mesh[0], 0.0, 1e-14, "published x start");
    assert_close(mesh[mesh.len() - 1], LENGTH_M, 1e-14, "published x end");
    println!(
        "BVP-01 first integrals: thermal_abs={max_thermal_first_integral_error:.6e} W/m2 species_abs={max_species_first_integral_error:.6e} kg/(m2 s)"
    );
    max_abs_error
}

#[test]
fn bvp01_zero_rate_setup_has_expected_mass_transport_and_boundary_scaling() {
    let mut reactor = analytical_transport_fixture(TransportCase {
        pe_q: 1.0,
        pe_d: [1.0, 1.0],
        q_right_w_per_m2: 0.0,
        j_right_kg_per_m2_s: [0.0, 0.0],
    });
    reactor
        .setup_bvp()
        .expect("zero-rate transport fixture should complete production setup");

    assert_eq!(reactor.kindata.substances, ["X0", "X1"]);
    assert_eq!(
        reactor.solver.unknowns,
        ["Teta", "q", "C0", "J0", "C1", "J1"]
    );
    assert_close(
        reactor.mean_molar_mass_kg_per_mol(),
        MOLAR_MASS_KG_PER_MOL,
        2e-7,
        "configured mean molar mass consistent with X0/X1 formula masses",
    );
    assert_close(
        reactor
            .ideal_gas_density()
            .expect("valid ideal-gas density"),
        DENSITY_KG_PER_M3,
        2e-7,
        "reference density",
    );
    assert_close(reactor.thermal_peclet_number(), 1.0, 2e-7, "thermal Peclet");
    for (index, pe) in reactor.mass_peclet_numbers().iter().enumerate() {
        assert_close(*pe, 1.0, 2e-7, &format!("mass Peclet[{index}]"));
        assert_close(
            *reactor
                .transport_cache()
                .get(SPECIES[index])
                .expect("transport coefficient should be cached after setup"),
            1e-4,
            2e-11,
            &format!("K[{index}] = rho*D"),
        );
    }

    let molar_masses = reactor
        .kindata
        .stecheodata
        .vec_of_molmasses
        .as_ref()
        .expect("kinetic setup should prepare molar masses");
    assert_eq!(molar_masses.len(), 2);
    for (index, mass_g_per_mol) in molar_masses.iter().enumerate() {
        assert_close(
            *mass_g_per_mol / 1000.0,
            MOLAR_MASS_KG_PER_MOL,
            2e-7,
            &format!("species molar mass[{index}]"),
        );
    }

    let variables: Vec<&str> = reactor.solver.unknowns.iter().map(String::as_str).collect();
    let state = [0.0, 0.17, 0.6, 0.002, 0.4, -0.003];
    let reaction_rate = reactor
        .map_eq_rate
        .get("X0=>X1")
        .expect("reaction rate should be part of assembled system")
        .eval_expression(&variables, &state);
    assert_eq!(
        reaction_rate, 0.0,
        "zero pre-exponential must make rate zero"
    );
    assert_eq!(
        reactor.heat_release.eval_expression(&variables, &state),
        0.0,
        "zero heat/rate must make the thermal source exactly zero"
    );

    let rhs: Vec<f64> = reactor
        .solver
        .eq_system
        .iter()
        .map(|equation| equation.eval_expression(&variables, &state))
        .collect();
    assert_close(
        rhs[1],
        reactor.thermal_peclet_number() * state[1],
        1e-14,
        "q RHS without source",
    );
    for (species_index, flux_index) in [(2usize, 3usize), (4, 5)] {
        assert_close(
            rhs[flux_index],
            reactor.mass_peclet_numbers()[if species_index == 2 { 0 } else { 1 }]
                * state[flux_index],
            1e-14,
            "species flux RHS without source",
        );
    }

    let expected_bc = [
        ("Teta", (0, 0.0)),
        ("q", (1, 0.0)),
        ("C0", (0, 0.6)),
        ("J0", (1, 0.0)),
        ("C1", (0, 0.4)),
        ("J1", (1, 0.0)),
    ];
    for (variable, (side, value)) in expected_bc {
        let actual = reactor
            .solver
            .BorderConditions
            .get(variable)
            .expect("all expected boundary states are present");
        assert_eq!(actual.0, side, "boundary side for {variable}");
        assert_close(actual.1, value, 0.0, &format!("solver-space BC {variable}"));
    }
}

#[test]
fn bvp01_analytic_transport_reference_handles_zero_and_high_peclet_limits() {
    // This test is deliberately independent of the production BVP equations.
    assert_close(integral_f(0.0, 0.37), 0.37, 1e-15, "F(0,z) limit");
    assert_close(integral_f(1.0, 1.0), 1.0 - (-1.0f64).exp(), 1e-15, "F(1,1)");
    assert_close(
        integral_f(30.0, 1.0),
        (1.0 - (-30.0f64).exp()) / 30.0,
        1e-15,
        "F(30,1)",
    );
}

fn solve_damped_with_strategy(
    case: TransportCase,
    intervals: usize,
    strategy_params: Option<SolverParams>,
    telemetry: Option<BvpApiTelemetryMode>,
) -> SimpleReactorTask {
    let mut reactor = analytical_transport_fixture(case);
    reactor.setup_bvp().expect("Damped fixture should prepare");
    reactor.solver.arg_name = "z".to_string();
    reactor.solver.x_range = (0.0, 1.0);
    if let Some(mode) = telemetry {
        reactor.solver = reactor.solver.with_bvp_telemetry_mode(mode);
    }
    let names = ["Teta", "q", "C0", "J0", "C1", "J1"];
    let config = NrbvpHandoffConfig::new(
        basic_initial_guess(intervals),
        0.0,
        1.0,
        intervals,
        "forward".to_string(),
        "Damped".to_string(),
        strategy_params,
        None,
        "Banded".to_string(),
        1e-8,
        Some(
            names
                .into_iter()
                .map(|name| (name.to_string(), 1e-8))
                .collect(),
        ),
        100,
        Some(HashMap::from([
            ("Teta".to_string(), (-10.0, 10.0)),
            ("q".to_string(), (-10.0, 10.0)),
            ("C0".to_string(), (-1.0, 2.0)),
            ("J0".to_string(), (-10.0, 10.0)),
            ("C1".to_string(), (-1.0, 2.0)),
            ("J1".to_string(), (-10.0, 10.0)),
        ])),
        Some("off".to_string()),
    )
    .with_solver_backend_config(ReactorBvpSolverConfig {
        symbolic_backend: ReactorBvpSymbolicBackend::ExprLegacy,
        ..ReactorBvpSolverConfig::default_lambdify()
    });
    reactor
        .solver
        .solve_with_bvp_api(config)
        .expect("Damped should solve source-free transport fixture");
    reactor
}

fn solve_damped(case: TransportCase, intervals: usize) -> SimpleReactorTask {
    solve_damped_with_strategy(case, intervals, None, None)
}

fn solve_frozen(case: TransportCase, nodes: usize) -> SimpleReactorTask {
    let mut reactor = analytical_transport_fixture(case);
    reactor.setup_bvp().expect("Frozen fixture should prepare");
    reactor.solver.arg_name = "z".to_string();
    reactor.solver.x_range = (0.0, 1.0);
    reactor
        .solver
        .solve_BVPfrozen(
            basic_initial_guess(nodes),
            FrozenSolverOptions::default()
                .with_tolerance(1e-8)
                .with_max_iterations(500),
        )
        .expect("Frozen should solve source-free transport fixture");
    reactor
}

fn solve_correction_newton(case: TransportCase, intervals: usize) -> SimpleReactorTask {
    let mut reactor = analytical_transport_fixture(case);
    reactor
        .setup_bvp()
        .expect("CorrectionNewton fixture should prepare");
    reactor.solver.arg_name = "z".to_string();
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
    options.bounds = HashMap::from([
        ("Teta".to_string(), (-10.0, 10.0)),
        ("q".to_string(), (-10.0, 10.0)),
        ("C0".to_string(), (-1.0, 2.0)),
        ("J0".to_string(), (-10.0, 10.0)),
        ("C1".to_string(), (-1.0, 2.0)),
        ("J1".to_string(), (-10.0, 10.0)),
    ]);
    reactor
        .solver
        .solve_with_hybrid_settings(
            BvpSolverSettings::CorrectionNewton(options),
            basic_initial_guess(intervals),
            0.0,
            1.0,
            intervals,
        )
        .expect("CorrectionNewton should solve source-free transport fixture");
    reactor
}

#[test]
fn bvp01_explicit_zero_fluxes_publish_constant_profiles() {
    let case = TransportCase {
        pe_q: 1.0,
        pe_d: [1.0, 1.0],
        q_right_w_per_m2: 0.0,
        j_right_kg_per_m2_s: [0.0, 0.0],
    };
    let mut reactor = solve_sci(case, 17, 128, 1e-8);
    reactor
        .postprocessing()
        .expect("constant-profile dimensional publication should succeed");
    let solution = reactor.solver.get_solution().unwrap();
    for (node, row) in solution.row_iter().enumerate() {
        assert_close(row[0], T_REFERENCE_K, 1e-7, &format!("T at node {node}"));
        assert_close(row[1], 0.0, 1e-9, &format!("q at node {node}"));
        assert_close(row[2], 0.6, 1e-9, &format!("X0 at node {node}"));
        assert_close(row[3], 0.0, 1e-9, &format!("J0 at node {node}"));
        assert_close(row[4], 0.4, 1e-9, &format!("X1 at node {node}"));
        assert_close(row[5], 0.0, 1e-9, &format!("J1 at node {node}"));
    }
}

#[test]
fn bvp01_sci_moderate_transport_matches_analytic_profile_and_dimensional_output() {
    let case = TransportCase::moderate_equal_transport();
    let mut reactor = solve_sci(case, 25, 256, 1e-8);

    let scaled_solution = reactor
        .solver
        .get_solution()
        .expect("solver should publish a solution")
        .clone();
    let scaled_mesh = reactor
        .solver
        .x_mesh
        .as_ref()
        .expect("solver should publish a mesh")
        .clone();
    assert_eq!(scaled_solution.ncols(), 6);
    assert_eq!(scaled_solution.nrows(), scaled_mesh.len());
    assert!(scaled_solution.iter().all(|value| value.is_finite()));

    let mut max_scaled_profile_error = 0.0_f64;
    for (node, z) in scaled_mesh.iter().enumerate() {
        let expected = analytic_physical_state(case, *z);
        for state in 0..6 {
            let scaled_expected = match state {
                0 => (expected[0] - T_REFERENCE_K) / T_SCALE_K,
                1 => expected[1] * LENGTH_M / T_SCALE_K,
                3 | 5 => expected[state] * LENGTH_M,
                _ => expected[state],
            };
            max_scaled_profile_error = max_scaled_profile_error
                .max((scaled_solution[(node, state)] - scaled_expected).abs());
        }
    }
    assert!(
        max_scaled_profile_error < 2e-5,
        "scaled profile max error={max_scaled_profile_error:.3e}"
    );

    reactor
        .postprocessing()
        .expect("physical publication should succeed");
    let physical_solution = reactor.solver.get_solution().unwrap();
    let physical_mesh = reactor.solver.x_mesh.as_ref().unwrap();
    assert_close(physical_mesh[0], 0.0, 1e-14, "published x start");
    assert_close(
        physical_mesh[physical_mesh.len() - 1],
        LENGTH_M,
        1e-14,
        "published x end",
    );
    let mut max_physical_profile_error = 0.0_f64;
    for (node, x) in physical_mesh.iter().enumerate() {
        let expected = analytic_physical_state(case, x / LENGTH_M);
        for state in 0..6 {
            max_physical_profile_error = max_physical_profile_error
                .max((physical_solution[(node, state)] - expected[state]).abs());
        }
    }
    assert!(
        max_physical_profile_error < 2e-3,
        "dimensional profile max error={max_physical_profile_error:.3e}"
    );
    assert_close(
        physical_solution[(physical_solution.nrows() - 1, 0)],
        361.3620558829,
        2e-3,
        "outlet temperature",
    );
    assert_close(
        physical_solution[(physical_solution.nrows() - 1, 2)],
        0.7264241118,
        2e-5,
        "outlet X0 mass fraction",
    );
    assert_close(
        physical_solution[(physical_solution.nrows() - 1, 4)],
        0.2735758882,
        2e-5,
        "outlet X1 mass fraction",
    );
}

#[test]
fn bvp01_sci_cases_match_reference_including_unequal_transport_and_boundary_layer() {
    let cases = [
        (
            "equal Pe=1",
            TransportCase::moderate_equal_transport(),
            2e-3,
        ),
        ("unequal Pe", TransportCase::unequal_transport(), 2e-3),
        ("Pe=30 layer", TransportCase::stiff_boundary_layer(), 2e-2),
    ];

    for (name, case, tolerance) in cases {
        let (seed_nodes, max_nodes, solver_tolerance) = if case.pe_q > 10.0 {
            (17, 128, 1e-5)
        } else {
            (33, 256, 1e-8)
        };
        let mut reactor = solve_sci(case, seed_nodes, max_nodes, solver_tolerance);
        assert_published_profile_matches(case, &mut reactor, tolerance);
        let solution = reactor.solver.get_solution().unwrap();
        let outlet = solution.nrows() - 1;
        assert_close(
            solution[(outlet, 0)],
            analytic_physical_state(case, 1.0)[0],
            tolerance,
            &format!("{name}: outlet temperature"),
        );
        assert_close(
            solution[(outlet, 2)],
            analytic_physical_state(case, 1.0)[2],
            tolerance,
            &format!("{name}: outlet X0"),
        );
        assert_close(
            solution[(outlet, 4)],
            analytic_physical_state(case, 1.0)[4],
            tolerance,
            &format!("{name}: outlet X1"),
        );
        if name == "unequal Pe" {
            assert_close(
                solution[(outlet, 2)] + solution[(outlet, 4)],
                1.0477302437,
                2e-5,
                "unequal-diffusivity outlet fraction sum is not constrained to one",
            );
        }
        println!(
            "BVP-01 {name}: nodes={} temperature_out={:.9} X0_out={:.9} X1_out={:.9}",
            reactor.solver.x_mesh.as_ref().unwrap().len(),
            solution[(outlet, 0)],
            solution[(outlet, 2)],
            solution[(outlet, 4)]
        );
    }
}

#[test]
fn bvp01_damped_and_frozen_match_the_same_moderate_analytic_case() {
    let case = TransportCase::moderate_equal_transport();
    let intervals = 256;
    let state_count = 6;

    let mut damped = solve_damped(case, intervals);
    assert_eq!(damped.solver.get_solution().unwrap().ncols(), state_count);
    let damped_error = assert_published_profile_matches(case, &mut damped, 1.0);
    println!("BVP-01 Damped Pe=1 grid={intervals}: max physical profile error={damped_error:.6e}");

    let mut frozen = analytical_transport_fixture(case);
    frozen.setup_bvp().expect("Frozen fixture should prepare");
    frozen.solver.arg_name = "z".to_string();
    frozen.solver.x_range = (0.0, 1.0);
    frozen
        .solver
        .solve_BVPfrozen(
            basic_initial_guess(intervals + 1),
            FrozenSolverOptions::default()
                .with_tolerance(1e-8)
                .with_max_iterations(500),
        )
        .expect("Frozen should solve the same moderate transport case");
    assert_eq!(frozen.solver.get_solution().unwrap().ncols(), state_count);
    let frozen_error = assert_published_profile_matches(case, &mut frozen, 1.0);
    println!("BVP-01 Frozen Pe=1 grid={intervals}: max physical profile error={frozen_error:.6e}");
}

#[test]
fn bvp01_correction_newton_matches_analytic_profile_and_retains_telemetry() {
    let case = TransportCase::moderate_equal_transport();
    let mut coarse = solve_correction_newton(case, 128);
    let coarse_error = published_profile_max_error(case, &mut coarse);
    let intervals = 512;
    let mut reactor = solve_correction_newton(case, intervals);
    let error = published_profile_max_error(case, &mut reactor);
    assert!(
        error < coarse_error,
        "fixed-mesh refinement should reduce CorrectionNewton profile error: coarse={coarse_error:.3e}, fine={error:.3e}"
    );
    assert!(
        error < 2e-4,
        "CorrectionNewton physical profile error={error:.3e}"
    );

    let telemetry = reactor
        .solver
        .last_telemetry
        .as_ref()
        .expect("KiThe facade should retain CorrectionNewton telemetry");
    let BvpSolverTelemetrySnapshot::CorrectionNewton(Some(telemetry)) = telemetry else {
        panic!("KiThe facade should preserve the CorrectionNewton telemetry variant");
    };
    assert!(telemetry.residual_evaluations > 0);
    assert!(telemetry.jacobian_evaluations > 0);
    assert!(telemetry.linear_solves > 0);
    assert!(telemetry.callback_time > std::time::Duration::ZERO);
    println!(
        "BVP-01 CorrectionNewton Pe=1 coarse_intervals=128 coarse_error={coarse_error:.6e} intervals={intervals} profile_error={error:.6e} residuals={} jacobians={} linear_solves={} callback_ms={:.3}",
        telemetry.residual_evaluations,
        telemetry.jacobian_evaluations,
        telemetry.linear_solves,
        telemetry.callback_time.as_secs_f64() * 1e3
    );
}

#[test]
fn bvp01_damped_pe30_mesh_refinement_reduces_analytic_profile_error() {
    let case = TransportCase::stiff_boundary_layer();
    let mut errors = Vec::new();
    for intervals in [32, 64, 128, 256, 512] {
        let mut reactor = solve_damped(case, intervals);
        let error = published_profile_max_error(case, &mut reactor);
        println!(
            "BVP-01 Damped Pe=30 intervals={intervals}: max physical profile error={error:.6e}"
        );
        errors.push(error);
    }
    assert!(
        errors[1] < errors[0],
        "64-interval error should improve over 32: {errors:?}"
    );
    assert!(
        errors[2] < errors[1],
        "128-interval error should improve over 64: {errors:?}"
    );

    let mut frozen_errors = Vec::new();
    for intervals in [32, 64, 128, 256, 512] {
        let mut reactor = solve_frozen(case, intervals + 1);
        let error = published_profile_max_error(case, &mut reactor);
        println!(
            "BVP-01 Frozen Pe=30 intervals={intervals}: max physical profile error={error:.6e}"
        );
        frozen_errors.push(error);
    }
    assert!(
        frozen_errors[1] < frozen_errors[0],
        "Frozen 64-interval error should improve over 32: {frozen_errors:?}"
    );
    assert!(
        frozen_errors[2] < frozen_errors[1],
        "Frozen 128-interval error should improve over 64: {frozen_errors:?}"
    );
    for (family, errors) in [("Damped", errors), ("Frozen", frozen_errors)] {
        assert!(
            errors[3] < errors[2] && errors[4] < errors[3],
            "{family} errors should keep decreasing through 512 intervals: {errors:?}"
        );
    }
}

#[test]
fn bvp01_damped_pe30_bvp_api_refinement_strategy_matrix() {
    let case = TransportCase::stiff_boundary_layer();
    let seed_intervals = 32;
    let max_refinements = 2;
    let mut fixed_grid = solve_damped(case, seed_intervals);
    let fixed_grid_error = published_profile_max_error(case, &mut fixed_grid);
    let strategies = [
        ("DoublePoints", GridRefinementMethod::DoublePoints, true),
        ("Pearson", GridRefinementMethod::Pearson(1e-3, 1.4), false),
        (
            "GrcarSmooke",
            GridRefinementMethod::GrcarSmooke(1e-3, 1e-2, 1.4),
            false,
        ),
        (
            "TwoPoint",
            GridRefinementMethod::TwoPoint(0.02, 0.02, 2.5),
            false,
        ),
    ];

    for (name, grid_method, expect_error_reduction) in strategies {
        let strategy = SolverParams {
            adaptive: Some(AdaptiveGridConfig {
                version: 1,
                max_refinements,
                grid_method,
            }),
            ..SolverParams::default()
        };
        let mut reactor = solve_damped_with_strategy(
            case,
            seed_intervals,
            Some(strategy),
            Some(BvpApiTelemetryMode::Counters),
        );

        let telemetry = reactor
            .solver
            .last_telemetry
            .as_ref()
            .expect("BVP_api should retain Damped telemetry");
        let BvpSolverTelemetrySnapshot::Damped(telemetry) = telemetry else {
            panic!("BVP_api should return Damped telemetry for Damped route");
        };
        assert_eq!(telemetry.telemetry_mode, BvpTelemetryMode::Counters);
        let refinement_count = telemetry.counters.grid_refinements;
        assert!(
            refinement_count > 0,
            "{name} should report an actual mesh refinement"
        );
        assert!(refinement_count <= max_refinements as u64);

        let final_nodes = reactor
            .solver
            .x_mesh
            .as_ref()
            .expect("refined mesh published")
            .len();
        assert!(
            final_nodes > seed_intervals + 1,
            "{name} should add nodes: seed={} final={}",
            seed_intervals + 1,
            final_nodes
        );
        let error = published_profile_max_error(case, &mut reactor);
        println!(
            "BVP-01 Damped adaptive Pe=30 strategy={name} seed_nodes={} final_nodes={} refinements={} fixed_seed_error={fixed_grid_error:.6e} adaptive_error={error:.6e}",
            seed_intervals + 1,
            final_nodes,
            refinement_count
        );
        assert!(error.is_finite());
        if expect_error_reduction {
            assert!(
                error < fixed_grid_error,
                "{name} should improve the fixed-seed profile error: adaptive={error:.6e}, fixed={fixed_grid_error:.6e}"
            );
        } else {
            assert!(
                error <= fixed_grid_error * (1.0 + 1e-8),
                "{name} should not worsen the fixed-seed profile error: adaptive={error:.6e}, fixed={fixed_grid_error:.6e}"
            );
        }
    }
}

#[test]
fn bvp01_damped_pe30_refinement_parameter_sweep_reports_physical_accuracy() {
    let case = TransportCase::stiff_boundary_layer();
    let seed_intervals = 32;
    let candidates = [
        ("Pearson-1e-2", GridRefinementMethod::Pearson(1e-2, 1.4)),
        ("Pearson-5e-2", GridRefinementMethod::Pearson(5e-2, 1.4)),
        ("Pearson-1e-1", GridRefinementMethod::Pearson(1e-1, 1.4)),
        (
            "Grcar-1e-2-5e-2",
            GridRefinementMethod::GrcarSmooke(1e-2, 5e-2, 1.4),
        ),
        (
            "Grcar-5e-2-1e-1",
            GridRefinementMethod::GrcarSmooke(5e-2, 1e-1, 1.4),
        ),
        (
            "Grcar-1e-1-2e-1",
            GridRefinementMethod::GrcarSmooke(1e-1, 2e-1, 1.4),
        ),
        (
            "TwoPoint-5e-2-5e-2",
            GridRefinementMethod::TwoPoint(5e-2, 5e-2, 2.5),
        ),
        (
            "TwoPoint-1e-1-1e-1",
            GridRefinementMethod::TwoPoint(1e-1, 1e-1, 2.5),
        ),
        (
            "TwoPoint-2e-2-2e-2-C1.4",
            GridRefinementMethod::TwoPoint(2e-2, 2e-2, 1.4),
        ),
        (
            "TwoPoint-2e-1-2e-1",
            GridRefinementMethod::TwoPoint(2e-1, 2e-1, 2.5),
        ),
        (
            "TwoPoint-5e-1-5e-1",
            GridRefinementMethod::TwoPoint(5e-1, 5e-1, 2.5),
        ),
        (
            "TwoPoint-1e-1-1e-1-C1.4",
            GridRefinementMethod::TwoPoint(1e-1, 1e-1, 1.4),
        ),
    ];
    let mut rows = Vec::new();

    for (name, grid_method) in candidates {
        let params = SolverParams {
            adaptive: Some(AdaptiveGridConfig {
                version: 1,
                max_refinements: 4,
                grid_method,
            }),
            ..SolverParams::default()
        };
        let mut reactor = solve_damped_with_strategy(
            case,
            seed_intervals,
            Some(params),
            Some(BvpApiTelemetryMode::Counters),
        );
        let telemetry = reactor
            .solver
            .last_telemetry
            .as_ref()
            .expect("Damped parameter sweep should retain native telemetry");
        let BvpSolverTelemetrySnapshot::Damped(telemetry) = telemetry else {
            panic!("parameter sweep must use the production Damped route");
        };
        let refinements = telemetry.counters.grid_refinements;
        let (error, worst_node, worst_state, worst_z) =
            published_profile_max_error_detail(case, &mut reactor);
        let nodes = reactor.solver.x_mesh.as_ref().unwrap().len();
        assert!(
            error.is_finite(),
            "{name} produced a non-finite oracle error"
        );
        if refinements == 0 {
            assert_eq!(
                nodes,
                seed_intervals + 1,
                "{name} changed nodes without reporting a refinement"
            );
        } else {
            assert!(nodes > seed_intervals + 1, "{name} did not add mesh nodes");
        }
        rows.push((
            name,
            nodes,
            refinements,
            error,
            worst_node,
            worst_state,
            worst_z,
        ));
    }

    rows.sort_by(|left, right| left.3.total_cmp(&right.3));
    for (name, nodes, refinements, error, worst_node, worst_state, worst_z) in &rows {
        println!(
            "BVP-01 Pe=30 tuned refinement={name} nodes={nodes} refinements={refinements} max_physical_error={error:.6e} worst=(node={worst_node}, state={worst_state}, z={worst_z:.6e})"
        );
    }
    let fixed_seed_error = {
        let mut fixed = solve_damped(case, seed_intervals);
        published_profile_max_error(case, &mut fixed)
    };
    println!(
        "BVP-01 Pe=30 fixed seed nodes={} max_physical_error={fixed_seed_error:.6e}",
        seed_intervals + 1
    );
    assert!(
        rows.iter().any(|row| row.3 < fixed_seed_error),
        "no tuned strategy improved the fixed-seed error: candidates={rows:?}, seed={fixed_seed_error:.6e}"
    );
    let tuned_two_point = rows
        .iter()
        .find(|row| row.0 == "TwoPoint-2e-2-2e-2-C1.4")
        .expect("the selected tuned TwoPoint configuration must be measured");
    assert_eq!(tuned_two_point.1, 72);
    assert_eq!(tuned_two_point.2, 4);
    assert!(
        tuned_two_point.3 < 2.7 && tuned_two_point.3 < fixed_seed_error * 0.1,
        "tuned TwoPoint must materially reduce the physical error: tuned={tuned_two_point:?}, seed={fixed_seed_error:.6e}"
    );
}
