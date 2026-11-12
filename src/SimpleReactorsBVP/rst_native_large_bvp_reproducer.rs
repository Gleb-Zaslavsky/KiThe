//! Standalone RST BVP reproducer. This module intentionally uses only RST's
//! symbolic `Expr` and `BVP_api`; state names are renamed to abstract `u0..u25`.

use std::collections::HashMap;

use RustedSciThe::numerical::BVP_Damp::NR_Damp_solver_damped::{
    AdaptiveGridConfig, DampedSolverOptions, SolverParams,
};
use RustedSciThe::numerical::BVP_Damp::NR_Damp_solver_frozen::FrozenSolverOptions;
use RustedSciThe::numerical::BVP_Damp::grid_api::GridRefinementMethod;
use RustedSciThe::numerical::BVP_Damp::telemetry::BvpSolveStatus;
use RustedSciThe::numerical::BVP_api::{
    BvpBoundaryCondition, BvpBoundarySide, BvpProblem, BvpSciSparseFactorization, BvpSolver,
    BvpSolverConfig, BvpSolverTelemetrySnapshot,
};
use RustedSciThe::numerical::BVP_sci::new::{
    BvpSciAssembly, BvpSciExecution, BvpSciMatrixLayout, BvpSciOptions, BvpSciTelemetryMode,
};
use RustedSciThe::symbolic::bvp::legacy::BvpSymbolicAssemblyBackend;
use RustedSciThe::symbolic::codegen::codegen_backend_selection::BackendSelectionPolicy;
use RustedSciThe::symbolic::symbolic_engine::Expr;
use nalgebra::DMatrix;

const STATE_COUNT: usize = 26;
const SCI_SEED_NODES: usize = 33;
const FROZEN_NODES: usize = 65;
const DAMPED_INTERVALS: usize = 33;

// RHS vector for u' = F(s, u), in state order u0..u25. The source equations
// below are ordinary scalar expressions; the renaming routine only maps their
// identifiers to abstract state coordinates before Expr parsing.
const RHS: [&str; STATE_COUNT] = [
    "q / 0.05",
    "q * 49.10419167354888 - (-103030398415.33754 * C1 * C7 * exp(-6000.220351214819 / (300 * Teta + 1200)) + -4128633572188.5405 * C1 * C10 * exp(-19702.140967043542 / (300 * Teta + 1200)) + 2280232215.9957895 * C1 * C11 * exp(560.1145056531153 / (300 * Teta + 1200)) * (300 * Teta + 1200) ^ 0.73 + 33299161969.77031 * C1 * C9 * exp(-909.8715419773874 / (300 * Teta + 1200)) + 712690917.1313457 * C1 * C3 * exp(-1210.310319942266 / (300 * Teta + 1200)) * (300 * Teta + 1200) ^ 1.16) * 0.00000033333333333333335",
    "J0 / 0.000004844744729497463",
    "412.8184479613349 * J0 - 0.0001 * (28011.995829600797 * C1 * C7 * exp(-6000.220351214819 / (300 * Teta + 1200)) + 288944.2885029702 * C1 * C10 * exp(-19702.140967043542 / (300 * Teta + 1200)) + 423.17015412771735 * C1 * C11 * exp(560.1145056531153 / (300 * Teta + 1200)) * (300 * Teta + 1200) ^ 0.73 + 56060.951656390236 * C1 * C9 * exp(-909.8715419773874 / (300 * Teta + 1200)) + 624.6888893573005 * C1 * C3 * exp(-1210.310319942266 / (300 * Teta + 1200)) * (300 * Teta + 1200) ^ 1.16)",
    "J1 / 0.000004844744729497463",
    "412.8184479613349 * J1 - 0.0001 * (-1239.9012171458817 * C1 * C11 * exp(560.1145056531153 / (300 * Teta + 1200)) * (300 * Teta + 1200) ^ 0.73 + -164260.26626664176 * C1 * C9 * exp(-909.8715419773874 / (300 * Teta + 1200)) + -1830.3571428571438 * C1 * C3 * exp(-1210.310319942266 / (300 * Teta + 1200)) * (300 * Teta + 1200) ^ 1.16 + -82075.98618433015 * C1 * C7 * exp(-6000.220351214819 / (300 * Teta + 1200)) + -846615.4134633419 * C1 * C10 * exp(-19702.140967043542 / (300 * Teta + 1200)))",
    "J2 / 0.000004844744729497463",
    "412.8184479613349 * J2 - 0.0001 * (104440.78855307552 * C1 * C9 * exp(-909.8715419773874 / (300 * Teta + 1200)) + 1163.7868833320908 * C1 * C3 * exp(-1210.310319942266 / (300 * Teta + 1200)) * (300 * Teta + 1200) ^ 1.16 + 52185.966291128585 * C1 * C7 * exp(-6000.220351214819 / (300 * Teta + 1200)) + 538299.2697684199 * C1 * C10 * exp(-19702.140967043542 / (300 * Teta + 1200)) + 788.3602272775099 * C1 * C11 * exp(560.1145056531153 / (300 * Teta + 1200)) * (300 * Teta + 1200) ^ 0.73)",
    "J3 / 0.000004844744729497463",
    "412.8184479613349 * J3 + 0.004188137016775251 * C1 * C3 * exp(-1210.310319942266 / (300 * Teta + 1200)) * (300 * Teta + 1200) ^ 1.16",
    "J4 / 0.000004844744729497463",
    "412.8184479613349 * J4 + -0.008376274033550502 * C1 * C3 * exp(-1210.310319942266 / (300 * Teta + 1200)) * (300 * Teta + 1200) ^ 1.16",
    "J5 / 0.000004844744729497463",
    "412.8184479613349 * J5 + -0.05070442518530647 * C1 * C11 * exp(560.1145056531153 / (300 * Teta + 1200)) * (300 * Teta + 1200) ^ 0.73",
    "J6 / 0.000004844744729497463",
    "412.8184479613349 * J6 + -6.337213343186177 * C1 * C7 * exp(-6000.220351214819 / (300 * Teta + 1200))",
    "J7 / 0.000004844744729497463",
    "412.8184479613349 * J7 - 0.0001 * (-61494.109368260986 * C1 * C7 * exp(-6000.220351214819 / (300 * Teta + 1200)) + 634312.9488745617 * C1 * C10 * exp(-19702.140967043542 / (300 * Teta + 1200)))",
    "J8 / 0.000004844744729497463",
    "412.8184479613349 * J8",
    "J9 / 0.000004844744729497463",
    "412.8184479613349 * J9 + 5.965541506821332 * C1 * C9 * exp(-909.8715419773874 / (300 * Teta + 1200))",
    "J10 / 0.000004844744729497463",
    "412.8184479613349 * J10 + 61.494109368260986 * C1 * C10 * exp(-19702.140967043542 / (300 * Teta + 1200))",
    "J11 / 0.000004844744729497463",
    "412.8184479613349 * J11 - 0.0001 * (-478.67341611241034 * C1 * C11 * exp(560.1145056531153 / (300 * Teta + 1200)) * (300 * Teta + 1200) ^ 0.73 + 63413.94112538935 * C1 * C9 * exp(-909.8715419773874 / (300 * Teta + 1200)))",
];

const LEFT_VALUES: [f64; 12] = [
    0.0,
    0.018405139435376694,
    0.0,
    3.5094811081461935e-7,
    0.0,
    0.01254430598477255,
    0.0,
    0.002298292330465739,
    0.7656249768090929,
    5.570256770757038e-6,
    0.20052924374725334,
    0.0005921204881571657,
];

fn state_index(token: &str) -> Option<usize> {
    match token {
        "Teta" => Some(0),
        "q" => Some(1),
        _ => {
            let (prefix, suffix) = token.split_at(1);
            let index = suffix.parse::<usize>().ok()?;
            if index >= 12 {
                return None;
            }
            match prefix {
                "C" => Some(2 + 2 * index),
                "J" => Some(3 + 2 * index),
                _ => None,
            }
        }
    }
}

fn rename_state_tokens(source: &str) -> String {
    let bytes = source.as_bytes();
    let mut output = String::with_capacity(source.len() + 8);
    let mut cursor = 0;
    while cursor < bytes.len() {
        if bytes[cursor].is_ascii_alphabetic() || bytes[cursor] == b'_' {
            let start = cursor;
            cursor += 1;
            while cursor < bytes.len()
                && (bytes[cursor].is_ascii_alphanumeric() || bytes[cursor] == b'_')
            {
                cursor += 1;
            }
            let token = &source[start..cursor];
            if let Some(state) = state_index(token) {
                output.push_str(&format!("u{state}"));
            } else {
                output.push_str(token);
            }
        } else {
            output.push(bytes[cursor] as char);
            cursor += 1;
        }
    }
    output
}

fn equations() -> Vec<Expr> {
    RHS.iter()
        .map(|rhs| Expr::parse_expression(&rename_state_tokens(rhs)))
        .collect()
}

fn boundary_conditions() -> Vec<BvpBoundaryCondition> {
    let mut conditions = Vec::with_capacity(STATE_COUNT);
    conditions.push(BvpBoundaryCondition {
        state: 0,
        side: BvpBoundarySide::Left,
        value: 0.0,
    });
    conditions.push(BvpBoundaryCondition {
        state: 1,
        side: BvpBoundarySide::Right,
        value: 0.0,
    });
    for (index, value) in LEFT_VALUES.into_iter().enumerate() {
        conditions.push(BvpBoundaryCondition {
            state: 2 + 2 * index,
            side: BvpBoundarySide::Left,
            value,
        });
        conditions.push(BvpBoundaryCondition {
            state: 3 + 2 * index,
            side: BvpBoundarySide::Right,
            value: 0.0,
        });
    }
    conditions
}

fn mesh(nodes: usize) -> Vec<f64> {
    (0..nodes)
        .map(|index| index as f64 / (nodes - 1) as f64)
        .collect()
}

fn problem(nodes: usize) -> BvpProblem {
    let states = (0..STATE_COUNT).map(|index| format!("u{index}")).collect();
    BvpProblem::new(
        equations(),
        states,
        "s".to_string(),
        mesh(nodes),
        boundary_conditions(),
    )
    .expect("frozen numerical BVP fixture should validate")
}

fn inlet_guess(nodes: usize) -> DMatrix<f64> {
    DMatrix::from_fn(STATE_COUNT, nodes, |state, _| match state {
        0 | 1 => 0.0,
        odd if odd % 2 == 1 => 0.0,
        even => LEFT_VALUES[(even - 2) / 2],
    })
}

fn damped_state_bounds(state: usize) -> (f64, f64) {
    if state == 0 {
        (-1.0, 3.0)
    } else if state == 1 || state % 2 == 1 {
        (-1e4, 1e4)
    } else {
        (-0.1, 1.1)
    }
}

fn assert_damped_profile_within_configured_bounds(solver: &BvpSolver) {
    let view = solver
        .solution_view()
        .expect("Damped should publish a solution");
    for node in 0..view.mesh().len() {
        for state in 0..STATE_COUNT {
            let value = view.value(node, state).expect("valid state index");
            let (lower, upper) = damped_state_bounds(state);
            let slack = 32.0 * f64::EPSILON * lower.abs().max(upper.abs()).max(1.0);
            assert!(
                value >= lower - slack && value <= upper + slack,
                "Damped state u{state}={value} at node {node} escaped [{lower}, {upper}]"
            );
        }
    }
}

fn calculate_profile_metrics(solver: &BvpSolver) -> (usize, f64, f64, f64) {
    let view = solver.solution_view().expect("solver must retain a result");
    let node_count = view.mesh().len();
    let mut min_even_coordinate = f64::INFINITY;
    let mut max_abs_coordinate = 0.0_f64;
    let mut max_even_coordinate_sum_error = 0.0_f64;
    for node in 0..node_count {
        let mut even_coordinate_sum = 0.0;
        for state in 0..STATE_COUNT {
            max_abs_coordinate = max_abs_coordinate.max(view.value(node, state).unwrap().abs());
        }
        for species in 0..12 {
            let coordinate = view.value(node, 2 + 2 * species).unwrap();
            min_even_coordinate = min_even_coordinate.min(coordinate);
            even_coordinate_sum += coordinate;
        }
        max_even_coordinate_sum_error =
            max_even_coordinate_sum_error.max((even_coordinate_sum - 1.0).abs());
    }
    (
        node_count,
        min_even_coordinate,
        max_abs_coordinate,
        max_even_coordinate_sum_error,
    )
}

#[test]
#[ignore = "standalone Expr/BVP_api Sci route is intentionally expensive in Debug; run explicitly"]
fn rst_native_sci_route_accepts_the_frozen_equation_system() {
    let options = BvpSciOptions::default()
        .with_execution(BvpSciExecution::Lambdify)
        .with_assembly(BvpSciAssembly::ExprLegacy)
        .with_matrix_layout(BvpSciMatrixLayout::Dense)
        .with_tolerance(1e-5)
        .with_limits(256, 100, 20)
        .with_telemetry(BvpSciTelemetryMode::Off);
    let mut solver = BvpSolver::new(
        problem(SCI_SEED_NODES),
        BvpSolverConfig::Sci {
            options,
            sparse_factorization: BvpSciSparseFactorization::StandardFaer,
            aot_config: None,
            initial_guess: inlet_guess(SCI_SEED_NODES),
        },
    )
    .expect("RST Sci route should prepare");
    solver.solve().expect("RST Sci route should solve");
    let (nodes, min_even, max_abs, group_sum_error) = calculate_profile_metrics(&solver);
    assert_eq!(nodes, 78);
    assert!(min_even >= -1e-5);
    assert!(group_sum_error < 1e-3);
    assert!(max_abs.is_finite() && max_abs < 1e4);
    println!(
        "Sci: nodes={nodes} min_even_state={min_even:.6e} max_abs_state={max_abs:.6e} even_state_sum_error={group_sum_error:.6e}"
    );
}

#[test]
fn rst_native_frozen_dense_route_returns_a_nonphysical_profile() {
    let mut solver = BvpSolver::new(
        problem(FROZEN_NODES),
        BvpSolverConfig::Frozen {
            options: FrozenSolverOptions::dense_naive()
                .with_tolerance(1e-5)
                .with_max_iterations(500),
            initial_guess: inlet_guess(FROZEN_NODES - 1),
        },
    )
    .expect("RST Frozen route should prepare");
    solver.solve().expect("Frozen returns its current iterate");
    let (nodes, min_even, max_abs, group_sum_error) = calculate_profile_metrics(&solver);
    assert_eq!(nodes, FROZEN_NODES);
    assert!(min_even < -1e-5);
    assert!(max_abs.is_finite() && max_abs < 1e4);
    assert!(group_sum_error < 1e-3);
    println!(
        "Frozen dense_naive: nodes={nodes} min_even_state={min_even:.6e} max_abs_state={max_abs:.6e} even_state_sum_error={group_sum_error:.6e}"
    );
}

fn damped_options(adaptive: bool) -> DampedSolverOptions {
    let backend = RustedSciThe::numerical::BVP_Damp::generated_solver_handoff::GeneratedBackendConfig::sparse_defaults()
        .with_backend_policy_override(Some(BackendSelectionPolicy::LambdifyOnly))
        .with_symbolic_assembly_backend(BvpSymbolicAssemblyBackend::ExprLegacy);
    let strategy = SolverParams {
        max_jac: Some(1),
        max_damp_iter: Some(12),
        adaptive: adaptive.then_some(AdaptiveGridConfig {
            version: 1,
            max_refinements: 1,
            grid_method: GridRefinementMethod::DoublePoints,
        }),
        ..SolverParams::default()
    };
    let mut bounds = HashMap::with_capacity(STATE_COUNT);
    let mut tolerances = HashMap::with_capacity(STATE_COUNT);
    for state in 0..STATE_COUNT {
        let name = format!("u{state}");
        let bound = damped_state_bounds(state);
        bounds.insert(name.clone(), bound);
        tolerances.insert(name, 1e-5);
    }
    DampedSolverOptions::sparse_damped()
        .with_generated_backend_config(backend)
        .with_scheme_name("trapezoid")
        .with_strategy_params(Some(strategy))
        .with_abs_tolerance(1e-5)
        .with_rel_tolerance(tolerances)
        .with_max_iterations(150)
        .with_bounds(bounds)
}

#[test]
fn rst_native_damped_fixed_grid_converges_with_trapezoid() {
    let mut solver = BvpSolver::new(
        problem(DAMPED_INTERVALS + 1),
        BvpSolverConfig::Damped {
            options: damped_options(false),
            initial_guess: inlet_guess(DAMPED_INTERVALS),
        },
    )
    .expect("RST Damped route should prepare");
    let solve_ok = solver.solve().is_ok();
    let BvpSolverTelemetrySnapshot::Damped(snapshot) = solver.telemetry_snapshot() else {
        panic!("expected Damped telemetry");
    };
    println!(
        "Damped fixed: solve_ok={} status={:?} iterations={} damping_trials={} linear_solves={}",
        solve_ok,
        snapshot.last_solve_status,
        snapshot.counters.iterations,
        snapshot.counters.damping_trials,
        snapshot.counters.linear_solves
    );
    assert!(solve_ok, "fixed-grid Damped solve should converge");
    assert_eq!(snapshot.last_solve_status, BvpSolveStatus::Converged);
    assert!(snapshot.counters.iterations > 0);
    assert_damped_profile_within_configured_bounds(&solver);
    let (nodes, _min_even, max_abs, group_sum_error) = calculate_profile_metrics(&solver);
    assert_eq!(nodes, DAMPED_INTERVALS + 1);
    assert!(max_abs.is_finite());
    assert!(group_sum_error < 1e-3);
}

#[test]
fn rst_native_damped_double_points_refinement_converges() {
    let mut solver = BvpSolver::new(
        problem(DAMPED_INTERVALS + 1),
        BvpSolverConfig::Damped {
            options: damped_options(true),
            initial_guess: inlet_guess(DAMPED_INTERVALS),
        },
    )
    .expect("RST Damped route should prepare");
    let solve_ok = solver.solve().is_ok();
    let BvpSolverTelemetrySnapshot::Damped(snapshot) = solver.telemetry_snapshot() else {
        panic!("expected Damped telemetry");
    };
    println!(
        "Damped DoublePoints: solve_ok={} status={:?} iterations={} refinements={}",
        solve_ok,
        snapshot.last_solve_status,
        snapshot.counters.iterations,
        snapshot.counters.grid_refinements
    );
    assert!(solve_ok, "DoublePoints Damped solve should converge");
    assert_eq!(snapshot.last_solve_status, BvpSolveStatus::Converged);
    assert!(snapshot.counters.grid_refinements > 0);
    assert_damped_profile_within_configured_bounds(&solver);
    let (nodes, min_even, max_abs, group_sum_error) = calculate_profile_metrics(&solver);
    assert!(nodes > DAMPED_INTERVALS + 1);
    assert!(min_even >= -1e-5);
    assert!(max_abs.is_finite());
    assert!(group_sum_error < 1e-3);
}
