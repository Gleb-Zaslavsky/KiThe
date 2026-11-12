//! Numerical-only RST reproducer for Damped initial-guess packing.
//!
//! This test uses no KiThe reactor, chemistry, or physical-unit API. It builds
//! a six-state symbolic BVP directly from `Expr` and the public `BVP_api`.

use RustedSciThe::numerical::BVP_Damp::NR_Damp_solver_damped::{DampedSolverOptions, SolverParams};
use RustedSciThe::numerical::BVP_Damp::telemetry::{BvpLoggingConfig, BvpLoggingMode};
use RustedSciThe::numerical::BVP_api::{
    BvpBoundaryCondition, BvpBoundarySide, BvpProblem, BvpSolver, BvpSolverConfig,
};
use RustedSciThe::symbolic::bvp::legacy::BvpSymbolicAssemblyBackend;
use RustedSciThe::symbolic::codegen::codegen_backend_selection::BackendSelectionPolicy;
use RustedSciThe::symbolic::symbolic_engine::Expr;
use nalgebra::DMatrix;
use std::collections::HashMap;

const INTERVALS: usize = 301;
const XI_MAX: f64 = 0.542_707_973_133_153_7;
const DELTA: f64 = 0.1;
const BETA: f64 = 0.05;
const GAMMA: f64 = 0.01;
const LAMBDA: f64 = 0.1;
const RHO_D: f64 = 1.0e-4;

fn hot_seed_fixture() -> Vec<[f64; 7]> {
    include_str!("fixtures/rst_damped_initial_guess_seed_301.csv")
        .lines()
        .skip(1)
        .filter(|line| !line.trim().is_empty())
        .map(|line| {
            let values = line
                .split(',')
                .map(|field| field.parse::<f64>().expect("numeric seed fixture"))
                .collect::<Vec<_>>();
            values
                .try_into()
                .expect("fixture row must contain x and six states")
        })
        .collect()
}

fn build_problem_and_seed() -> (BvpProblem, DMatrix<f64>, DMatrix<f64>, DMatrix<f64>) {
    let fixture = hot_seed_fixture();
    assert_eq!(fixture.len(), INTERVALS + 1);

    let xi = format!("{XI_MAX:.17e}");
    let delta = format!("{DELTA:.17e}");
    let beta = format!("{BETA:.17e}");
    let gamma = format!("{GAMMA:.17e}");
    let lambda = format!("{LAMBDA:.17e}");
    let rho_d = format!("{RHO_D:.17e}");
    let rate = format!("{delta}*exp(u0/(1+{beta}*u0))*u2/(1+{beta}*u0)");
    let q_source_scale = format!("({lambda}*{xi}*{xi})");
    let flux_source_scale = format!("({rho_d}*{gamma}*{xi}*{xi})");
    let equations = [
        format!("u1/({lambda})"),
        format!("{xi}*u1-{q_source_scale}*({rate})"),
        format!("u3/({rho_d})"),
        format!("{xi}*u3+{flux_source_scale}*({rate})"),
        format!("u5/({rho_d})"),
        format!("{xi}*u5-{flux_source_scale}*({rate})"),
    ]
    .into_iter()
    .map(|equation| Expr::parse_expression(&equation))
    .collect();
    let state_names = (0..6).map(|index| format!("u{index}")).collect();
    let mesh = (0..=INTERVALS)
        .map(|node| node as f64 / INTERVALS as f64)
        .collect();
    let boundary_conditions = [
        (0, BvpBoundarySide::Left, 0.0),
        (1, BvpBoundarySide::Right, 0.0),
        (2, BvpBoundarySide::Left, 1.0),
        (3, BvpBoundarySide::Right, 0.0),
        (4, BvpBoundarySide::Left, 0.0),
        (5, BvpBoundarySide::Right, 0.0),
    ]
    .into_iter()
    .map(|(state, side, value)| BvpBoundaryCondition { state, side, value })
    .collect();
    let problem = BvpProblem::new(equations, state_names, "x", mesh, boundary_conditions)
        .expect("valid standalone symbolic BVP");

    let full_seed = DMatrix::from_fn(6, fixture.len(), |state, node| {
        let [_, u0, u1, u2, u3, u4, u5] = fixture[node];
        [u0, u1, u2, u3, u4, u5][state]
    });
    let mut reduced_values = Vec::with_capacity(6 * INTERVALS);
    for node in 0..fixture.len() {
        for state in 0..6 {
            let fixed = (node == 0 && matches!(state, 0 | 2 | 4))
                || (node == INTERVALS && matches!(state, 1 | 3 | 5));
            if !fixed {
                reduced_values.push(full_seed[(state, node)]);
            }
        }
    }
    assert_eq!(reduced_values.len(), 6 * INTERVALS);
    let compact_seed = DMatrix::from_column_slice(6, INTERVALS, &reduced_values);
    let state_profile_seed = DMatrix::from_fn(6, INTERVALS, |state, interval| {
        let node = if matches!(state, 0 | 2 | 4) {
            interval + 1
        } else {
            interval
        };
        full_seed[(state, node)]
    });
    (problem, full_seed, compact_seed, state_profile_seed)
}

fn independent_trapezoid_residuals(full_seed: &DMatrix<f64>) -> [f64; 6] {
    let h = 1.0 / INTERVALS as f64;
    let mut sums = [0.0; 6];
    for interval in 0..INTERVALS {
        let node_rhs = |node: usize| {
            let u0 = full_seed[(0, node)];
            let u1 = full_seed[(1, node)];
            let u2 = full_seed[(2, node)];
            let u3 = full_seed[(3, node)];
            let u5 = full_seed[(5, node)];
            let rate = DELTA * (u0 / (1.0 + BETA * u0)).exp() * u2 / (1.0 + BETA * u0);
            [
                u1 / LAMBDA,
                XI_MAX * u1 - LAMBDA * XI_MAX.powi(2) * rate,
                u3 / RHO_D,
                XI_MAX * u3 + RHO_D * GAMMA * XI_MAX.powi(2) * rate,
                u5 / RHO_D,
                XI_MAX * u5 - RHO_D * GAMMA * XI_MAX.powi(2) * rate,
            ]
        };
        let rhs0 = node_rhs(interval);
        let rhs1 = node_rhs(interval + 1);
        for state in 0..6 {
            let defect = full_seed[(state, interval + 1)]
                - full_seed[(state, interval)]
                - 0.5 * h * (rhs0[state] + rhs1[state]);
            sums[state] += defect * defect;
        }
    }
    sums.map(f64::sqrt)
}

fn initial_residual_norm(
    trace: &RustedSciThe::numerical::BVP_Damp::telemetry::BvpDampedDecisionTrace,
) -> Option<f64> {
    use RustedSciThe::numerical::BVP_Damp::telemetry::BvpDampedDecisionEvent;

    trace.events.iter().find_map(|event| match event {
        BvpDampedDecisionEvent::IterationDecision {
            base_residual_norm, ..
        } => *base_residual_norm,
        _ => None,
    })
}

#[test]
fn rst_damped_state_major_guess_preserves_hot_root() {
    let (problem, full_seed, compact_seed, state_profile_seed) = build_problem_and_seed();
    let seed_max = (0..full_seed.ncols())
        .map(|node| full_seed[(0, node)])
        .fold(f64::NEG_INFINITY, f64::max);
    let independent_residuals = independent_trapezoid_residuals(&full_seed);
    let independent_residual = independent_residuals
        .iter()
        .map(|value| value * value)
        .sum::<f64>()
        .sqrt();
    assert!(seed_max > 90.0, "fixture must be a hot-branch seed");
    assert!(independent_residual.is_finite());

    let bounds = HashMap::from([
        ("u0".to_owned(), (-19.9, 120.0)),
        ("u1".to_owned(), (-1.0e6, 1.0e6)),
        ("u2".to_owned(), (-1.0e12, 1.0e12)),
        ("u3".to_owned(), (-1.0e12, 1.0e12)),
        ("u4".to_owned(), (-1.0e12, 1.0e12)),
        ("u5".to_owned(), (-1.0e12, 1.0e12)),
    ]);
    let relative_tolerances = (0..6).map(|index| (format!("u{index}"), 1.0e-8)).collect();
    let base_options = DampedSolverOptions::sparse_damped();
    let generated = base_options
        .generated_backend_config
        .clone()
        .with_backend_policy_override(Some(BackendSelectionPolicy::LambdifyOnly))
        .with_symbolic_assembly_backend(BvpSymbolicAssemblyBackend::ExprLegacy);
    let options = base_options
        .with_generated_backend_config(generated)
        .trapezoid_derivative()
        .with_strategy_params(Some(SolverParams {
            max_jac: Some(1),
            max_damp_iter: Some(20),
            adaptive: None,
            ..SolverParams::default()
        }))
        .with_abs_tolerance(1.0e-8)
        .with_rel_tolerance(relative_tolerances)
        .with_max_iterations(150)
        .with_bounds(bounds)
        .with_bvp_logging_config(
            BvpLoggingConfig::new(BvpLoggingMode::Detailed).with_max_events(4096),
        );
    let mut compact_solver = BvpSolver::new(
        problem.clone(),
        BvpSolverConfig::Damped {
            options: options.clone(),
            initial_guess: compact_seed,
        },
    )
    .expect("RST should prepare the compact-state standalone problem");
    let compact_result = compact_solver.solve().map(|_| ());
    let compact_trace = compact_solver
        .damped_decision_trace_snapshot()
        .expect("BVP_api exposes the Damped decision trace");

    let mut profile_solver = BvpSolver::new(
        problem,
        BvpSolverConfig::Damped {
            options,
            initial_guess: state_profile_seed,
        },
    )
    .expect("RST should prepare the state-profile standalone problem");
    let profile_result = profile_solver.solve().map(|_| ());
    let profile_trace = profile_solver
        .damped_decision_trace_snapshot()
        .expect("BVP_api exposes the Damped decision trace");

    let compact_initial_residual = initial_residual_norm(&compact_trace)
        .expect("compact-state trace contains the initial residual");
    let profile_initial_residual = initial_residual_norm(&profile_trace)
        .expect("state-major trace contains the initial residual");
    println!(
        "RST state-major guess packing regression (rechecked on 0.4.35): n={INTERVALS} states=6 seed_u0_max={seed_max:.12e} independent_full_residual_l2={independent_residual:.12e} row_residual_l2={independent_residuals:?} compact_permutation_initial_residual={compact_initial_residual:.12e} state_major_initial_residual={profile_initial_residual:.12e} compact_solve={compact_result:?} state_major_solve={profile_result:?}"
    );
    for event in &profile_trace.events {
        println!("{event:#?}");
    }
    let compact_u0_range = if let Some(view) = compact_solver.solution_view() {
        let u0_range =
            (0..view.mesh().len()).fold((f64::INFINITY, f64::NEG_INFINITY), |(low, high), node| {
                let value = view.value(node, 0).expect("u0");
                (low.min(value), high.max(value))
            });
        println!(
            "compact-permutation published u0 range=[{:.12e},{:.12e}]",
            u0_range.0, u0_range.1
        );
        Some(u0_range)
    } else {
        None
    };
    let profile_u0_range = if let Some(view) = profile_solver.solution_view() {
        let u0_range =
            (0..view.mesh().len()).fold((f64::INFINITY, f64::NEG_INFINITY), |(low, high), node| {
                let value = view.value(node, 0).expect("u0");
                (low.min(value), high.max(value))
            });
        println!(
            "state-major published u0 range=[{:.12e},{:.12e}]",
            u0_range.0, u0_range.1
        );
        Some(u0_range)
    } else {
        None
    };

    assert_eq!(compact_trace.dropped_events, 0);
    assert_eq!(profile_trace.dropped_events, 0);
    assert!(
        compact_result.is_ok(),
        "compact-permutation control solve: {compact_result:?}"
    );
    assert!(
        profile_result.is_ok(),
        "documented state-major solve: {profile_result:?}"
    );
    assert!(profile_initial_residual < 1e-8);
    assert!(compact_initial_residual > 1e3);
    assert!(profile_u0_range.expect("state-major solution").1 > 90.0);
    assert!(compact_u0_range.expect("compact-permutation solution").1 < 0.02);
}
