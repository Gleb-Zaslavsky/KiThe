//! Six-state RST/SciPy parity for the standalone MF-HIGH numerical system.

use RustedSciThe::numerical::BVP_api::{
    BvpBoundaryCondition, BvpBoundarySide, BvpProblem, BvpSolver, BvpSolverConfig,
    BvpSolverTelemetrySnapshot,
};
use RustedSciThe::numerical::BVP_sci::new::{
    BvpSciAssembly, BvpSciExecution, BvpSciMatrixLayout, BvpSciOptions, BvpSciTelemetryMode,
};
use RustedSciThe::symbolic::symbolic_engine::Expr;
use nalgebra::DMatrix;

const LAMBDA: f64 = 0.1;
const GAMMA: f64 = 0.01;
const BETA: f64 = 0.05;
const DIFFUSIVITY: f64 = 1e-4;

pub(super) fn fixture_rows(csv: &str) -> Vec<[f64; 7]> {
    csv.lines()
        .skip(1)
        .filter(|line| !line.trim().is_empty())
        .map(|line| {
            line.split(',')
                .map(|field| field.parse::<f64>().expect("numeric six-state fixture"))
                .collect::<Vec<_>>()
                .try_into()
                .expect("fixture row: xi plus six states")
        })
        .collect()
}

pub(super) fn build_problem(seed: &[[f64; 7]], delta: f64, pe: f64) -> Result<BvpProblem, String> {
    let denominator = format!("(1+{BETA:.17e}*theta)");
    let exponential = format!("exp(theta/{denominator})");
    let source = format!("{delta:.17e}*{exponential}*C0/{denominator}");
    let equations = [
        format!("q/{LAMBDA:.17e}"),
        format!("{pe:.17e}*q-{:.17e}*{source}", LAMBDA * pe * pe),
        format!("J0/{DIFFUSIVITY:.17e}"),
        format!(
            "{pe:.17e}*J0+{:.17e}*{source}",
            DIFFUSIVITY * GAMMA * pe * pe
        ),
        format!("J1/{DIFFUSIVITY:.17e}"),
        format!(
            "{pe:.17e}*J1-{:.17e}*{source}",
            DIFFUSIVITY * GAMMA * pe * pe
        ),
    ]
    .into_iter()
    .map(|equation| Expr::parse_expression(&equation))
    .collect();
    let mesh = seed.iter().map(|row| row[0]).collect::<Vec<_>>();
    let boundaries = [
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
    BvpProblem::new(
        equations,
        vec!["theta", "q", "C0", "J0", "C1", "J1"]
            .into_iter()
            .map(str::to_owned)
            .collect(),
        "z",
        mesh.clone(),
        boundaries,
    )
    .map_err(|error| format!("standalone BvpProblem: {error}"))
}

fn solve_rst(
    seed: &[[f64; 7]],
    delta: f64,
    pe: f64,
    max_nodes: usize,
) -> Result<(Vec<[f64; 7]>, String), String> {
    let problem = build_problem(seed, delta, pe)?;
    let options = BvpSciOptions::default()
        .with_execution(BvpSciExecution::Lambdify)
        .with_assembly(BvpSciAssembly::ExprLegacy)
        .with_matrix_layout(BvpSciMatrixLayout::Sparse)
        .with_tolerance(2e-7)
        .with_boundary_tolerance(2e-7)
        .with_limits(max_nodes, 100, 20)
        .with_telemetry(BvpSciTelemetryMode::Detailed);
    let initial_guess = DMatrix::from_fn(6, seed.len(), |state, node| seed[node][state + 1]);
    let mut solver = BvpSolver::new(
        problem,
        BvpSolverConfig::Sci {
            options,
            sparse_factorization:
                RustedSciThe::numerical::BVP_api::BvpSciSparseFactorization::StandardFaer,
            aot_config: None,
            initial_guess,
        },
    )
    .map_err(|error| format!("six-state BVP_sci construction: {error}"))?;
    solver
        .solve()
        .map_err(|error| format!("six-state BVP_sci solve: {error:?}"))?;
    let view = solver.solution_view().expect("successful solution view");
    let sampled = seed
        .iter()
        .map(|row| {
            let z = row[0];
            let upper = view
                .mesh()
                .partition_point(|&value| value < z)
                .clamp(1, view.mesh().len() - 1);
            let lower = upper - 1;
            let fraction = (z - view.mesh()[lower]) / (view.mesh()[upper] - view.mesh()[lower]);
            let mut values = [z; 7];
            for state in 0..6 {
                let left = view.value(lower, state).expect("left state");
                let right = view.value(upper, state).expect("right state");
                values[state + 1] = left + fraction * (right - left);
            }
            values
        })
        .collect();
    let telemetry = match solver.telemetry_snapshot() {
        BvpSolverTelemetrySnapshot::Sci(snapshot) => format!(
            "nodes={} iterations={} refinements={} rhs={} jacobian={} solve_ms={:?}",
            view.mesh().len(),
            snapshot.newton_iterations,
            snapshot.mesh_refinements,
            snapshot.rhs_calls,
            snapshot.jacobian_calls,
            snapshot.full_solve_ms,
        ),
        other => format!("nodes={} telemetry={other:?}", view.mesh().len()),
    };
    Ok((sampled, telemetry))
}

fn assert_matches_kithe_rhs(case_index: usize) {
    let case = super::cases()[case_index];
    let mut reactor = super::reactor_for(case);
    reactor
        .setup_bvp()
        .unwrap_or_else(|error| panic!("{} KiThe setup: {error}", case.name));
    let pe = case.xi_max;
    let names = reactor
        .solver
        .unknowns
        .iter()
        .map(String::as_str)
        .collect::<Vec<_>>();
    let actual = reactor
        .solver
        .eq_system
        .iter()
        .map(|equation| equation.lambdify_borrowed_thread_safe(&names))
        .map(|evaluate| evaluate(&[20.0, 0.4, 0.73, 1.2e-5, 0.27, -1.2e-5]))
        .collect::<Vec<_>>();
    let source = case.delta * (20.0_f64 / (1.0 + BETA * 20.0)).exp() * 0.73 / (1.0 + BETA * 20.0);
    let expected = [
        0.4 / LAMBDA,
        pe * 0.4 - LAMBDA * pe * pe * source,
        1.2e-5 / DIFFUSIVITY,
        pe * 1.2e-5 + DIFFUSIVITY * GAMMA * pe * pe * source,
        -1.2e-5 / DIFFUSIVITY,
        -pe * 1.2e-5 - DIFFUSIVITY * GAMMA * pe * pe * source,
    ];
    for (index, (actual, expected)) in actual.into_iter().zip(expected).enumerate() {
        assert!(
            (actual - expected).abs() <= 1e-11 * expected.abs().max(1.0),
            "production/standalone RHS mismatch at state {index}: {actual:.17e} vs {expected:.17e}"
        );
    }
}

#[test]
fn mf1963_six_state_rst_sci_matches_standalone_scipy_fixture() {
    assert_matches_kithe_rhs(1);
    let seed = fixture_rows(include_str!("fixtures/mf1963_mf_high_six_state_seed.csv"));
    let scipy = fixture_rows(include_str!(
        "fixtures/mf1963_mf_high_six_state_scipy_solution.csv"
    ));
    assert_eq!(seed.len(), scipy.len());

    let (rst, telemetry) =
        solve_rst(&seed, 0.1, 0.5427079731331537, 20_000).expect("six-state MF-HIGH RST solve");
    let mut max_abs = [0.0_f64; 6];
    for (actual, expected) in rst.iter().zip(&scipy) {
        assert_eq!(actual[0], expected[0], "shared mesh coordinates");
        for state in 0..6 {
            max_abs[state] = max_abs[state].max((actual[state + 1] - expected[state + 1]).abs());
        }
    }
    println!("MF-04B six-state MF-HIGH: max_abs={max_abs:?} {telemetry}");
    for (state, tolerance) in [1e-7, 1e-7, 1e-7, 1e-7, 1e-7, 1e-7].into_iter().enumerate() {
        assert!(
            max_abs[state] < tolerance,
            "state {state} parity: {:?}",
            max_abs
        );
    }

    assert!(rst[0][1].abs() < 2e-7, "left theta boundary");
    assert!(rst.last().unwrap()[2].abs() < 2e-7, "right q boundary");
    assert!((rst[0][3] - 1.0).abs() < 2e-7, "left C0 boundary");
    for state in [4, 6] {
        assert!(
            rst.last().unwrap()[state].abs() < 2e-7,
            "right flux boundary"
        );
    }
    assert!(rst[0][5].abs() < 2e-7, "left C1 boundary");
}

#[test]
fn mf1963_six_state_rst_sci_matches_scipy_on_low_clustered_mesh() {
    assert_matches_kithe_rhs(0);
    let seed = fixture_rows(include_str!("fixtures/mf1963_mf_low_six_state_seed.csv"));
    let scipy = fixture_rows(include_str!(
        "fixtures/mf1963_mf_low_six_state_scipy_solution.csv"
    ));
    let (rst, telemetry) = solve_rst(&seed, 1e-4, 12017.324603855763, 6_000)
        .expect("clustered six-state MF-LOW RST solve");
    assert_eq!(rst.len(), scipy.len());

    let mut max_abs = [0.0_f64; 6];
    for (actual, expected) in rst.iter().zip(&scipy) {
        assert_eq!(actual[0], expected[0], "shared mesh coordinates");
        for state in 0..6 {
            max_abs[state] = max_abs[state].max((actual[state + 1] - expected[state + 1]).abs());
        }
    }
    println!("MF-04B six-state MF-LOW clustered: max_abs={max_abs:?} {telemetry}");
    for (state, tolerance) in [1e-7, 1e-7, 1e-7, 1e-7, 1e-7, 1e-7].into_iter().enumerate() {
        assert!(
            max_abs[state] < tolerance,
            "MF-LOW state {state} parity: {:?}",
            max_abs
        );
    }
    assert!(rst[0][1].abs() < 2e-7, "left theta boundary");
    assert!(rst.last().unwrap()[2].abs() < 2e-7, "right q boundary");
    assert!((rst[0][3] - 1.0).abs() < 2e-7, "left C0 boundary");
}

#[test]
fn mf1963_six_state_rst_sci_low_uniform_mesh_hits_node_limit() {
    let uniform = fixture_rows(include_str!(
        "fixtures/mf1963_mf_low_uniform_six_state_seed.csv"
    ));
    let error = solve_rst(&uniform, 1e-4, 12017.324603855763, 6_000)
        .expect_err("uniform six-state MF-LOW mesh should hit the node cap");
    println!(
        "MF-04B six-state MF-LOW uniform: nodes={} RST_FAIL={error}",
        uniform.len()
    );
    assert!(
        error.to_ascii_lowercase().contains("mesh")
            || error.to_ascii_lowercase().contains("node")
            || error.to_ascii_lowercase().contains("limit"),
        "expected typed mesh/node cap failure, got {error}"
    );
}
