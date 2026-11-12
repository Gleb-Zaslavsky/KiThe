//! Direct numerical parity probes: RST BVP_api::Sci vs standalone SciPy.

use RustedSciThe::numerical::BVP_api::{
    BvpBoundaryCondition, BvpBoundarySide, BvpProblem, BvpSolver, BvpSolverConfig,
    BvpSolverTelemetrySnapshot,
};
use RustedSciThe::numerical::BVP_sci::new::{
    BvpSciAssembly, BvpSciExecution, BvpSciMatrixLayout, BvpSciOptions, BvpSciTelemetryMode,
};
use RustedSciThe::symbolic::symbolic_engine::Expr;
use nalgebra::DMatrix;

const BETA: f64 = 0.05;
const GAMMA: f64 = 0.01;

fn parse_rows(csv: &str, description: &str) -> Vec<[f64; 3]> {
    csv.lines()
        .skip(1)
        .filter(|line| !line.trim().is_empty())
        .map(|line| {
            let values = line
                .split(',')
                .map(|field| field.parse::<f64>().expect(description))
                .collect::<Vec<_>>();
            values.try_into().expect("fixture row: xi,theta,p")
        })
        .collect()
}

fn rst_solution(
    seed: &[[f64; 3]],
    delta: f64,
    max_nodes: usize,
) -> Result<(Vec<[f64; 3]>, String), String> {
    let delta = format!("{delta:.17e}");
    let beta = format!("{BETA:.17e}");
    let gamma = format!("{GAMMA:.17e}");
    let equations = [
        "p".to_owned(),
        format!("p-{delta}*exp(theta/(1+{beta}*theta))*(1-{gamma}*theta)/(1+{beta}*theta)"),
    ]
    .into_iter()
    .map(|equation| Expr::parse_expression(&equation))
    .collect();
    let mesh = seed.iter().map(|row| row[0]).collect::<Vec<_>>();
    let boundaries = vec![
        BvpBoundaryCondition {
            state: 0,
            side: BvpBoundarySide::Left,
            value: 0.0,
        },
        BvpBoundaryCondition {
            state: 1,
            side: BvpBoundarySide::Right,
            value: 0.0,
        },
    ];
    let problem = BvpProblem::new(
        equations,
        vec!["theta".into(), "p".into()],
        "xi",
        mesh.clone(),
        boundaries,
    )
    .map_err(|error| format!("BvpProblem construction: {error}"))?;
    let options = BvpSciOptions::default()
        .with_execution(BvpSciExecution::Lambdify)
        .with_assembly(BvpSciAssembly::ExprLegacy)
        .with_matrix_layout(BvpSciMatrixLayout::Sparse)
        .with_tolerance(2e-7)
        .with_boundary_tolerance(2e-7)
        .with_limits(max_nodes, 100, 20)
        .with_telemetry(BvpSciTelemetryMode::Detailed);
    let initial_guess = DMatrix::from_fn(2, seed.len(), |state, node| seed[node][state + 1]);
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
    .map_err(|error| format!("BvpSolver construction: {error}"))?;
    solver
        .solve()
        .map_err(|error| format!("RST BVP_sci solve: {error:?}"))?;
    let view = solver.solution_view().expect("successful solve has a view");
    let interpolate = |xi: f64, state: usize| {
        let output_mesh = view.mesh();
        let upper = output_mesh
            .partition_point(|&value| value < xi)
            .clamp(1, output_mesh.len() - 1);
        let lower = upper - 1;
        let fraction = (xi - output_mesh[lower]) / (output_mesh[upper] - output_mesh[lower]);
        let left = view.value(lower, state).expect("left state value");
        let right = view.value(upper, state).expect("right state value");
        left + fraction * (right - left)
    };
    let sampled = seed
        .iter()
        .map(|row| [row[0], interpolate(row[0], 0), interpolate(row[0], 1)])
        .collect();
    let telemetry = match solver.telemetry_snapshot() {
        BvpSolverTelemetrySnapshot::Sci(snapshot) => format!(
            "iterations={} refinements={} rhs={} jacobian={} solve_ms={:?}",
            snapshot.newton_iterations,
            snapshot.mesh_refinements,
            snapshot.rhs_calls,
            snapshot.jacobian_calls,
            snapshot.full_solve_ms,
        ),
        other => format!("unexpected={other:?}"),
    };
    Ok((sampled, format!("nodes={} {telemetry}", view.mesh().len())))
}

fn compare_pair(label: &str, seed_csv: &str, scipy_csv: &str, delta: f64, required: bool) {
    let seed = parse_rows(seed_csv, "numeric shared input fixture");
    let scipy = parse_rows(scipy_csv, "numeric frozen SciPy solution");
    assert_eq!(seed.len(), scipy.len());
    for (input, expected) in seed.iter().zip(&scipy) {
        assert_eq!(input[0], expected[0], "{label}: SciPy coordinates differ");
    }

    let (rst, telemetry) = match rst_solution(&seed, delta, 20_000) {
        Ok(solution) => solution,
        Err(error) if !required => {
            println!(
                "MF-04B {label}: RST_FAIL nodes={} error={error}",
                seed.len()
            );
            return;
        }
        Err(error) => panic!("MF-04B {label}: {error}"),
    };
    let mut max_theta_delta = 0.0_f64;
    let mut max_p_delta = 0.0_f64;
    for (actual, expected) in rst.iter().zip(&scipy) {
        max_theta_delta = max_theta_delta.max((actual[1] - expected[1]).abs());
        max_p_delta = max_p_delta.max((actual[2] - expected[2]).abs());
    }
    println!(
        "MF-04B {label}: max_delta=(theta={max_theta_delta:.9e},p={max_p_delta:.9e}) {telemetry}"
    );
    assert!(rst[0][1].abs() <= 2e-7, "{label}: left boundary");
    assert!(
        rst.last().unwrap()[2].abs() <= 2e-7,
        "{label}: right boundary"
    );
    assert!(max_theta_delta < 1e-7, "{label}: theta parity");
    assert!(max_p_delta < 1e-5, "{label}: p parity");
}

#[test]
fn mf1963_two_state_rst_sci_uniform_mf_low_mesh_hits_node_limit() {
    let uniform = parse_rows(
        include_str!("fixtures/mf1963_mf_low_uniform_two_state_seed.csv"),
        "numeric uniform MF-LOW fixture",
    );
    let error = rst_solution(&uniform, 1e-4, 6_000)
        .expect_err("uniform MF-LOW mesh should exhaust the matched 6000-node cap");
    println!(
        "MF-04B MF-LOW uniform mesh: initial_nodes={} cap=6000 RST_FAIL={error}",
        uniform.len(),
    );
    assert!(
        error.to_ascii_lowercase().contains("mesh")
            || error.to_ascii_lowercase().contains("node")
            || error.to_ascii_lowercase().contains("limit"),
        "expected a typed node/mesh limit failure, got: {error}"
    );
}

#[test]
fn mf1963_two_state_rst_sci_matches_scipy_on_shared_mf_high_fixture() {
    compare_pair(
        "MF-HIGH reference seed",
        include_str!("fixtures/mf1963_mf_high_two_state_seed.csv"),
        include_str!("fixtures/mf1963_mf_high_scipy_solution.csv"),
        0.1,
        true,
    );
}

#[test]
fn mf1963_two_state_rst_sci_matches_scipy_on_shared_mf_low_nonuniform_fixture() {
    compare_pair(
        "MF-LOW clustered reference seed",
        include_str!("fixtures/mf1963_mf_low_two_state_seed.csv"),
        include_str!("fixtures/mf1963_mf_low_scipy_solution.csv"),
        1e-4,
        true,
    );
}
