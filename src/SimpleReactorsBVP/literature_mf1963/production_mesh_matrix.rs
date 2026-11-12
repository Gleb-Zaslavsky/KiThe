//! MF-LOW solver stories through the real SimpleReactorsBVP production path.
//!
//! The RST-only matrix is kept in `rst_solver_selection`. This module checks
//! that the same clustered input survives KiThe equation construction,
//! explicit-mesh handoff, validation, and dimensional postprocessing.

use super::{
    Case, ReactorError, SimpleReactorTask, cases, interpolate, mf_discrete_trapezoid_defects,
    production_profile_errors, reactor_for, six_state_parity::fixture_rows,
};
use RustedSciThe::numerical::BVP_CorrectionNewton::{
    CorrectionNewtonMatrixLayout, CorrectionNewtonOptions,
};
use RustedSciThe::numerical::BVP_Damp::NR_Damp_solver_damped::{
    AdaptiveGridConfig, DampedSolverOptions, SolverParams,
};
use RustedSciThe::numerical::BVP_Damp::NR_Damp_solver_frozen::FrozenSolverOptions;
use RustedSciThe::numerical::BVP_Damp::grid_api::GridRefinementMethod;
use RustedSciThe::numerical::BVP_Damp::telemetry::{
    BvpLoggingConfig, BvpLoggingMode, BvpTelemetryMode,
};
use RustedSciThe::numerical::BVP_api::{
    BvpApiTelemetryMode, BvpBoundaryCondition, BvpBoundarySide, BvpProblem,
    BvpSciSparseFactorization, BvpSolver, BvpSolverSettings,
};
use RustedSciThe::numerical::BVP_sci::new::{
    BvpSciAssembly, BvpSciExecution, BvpSciMatrixLayout, BvpSciOptions, BvpSciTelemetryMode,
};
use RustedSciThe::numerical::logging::SolverLogLevel;
use RustedSciThe::symbolic::symbolic_engine::Expr;
use nalgebra::DMatrix;
use std::collections::HashMap;

const LOW_SEED: &str = include_str!("fixtures/mf1963_mf_low_six_state_seed.csv");
const RST_LOW_SCIPY_REFERENCE: &str =
    include_str!("fixtures/mf1963_mf_low_six_state_scipy_solution.csv");
const LOW_PE: f64 = 12017.324603855763;
const STATE_NAMES: [&str; 6] = ["Teta", "q", "C0", "J0", "C1", "J1"];

#[cfg(test)]
mod high_parity;
#[cfg(test)]
mod low_damped_grid_coefficients;
#[cfg(test)]
mod low_damped_tuning;

fn positive_teta_lower_bound() -> f64 {
    let zero_kelvin_teta: f64 = -300.0 / 15.0;
    zero_kelvin_teta + 64.0 * f64::EPSILON * zero_kelvin_teta.abs().max(1.0)
}

fn standalone_equations_for(case: Case) -> Vec<Expr> {
    let pe = case.xi_max;
    let source = format!(
        "{:.17e}*exp(Teta/(1+0.05*Teta))*C0/(1+0.05*Teta)",
        case.delta
    );
    vec![
        Expr::parse_expression("q/0.1"),
        Expr::parse_expression(&format!("{pe}*q-0.1*{pe}^2*({source})")),
        Expr::parse_expression("J0/1e-4"),
        Expr::parse_expression(&format!("{pe}*J0+1e-4*0.01*{pe}^2*({source})")),
        Expr::parse_expression("J1/1e-4"),
        Expr::parse_expression(&format!("{pe}*J1-1e-4*0.01*{pe}^2*({source})")),
    ]
}

fn problem_with_equations(mesh: Vec<f64>, equations: Vec<Expr>) -> BvpProblem {
    BvpProblem::new(
        equations,
        STATE_NAMES.iter().map(|name| (*name).to_owned()).collect(),
        "z",
        mesh,
        vec![
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
            BvpBoundaryCondition {
                state: 2,
                side: BvpBoundarySide::Left,
                value: 1.0,
            },
            BvpBoundaryCondition {
                state: 3,
                side: BvpBoundarySide::Right,
                value: 0.0,
            },
            BvpBoundaryCondition {
                state: 4,
                side: BvpBoundarySide::Left,
                value: 0.0,
            },
            BvpBoundaryCondition {
                state: 5,
                side: BvpBoundarySide::Right,
                value: 0.0,
            },
        ],
    )
    .expect("standalone MF-LOW system must be a valid RST problem")
}

fn solve_standalone(
    mesh: Vec<f64>,
    guess: DMatrix<f64>,
    settings: BvpSolverSettings,
) -> (Vec<f64>, DMatrix<f64>) {
    solve_standalone_for_case(low_case(), mesh, guess, settings)
}

fn solve_standalone_for_case(
    case: Case,
    mesh: Vec<f64>,
    guess: DMatrix<f64>,
    settings: BvpSolverSettings,
) -> (Vec<f64>, DMatrix<f64>) {
    solve_standalone_problem(mesh, guess, standalone_equations_for(case), settings)
}

fn solve_standalone_problem(
    mesh: Vec<f64>,
    guess: DMatrix<f64>,
    equations: Vec<Expr>,
    settings: BvpSolverSettings,
) -> (Vec<f64>, DMatrix<f64>) {
    let problem = problem_with_equations(mesh, equations);
    let mut solver = BvpSolver::new(problem, settings.with_node_aligned_initial_guess(guess))
        .expect("standalone RST solver should prepare");
    solver
        .solve()
        .expect("standalone RST solver should converge");
    let view = solver.solution_view().expect("standalone solution view");
    let mesh = view.mesh().to_vec();
    let values = DMatrix::from_fn(6, mesh.len(), |state, node| {
        view.value(node, state).expect("all six states are present")
    });
    (mesh, values)
}

fn assert_equation_parity(
    reactor: &super::SimpleReactorTask,
    mesh: &[f64],
    guess: &DMatrix<f64>,
    case: Case,
) {
    let canonical_equations = standalone_equations_for(case);
    assert_eq!(reactor.solver.eq_system.len(), canonical_equations.len());
    let mut max_scaled = [0.0_f64; 6];
    let mut worst = [(0usize, 0.0_f64, 0.0_f64); 6];
    for node in 0..mesh.len() {
        let state: Vec<_> = (0..6).map(|index| guess[(index, node)]).collect();
        for row in 0..6 {
            let production = reactor.solver.eq_system[row].eval_expression(&STATE_NAMES, &state);
            let standalone = canonical_equations[row].eval_expression(&STATE_NAMES, &state);
            let scaled = (production - standalone).abs() / (1.0 + standalone.abs());
            if scaled > max_scaled[row] {
                max_scaled[row] = scaled;
                worst[row] = (node, production, standalone);
            }
        }
    }
    println!(
        "MF equation parity: Pe_q={:.17e} max_scaled={max_scaled:?} worst={worst:?}",
        reactor.Pe_q,
    );
    if max_scaled.iter().copied().fold(0.0, f64::max) > 2e-11 {
        println!(
            "MF-LOW equations: production=[{}] standalone=[{}]",
            reactor
                .solver
                .eq_system
                .iter()
                .map(ToString::to_string)
                .collect::<Vec<_>>()
                .join(" | "),
            canonical_equations
                .iter()
                .map(ToString::to_string)
                .collect::<Vec<_>>()
                .join(" | "),
        );
    }
    assert!(
        max_scaled.iter().copied().fold(0.0, f64::max) <= 2e-11,
        "KiThe-generated equations differ from independent RST equations: {max_scaled:?}"
    );
}

fn assert_solution_parity(
    label: &str,
    standalone_mesh: &[f64],
    standalone_values: &DMatrix<f64>,
    comparison_mesh: &[f64],
    reactor: &super::SimpleReactorTask,
) -> [f64; 6] {
    let physical_mesh = reactor
        .solver
        .x_mesh
        .as_ref()
        .expect("published physical mesh");
    let physical_values = reactor
        .solver
        .get_solution()
        .expect("published physical solution");
    assert_eq!(physical_values.ncols(), 6);
    assert_eq!(physical_values.nrows(), physical_mesh.len());

    let length = reactor.L;
    let facade_mesh = physical_mesh.iter().map(|x| x / length).collect::<Vec<_>>();
    let facade_values = DMatrix::from_fn(6, facade_mesh.len(), |state, node| match state {
        0 => (physical_values[(node, state)] - 300.0) / 15.0,
        1 => physical_values[(node, state)] * length / 15.0,
        3 | 5 => physical_values[(node, state)] * length,
        _ => physical_values[(node, state)],
    });
    let matching_node = |mesh: &[f64], z: f64| {
        let index = mesh.partition_point(|value| *value < z);
        [index.checked_sub(1), (index < mesh.len()).then_some(index)]
            .into_iter()
            .flatten()
            .find(|candidate| {
                (mesh[*candidate] - z).abs() <= 16.0 * f64::EPSILON * z.abs().max(1.0)
            })
    };
    let mut max_state_scaled = [0.0_f64; 6];
    let mut max_mesh_scaled = 0.0_f64;
    for &z in comparison_mesh {
        let direct_node = matching_node(standalone_mesh, z)
            .unwrap_or_else(|| panic!("{label}: standalone route dropped seed node z={z:.17e}"));
        let facade_node = matching_node(&facade_mesh, z)
            .unwrap_or_else(|| panic!("{label}: facade route dropped seed node z={z:.17e}"));
        max_mesh_scaled = max_mesh_scaled
            .max((facade_mesh[facade_node] - standalone_mesh[direct_node]).abs() / (1.0 + z.abs()));
        for state in 0..6 {
            let facade_internal = facade_values[(state, facade_node)];
            let direct = standalone_values[(state, direct_node)];
            max_state_scaled[state] = max_state_scaled[state]
                .max((facade_internal - direct).abs() / (1.0 + direct.abs()));
        }
    }
    println!(
        "MF parity {label}: direct_nodes={} facade_nodes={} seed_mesh_scaled={max_mesh_scaled:.3e} max_state_scaled={max_state_scaled:?}",
        standalone_mesh.len(),
        facade_mesh.len(),
    );
    assert!(
        max_mesh_scaled <= 2e-12,
        "{label}: facade altered the supplied seed mesh"
    );
    for (state, error) in max_state_scaled.iter().copied().enumerate() {
        let tolerance = 2e-10;
        assert!(
            error <= tolerance,
            "{label}: facade altered RST state {} beyond {tolerance:e}: {max_state_scaled:?}",
            STATE_NAMES[state]
        );
    }
    max_state_scaled
}

fn low_case() -> Case {
    cases()[0]
}

fn clustered_input() -> (Vec<f64>, DMatrix<f64>) {
    let rows = fixture_rows(LOW_SEED);
    let mesh = rows.iter().map(|row| row[0]).collect::<Vec<_>>();
    let profile = DMatrix::from_fn(6, rows.len(), |state, node| rows[node][state + 1]);
    (mesh, profile)
}

fn sci_settings(max_nodes: usize) -> BvpSolverSettings {
    sci_settings_with_tolerance(max_nodes, 2e-7)
}

fn sci_settings_with_tolerance(max_nodes: usize, tolerance: f64) -> BvpSolverSettings {
    let options = BvpSciOptions::default()
        .with_execution(BvpSciExecution::Lambdify)
        .with_assembly(BvpSciAssembly::ExprLegacy)
        .with_matrix_layout(BvpSciMatrixLayout::Sparse)
        .with_tolerance(tolerance)
        .with_boundary_tolerance(tolerance)
        .with_limits(max_nodes, 100, 20)
        .with_telemetry(BvpSciTelemetryMode::Detailed);
    BvpSolverSettings::Sci {
        options,
        sparse_factorization: BvpSciSparseFactorization::StandardFaer,
        aot_config: None,
    }
}

fn correction_newton_settings() -> BvpSolverSettings {
    let states = STATE_NAMES;
    let mut options = CorrectionNewtonOptions::default();
    options.matrix_layout = CorrectionNewtonMatrixLayout::BandedAuto;
    options.absolute_tolerance = 2e-7;
    options.default_relative_tolerance = 2e-7;
    options.relative_tolerances = states
        .iter()
        .map(|state| ((*state).to_owned(), 2e-7))
        .collect();
    options.max_iterations = 100;
    options
        .bounds
        .insert("Teta".into(), (positive_teta_lower_bound(), 1e6));
    options.logging_level = SolverLogLevel::Off;
    BvpSolverSettings::CorrectionNewton(options)
}

fn damped_settings(grid_method: GridRefinementMethod, max_refinements: usize) -> BvpSolverSettings {
    damped_settings_with_tolerance(grid_method, max_refinements, 2e-7)
}

fn damped_settings_with_tolerance(
    grid_method: GridRefinementMethod,
    max_refinements: usize,
    tolerance: f64,
) -> BvpSolverSettings {
    let states = STATE_NAMES;
    let options = DampedSolverOptions::banded_damped()
        .with_banded_lambdify()
        .trapezoid_derivative()
        .with_abs_tolerance(tolerance)
        .with_rel_tolerance(
            states
                .iter()
                .map(|state| ((*state).to_owned(), tolerance))
                .collect(),
        )
        .with_bounds(
            states
                .iter()
                .map(|state| {
                    let lower = if *state == "Teta" {
                        positive_teta_lower_bound()
                    } else {
                        -1e6
                    };
                    ((*state).to_owned(), (lower, 1e6))
                })
                .collect(),
        )
        .with_max_iterations(100)
        .with_strategy_params(Some(SolverParams {
            max_jac: Some(3),
            max_damp_iter: Some(10),
            damp_factor: Some(0.5),
            adaptive: Some(AdaptiveGridConfig {
                version: 1,
                max_refinements,
                grid_method,
            }),
        }));
    BvpSolverSettings::Damped(options)
}

fn standalone_reference_errors(
    case: Case,
    mesh: &[f64],
    values: &DMatrix<f64>,
    length: f64,
) -> (f64, f64, f64, f64) {
    let reference = super::reference_profile(case);
    let (mut t_error, mut y_error, mut q_error, mut closure_error) =
        (0.0_f64, 0.0_f64, 0.0_f64, 0.0_f64);
    let mut q_scale = 0.0_f64;
    for node in 0..mesh.len() {
        let xi = mesh[node] * case.xi_max;
        let theta = super::interpolate(&reference, xi, |point| point.theta);
        let p = super::interpolate(&reference, xi, |point| point.p);
        let y_a = super::interpolate(&reference, xi, |point| point.y_a);
        let temperature = 300.0 + 15.0 * values[(0, node)];
        let heat_flux = values[(1, node)] * 15.0 / length;
        let mass_fraction_a = values[(2, node)];
        let mass_fraction_b = values[(4, node)];
        t_error = t_error.max((temperature - (300.0 + 15.0 * theta)).abs());
        y_error = y_error
            .max((mass_fraction_a - y_a).abs())
            .max((mass_fraction_b - (1.0 - y_a)).abs());
        q_error = q_error.max((heat_flux - 0.01 * 1000.0 * 15.0 * p).abs());
        closure_error = closure_error.max((mass_fraction_a + mass_fraction_b - 1.0).abs());
        q_scale = q_scale.max((0.01 * 1000.0 * 15.0 * p).abs());
    }
    (t_error, y_error, q_error / q_scale, closure_error)
}

fn passes_mflow_reference_gate(errors: (f64, f64, f64, f64), values: &DMatrix<f64>) -> bool {
    let min_c0 = values.row(2).iter().copied().fold(f64::INFINITY, f64::min);
    let min_c1 = values.row(4).iter().copied().fold(f64::INFINITY, f64::min);
    errors.0 < 0.5
        && errors.1 < 2e-3
        && errors.2 < 2e-3
        && errors.3 < 2e-8
        && min_c0 >= -1e-12
        && min_c1 >= -1e-12
}

#[test]
fn mf1963_mflow_reference_fixtures_report_cross_workflow_drift() {
    let reference = super::reference_profile(low_case());
    let mut max_theta = (0.0_f64, 0.0_f64);
    let mut max_slope = (0.0_f64, 0.0_f64);
    let mut max_species = (0.0_f64, 0.0_f64);

    for (line_index, line) in RST_LOW_SCIPY_REFERENCE.lines().enumerate() {
        if line_index == 0 || line.trim().is_empty() {
            continue;
        }
        let row = line
            .split(',')
            .map(|value| value.parse::<f64>().expect("numeric RST reference cell"))
            .collect::<Vec<_>>();
        assert_eq!(row.len(), 7);
        let xi = row[0] * LOW_PE;
        let theta_reference = interpolate(&reference, xi, |point| point.theta);
        let slope_reference = interpolate(&reference, xi, |point| point.p);
        let species_reference = interpolate(&reference, xi, |point| point.y_a);
        let theta_error = (row[1] - theta_reference).abs();
        let slope_error = (row[2] / (0.1 * LOW_PE) - slope_reference).abs();
        let species_error = (row[3] - species_reference).abs();
        if theta_error > max_theta.0 {
            max_theta = (theta_error, xi);
        }
        if slope_error > max_slope.0 {
            max_slope = (slope_error, xi);
        }
        if species_error > max_species.0 {
            max_species = (species_error, xi);
        }
    }

    println!(
        "MF-LOW reference fixture drift: max_abs_theta={:.12e} at xi={:.9} max_abs_p={:.12e} at xi={:.9} max_abs_YA={:.12e} at xi={:.9}",
        max_theta.0, max_theta.1, max_slope.0, max_slope.1, max_species.0, max_species.1,
    );
    assert!((1.17e-3..1.19e-3).contains(&max_theta.0));
    assert!((0.99e-3..1.02e-3).contains(&max_slope.0));
    assert!((1.17e-5..1.19e-5).contains(&max_species.0));
}

fn frozen_settings(
    strategy_params: Option<HashMap<String, Option<Vec<f64>>>>,
) -> BvpSolverSettings {
    let mut options = FrozenSolverOptions::banded_frozen()
        .with_banded_lambdify()
        .trapezoid_derivative()
        .with_tolerance(2e-7)
        .with_max_iterations(100)
        .with_bvp_telemetry_mode(BvpTelemetryMode::Detailed)
        .with_bvp_logging_config(
            BvpLoggingConfig::new(BvpLoggingMode::Detailed).with_max_events(128),
        );
    if let Some(params) = strategy_params {
        options = options.with_strategy_params(Some(params));
    }
    BvpSolverSettings::Frozen(options)
}

fn solve_production(
    case: Case,
    mesh: Vec<f64>,
    guess: DMatrix<f64>,
    settings: BvpSolverSettings,
) -> super::SimpleReactorTask {
    try_solve_production(case, mesh, guess, settings)
        .unwrap_or_else(|error| panic!("production BVP solve failed: {error}"))
}

fn try_solve_production(
    case: Case,
    mesh: Vec<f64>,
    guess: DMatrix<f64>,
    settings: BvpSolverSettings,
) -> Result<SimpleReactorTask, ReactorError> {
    let mut reactor = reactor_for(case);
    reactor.setup_bvp()?;
    reactor.solver = reactor
        .solver
        .with_bvp_telemetry_mode(BvpApiTelemetryMode::Detailed);
    reactor
        .solver
        .solve_with_bvp_api_on_mesh(mesh, guess, settings)?;
    reactor.postprocessing()?;
    assert!(
        reactor.solver.last_telemetry.is_some(),
        "facade telemetry override must publish a route snapshot"
    );
    Ok(reactor)
}

#[test]
#[ignore = "paired MF-LOW direct-RST versus production-facade parity matrix; includes expensive adaptive solves"]
fn mf1963_mflow_production_facade_matches_standalone_rst_solver_outputs() {
    let routes: Vec<(&str, Box<dyn Fn() -> BvpSolverSettings>)> = vec![
        ("BVP_sci/Sparse", Box::new(|| sci_settings(6000))),
        (
            "CorrectionNewton/BandedAuto",
            Box::new(correction_newton_settings),
        ),
        (
            "Damped/TwoPoint(0.01,0.01,2.0)",
            Box::new(|| damped_settings(GridRefinementMethod::TwoPoint(0.01, 0.01, 2.0), 2)),
        ),
        (
            "Damped/Pearson(1e-3,1.4)",
            Box::new(|| damped_settings(GridRefinementMethod::Pearson(1e-3, 1.4), 2)),
        ),
        (
            "Damped/GrcarSmooke(1e-3,1e-2,1.4)",
            Box::new(|| damped_settings(GridRefinementMethod::GrcarSmooke(1e-3, 1e-2, 1.4), 2)),
        ),
        (
            "Damped/DoublePoints/4",
            Box::new(|| damped_settings(GridRefinementMethod::DoublePoints, 4)),
        ),
        (
            "Frozen/at_high_morm=1e-6",
            Box::new(|| {
                let BvpSolverSettings::Frozen(options) = frozen_settings(Some(HashMap::from([(
                    "at_high_morm".to_owned(),
                    Some(vec![1e-6]),
                )]))) else {
                    unreachable!()
                };
                BvpSolverSettings::Frozen(options)
            }),
        ),
    ];

    let (mesh, guess) = clustered_input();
    let mut equation_probe = reactor_for(low_case());
    equation_probe
        .setup_bvp()
        .expect("production equation setup");
    assert_equation_parity(&equation_probe, &mesh, &guess, low_case());
    let production_equations = equation_probe.solver.eq_system.clone();

    for (label, settings) in routes {
        let (rst_mesh, rst_values) = solve_standalone(mesh.clone(), guess.clone(), settings());
        let (generated_mesh, generated_values) = solve_standalone_problem(
            mesh.clone(),
            guess.clone(),
            production_equations.clone(),
            settings(),
        );
        let reactor = solve_production(low_case(), mesh.clone(), guess.clone(), settings());
        assert_solution_parity(label, &generated_mesh, &generated_values, &mesh, &reactor);

        let errors = production_profile_errors(low_case(), &reactor);
        let direct_errors =
            standalone_reference_errors(low_case(), &rst_mesh, &rst_values, reactor.L);
        let generated_errors =
            standalone_reference_errors(low_case(), &generated_mesh, &generated_values, reactor.L);
        println!(
            "MF-LOW parity {label}: canonical RST dT={:.12e} K dY={:.12e} relative_dq={:.12e}; generated-equation RST dT={:.12e} K dY={:.12e} relative_dq={:.12e}; facade dT={:.12e} K dY={:.12e} relative_dq={:.12e}",
            direct_errors.0,
            direct_errors.1,
            direct_errors.2,
            generated_errors.0,
            generated_errors.1,
            generated_errors.2,
            errors.0,
            errors.1,
            errors.2,
        );
        let standalone_pass = passes_mflow_reference_gate(direct_errors, &rst_values);
        let generated_pass = passes_mflow_reference_gate(generated_errors, &generated_values);
        let facade_values = reactor.solver.get_solution().unwrap();
        let facade_pass = passes_mflow_reference_gate(errors, &facade_values.transpose());
        assert_eq!(
            facade_pass, standalone_pass,
            "{label}: production facade changed the canonical standalone physical acceptance classification"
        );
        assert_eq!(
            facade_pass, generated_pass,
            "{label}: facade changed the same-equation direct-RST physical acceptance classification"
        );
        for (metric, facade, generated, tolerance) in [
            ("temperature", errors.0, generated_errors.0, 1e-7),
            ("species", errors.1, generated_errors.1, 1e-10),
            ("heat flux", errors.2, generated_errors.2, 1e-10),
        ] {
            assert!(
                facade <= generated + tolerance,
                "{label}: facade worsened {metric} vs same-equation RST: facade={facade:.17e}, generated={generated:.17e}, allowance={tolerance:.3e}"
            );
        }
    }
}

#[test]
#[ignore = "MF-LOW 2049-node clustered production-facade solver comparison"]
fn mf1963_mflow_clustered_mesh_runs_through_production_facade_for_selected_solvers() {
    for (name, settings) in [
        ("BVP_sci/Sparse", sci_settings(30_000)),
        ("CorrectionNewton/BandedAuto", correction_newton_settings()),
        (
            "Damped/TwoPoint(0.01,0.01,2.0)",
            damped_settings(GridRefinementMethod::TwoPoint(0.01, 0.01, 2.0), 2),
        ),
        (
            "Damped/Pearson(1e-3,1.4)",
            damped_settings(GridRefinementMethod::Pearson(1e-3, 1.4), 2),
        ),
        (
            "Damped/GrcarSmooke(1e-3,1e-2,1.4)",
            damped_settings(GridRefinementMethod::GrcarSmooke(1e-3, 1e-2, 1.4), 2),
        ),
        (
            "Damped/DoublePoints/4",
            damped_settings(GridRefinementMethod::DoublePoints, 4),
        ),
    ] {
        let (mesh, guess) = clustered_input();
        let reactor = solve_production(low_case(), mesh, guess, settings);
        let solution = reactor.solver.get_solution().expect("published solution");
        let mesh = reactor.solver.x_mesh.as_ref().expect("physical mesh");
        println!(
            "MF-LOW production {name}: min(Y_A)={:.9e} min(Y_B)={:.9e}",
            solution
                .column(2)
                .iter()
                .copied()
                .fold(f64::INFINITY, f64::min),
            solution
                .column(4)
                .iter()
                .copied()
                .fold(f64::INFINITY, f64::min),
        );
        let errors = production_profile_errors(low_case(), &reactor);
        let defects = if name.starts_with("BVP_sci") {
            None
        } else {
            Some(mf_discrete_trapezoid_defects(low_case(), &reactor))
        };
        println!(
            "MF-LOW production {name}: dT={:.3e} K dY={:.3e} relative_dq={:.3e} closure={:.3e}",
            errors.0, errors.1, errors.2, errors.3,
        );
        assert!(solution.iter().all(|value| value.is_finite()));
        assert!(
            solution
                .column(0)
                .iter()
                .all(|teta| 300.0 + 15.0 * teta > 0.0)
        );
        let min_y_a = solution
            .column(2)
            .iter()
            .copied()
            .fold(f64::INFINITY, f64::min);
        let min_y_b = solution
            .column(4)
            .iter()
            .copied()
            .fold(f64::INFINITY, f64::min);
        if name.starts_with("CorrectionNewton") {
            // This is a converged fixed-mesh algebraic root, but not the
            // reference physical branch for this stiff LOW case.
            assert!(
                (-5e-8..-1e-9).contains(&min_y_a),
                "unexpected Y_A minimum {min_y_a:e}"
            );
            assert!(min_y_b >= -1e-12, "unexpected Y_B minimum {min_y_b:e}");
            assert!(
                (100.0..150.0).contains(&errors.0),
                "unexpected dT: {}",
                errors.0
            );
            assert!(
                (0.05..0.15).contains(&errors.1),
                "unexpected dY: {}",
                errors.1
            );
            assert!(
                (0.1..0.2).contains(&errors.2),
                "unexpected relative dq: {}",
                errors.2
            );
            assert!(errors.3 < 2e-8);
            let defects = defects.expect("fixed-mesh defect report");
            println!("MF-LOW production {name}: independent defects={defects:?}");
            assert!(defects.0 < 1e-5 && defects.1 < 1e-5 && defects.2 < 1e-5);
            println!(
                "MF-LOW production {name}: solver converged; physical-reference gate intentionally NOT met"
            );
        } else if name.starts_with("Damped/") {
            let defects = defects.expect("adaptive Damped defect report");
            println!("MF-LOW production {name}: independent defects={defects:?}");
            assert!(defects.0 < 1e-5 && defects.1 < 1e-5 && defects.2 < 1e-5);
            let physical_pass = min_y_a >= -1e-12
                && min_y_b >= -1e-12
                && errors.0 < 0.5
                && errors.1 < 2e-3
                && errors.2 < 2e-3
                && errors.3 < 2e-8;
            if name.contains("DoublePoints") {
                assert!(
                    physical_pass,
                    "MF-LOW DoublePoints should pass physical acceptance"
                );
                println!("MF-LOW production {name}: physical-reference gate PASS");
            } else {
                assert!(
                    !physical_pass,
                    "unexpected MF-LOW physical acceptance for {name}"
                );
                println!(
                    "MF-LOW production {name}: converged discrete root; physical-reference gate NOT met"
                );
            }
        } else {
            assert!(min_y_a >= -1e-12, "nonphysical Y_A minimum {min_y_a:e}");
            assert!(min_y_b >= -1e-12, "nonphysical Y_B minimum {min_y_b:e}");
            assert!(errors.0 < 0.5 && errors.1 < 2e-3 && errors.2 < 2e-3);
            assert!(errors.3 < 2e-8);
        }
        assert_eq!(solution.nrows(), mesh.len());
        assert_eq!(solution.ncols(), 6);
        println!(
            "MF-LOW production {name}: physical_nodes={} x=[{:.9e},{:.9e}] m T=[{:.6e},{:.6e}] K",
            solution.nrows(),
            mesh[0],
            mesh[mesh.len() - 1],
            solution
                .column(0)
                .iter()
                .copied()
                .fold(f64::INFINITY, f64::min),
            solution
                .column(0)
                .iter()
                .copied()
                .fold(f64::NEG_INFINITY, f64::max),
        );
    }
}

#[test]
#[ignore = "MF-HIGH production-facade Damped grid-strategy matrix"]
fn mf1963_mfhigh_damped_grid_strategy_matrix_matches_reference_gate() {
    let case = cases()[1];
    let nodes = 301;
    let mesh = (0..nodes)
        .map(|node| node as f64 / (nodes - 1) as f64)
        .collect::<Vec<_>>();
    let guess = super::reference_seed(case, nodes);
    for (label, method, max_refinements) in [
        (
            "TwoPoint(0.01,0.01,2.0)",
            GridRefinementMethod::TwoPoint(0.01, 0.01, 2.0),
            2,
        ),
        (
            "Pearson(1e-3,1.4)",
            GridRefinementMethod::Pearson(1e-3, 1.4),
            2,
        ),
        (
            "GrcarSmooke(1e-3,1e-2,1.4)",
            GridRefinementMethod::GrcarSmooke(1e-3, 1e-2, 1.4),
            2,
        ),
        ("DoublePoints/4", GridRefinementMethod::DoublePoints, 4),
    ] {
        let reactor = solve_production(
            case,
            mesh.clone(),
            guess.clone(),
            damped_settings(method, max_refinements),
        );
        let errors = production_profile_errors(case, &reactor);
        println!(
            "MF-HIGH production Damped/{label}: nodes={} dT={:.3e} K dY={:.3e} relative_dq={:.3e} closure={:.3e}",
            reactor.solver.get_solution().unwrap().nrows(),
            errors.0,
            errors.1,
            errors.2,
            errors.3,
        );
        assert!(errors.0 < 0.5 && errors.1 < 2e-3 && errors.2 < 2e-3);
        assert!(errors.3 < 2e-8);
    }
}

#[test]
#[ignore = "MF-HIGH production-facade Frozen refresh-policy matrix"]
fn mf1963_mfhigh_frozen_refresh_policy_matrix_matches_reference_gate() {
    let case = cases()[1];
    let nodes = 301;
    let mesh = (0..nodes)
        .map(|node| node as f64 / (nodes - 1) as f64)
        .collect::<Vec<_>>();
    let guess = super::reference_seed(case, nodes);
    let strategies = [
        ("default", None),
        (
            "at_high_morm=1e-6",
            Some(HashMap::from([(
                "at_high_morm".to_owned(),
                Some(vec![1e-6]),
            )])),
        ),
        (
            "every_m=1",
            Some(HashMap::from([("every_m".to_owned(), Some(vec![1.0]))])),
        ),
    ];

    let mut profiles = Vec::new();
    for (label, params) in strategies {
        let reactor = solve_production(case, mesh.clone(), guess.clone(), frozen_settings(params));
        let errors = production_profile_errors(case, &reactor);
        let profile = reactor
            .solver
            .get_solution()
            .expect("published Frozen solution");
        println!(
            "MF-HIGH production Frozen/{label}: nodes={} dT={:.3e} K dY={:.3e} relative_dq={:.3e} closure={:.3e}",
            profile.nrows(),
            errors.0,
            errors.1,
            errors.2,
            errors.3,
        );
        assert!(errors.0 < 0.5 && errors.1 < 2e-3 && errors.2 < 2e-3);
        assert!(errors.3 < 2e-8);
        profiles.push((label, profile.clone()));
    }

    let reference = &profiles[0].1;
    for (label, profile) in profiles.iter().skip(1) {
        let distance = (profile - reference).amax();
        println!("MF-HIGH production Frozen/{label}: max profile delta={distance:.3e}");
        assert!(
            distance < 2e-5,
            "Frozen refresh policy changed the accepted root"
        );
    }
}

#[test]
#[ignore = "MF-LOW production-facade Frozen refresh-on-correction configuration"]
fn mf1963_mflow_frozen_refresh_on_correction_runs_through_production_facade() {
    let (mesh, guess) = clustered_input();
    let reactor = solve_production(
        low_case(),
        mesh,
        guess,
        frozen_settings(Some(HashMap::from([(
            "at_high_morm".to_owned(),
            Some(vec![1e-6]),
        )]))),
    );
    let solution = reactor
        .solver
        .get_solution()
        .expect("published Frozen solution");
    let errors = production_profile_errors(low_case(), &reactor);
    let defects = mf_discrete_trapezoid_defects(low_case(), &reactor);
    let min_y_a = solution
        .column(2)
        .iter()
        .copied()
        .fold(f64::INFINITY, f64::min);
    let min_y_b = solution
        .column(4)
        .iter()
        .copied()
        .fold(f64::INFINITY, f64::min);
    println!(
        "MF-LOW production Frozen/at_high_morm=1e-6: nodes={} min(Y_A)={min_y_a:.3e} min(Y_B)={min_y_b:.3e} dT={:.3e} K dY={:.3e} relative_dq={:.3e} closure={:.3e} defects={defects:?}",
        solution.nrows(),
        errors.0,
        errors.1,
        errors.2,
        errors.3,
    );
    assert!(solution.iter().all(|value| value.is_finite()));
    assert!(
        solution
            .column(0)
            .iter()
            .all(|teta| 300.0 + 15.0 * teta > 0.0)
    );
    assert!(errors.3 < 2e-8);
    assert!(defects.0 < 1e-5 && defects.1 < 1e-5 && defects.2 < 1e-5);
    let physical_pass = min_y_a >= -1e-12
        && min_y_b >= -1e-12
        && errors.0 < 0.5
        && errors.1 < 2e-3
        && errors.2 < 2e-3;
    assert!(
        !physical_pass,
        "LOW Frozen must not be classified as reference-accurate"
    );
    assert!((-5e-8..-1e-9).contains(&min_y_a));
    assert!((100.0..150.0).contains(&errors.0));
    assert!((0.05..0.15).contains(&errors.1));
    assert!((0.1..0.2).contains(&errors.2));
    println!(
        "MF-LOW production Frozen/at_high_morm=1e-6: physical-reference gate {}",
        if physical_pass {
            "PASS"
        } else {
            "NOT MET (alternate discrete branch)"
        },
    );
}

#[test]
fn mf1963_mfhigh_uniform_mesh_runs_through_production_facade() {
    let case = cases()[1];
    let nodes = 301;
    let mesh = (0..nodes)
        .map(|node| node as f64 / (nodes - 1) as f64)
        .collect::<Vec<_>>();
    let guess = super::reference_seed(case, nodes);
    let reactor = solve_production(case, mesh, guess, correction_newton_settings());
    let errors = production_profile_errors(case, &reactor);
    println!(
        "MF-HIGH production CorrectionNewton/BandedAuto: dT={:.3e} K dY={:.3e} relative_dq={:.3e} closure={:.3e}",
        errors.0, errors.1, errors.2, errors.3,
    );
    assert!(errors.0 < 0.5 && errors.1 < 2e-3 && errors.2 < 2e-3);
    assert!(errors.3 < 2e-8);
}

#[test]
#[ignore = "MF-LOW 2049-node uniform control reaches the collocation node cap"]
fn mf1963_mflow_uniform_mesh_control_reaches_node_cap_through_production_facade() {
    let case = low_case();
    let nodes = 2049;
    let mesh = (0..nodes)
        .map(|node| node as f64 / (nodes - 1) as f64)
        .collect::<Vec<_>>();
    let guess = super::reference_seed(case, nodes);
    let mut reactor = reactor_for(case);
    reactor.setup_bvp().expect("MF-LOW physical BVP setup");
    let error = reactor
        .solver
        .solve_with_bvp_api_on_mesh(mesh, guess, sci_settings(6000))
        .expect_err("uniform MF-LOW control should reach the configured node cap");
    let diagnostic = format!("{error:?}");
    println!("MF-LOW production uniform-mesh control: {diagnostic}");
    assert!(
        diagnostic.to_ascii_lowercase().contains("node")
            || diagnostic.to_ascii_lowercase().contains("mesh")
            || diagnostic.to_ascii_lowercase().contains("limit"),
        "expected node-cap failure, got {diagnostic}"
    );
    assert!(reactor.solver.get_solution().is_none());
    assert!(reactor.solver.x_mesh.is_none());
}

#[test]
fn mf1963_explicit_mesh_rejection_preserves_published_solution() {
    let case = cases()[1];
    let nodes = 301;
    let mesh = (0..nodes)
        .map(|node| node as f64 / (nodes - 1) as f64)
        .collect::<Vec<_>>();
    let guess = super::reference_seed(case, nodes);
    let mut reactor = reactor_for(case);
    reactor.setup_bvp().expect("MF-HIGH physical BVP setup");
    reactor
        .solver
        .solve_with_bvp_api_on_mesh(mesh.clone(), guess.clone(), correction_newton_settings())
        .expect("initial production solve");
    reactor
        .postprocessing()
        .expect("initial dimensional postprocessing");
    let published_solution = reactor.solver.get_solution().unwrap().clone();
    let published_mesh = reactor.solver.x_mesh.as_ref().unwrap().clone();

    let mut invalid_guess = guess;
    invalid_guess[(0, 0)] = 1e-3;
    let error = reactor
        .solver
        .solve_with_bvp_api_on_mesh(mesh, invalid_guess, correction_newton_settings())
        .expect_err("guess violating theta(0)=0 must be rejected");
    let message = error.to_string();
    assert!(
        message.contains("violates `Teta` boundary condition"),
        "{error}"
    );
    assert_eq!(reactor.solver.get_solution().unwrap(), &published_solution);
    assert_eq!(reactor.solver.x_mesh.as_ref().unwrap(), &published_mesh);
}

#[test]
fn mf1963_explicit_mesh_rejects_boundary_inconsistent_profile_before_solver() {
    let case = low_case();
    let (mesh, mut guess) = clustered_input();
    guess[(0, 0)] = 1e-3;
    let mut reactor = reactor_for(case);
    reactor.setup_bvp().expect("MF-LOW physical BVP setup");
    let error = reactor
        .solver
        .solve_with_bvp_api_on_mesh(mesh, guess, correction_newton_settings())
        .expect_err("guess violating theta(0)=0 must be rejected");
    let message = error.to_string();
    assert!(
        message.contains("violates `Teta` boundary condition"),
        "{error}"
    );
    assert!(reactor.solver.get_solution().is_none());
    assert!(reactor.solver.x_mesh.is_none());
}

#[test]
fn mf1963_explicit_mesh_rejects_nonmonotone_mesh_before_solver() {
    let case = low_case();
    let (mut mesh, guess) = clustered_input();
    mesh.swap(10, 11);
    let mut reactor = reactor_for(case);
    reactor.setup_bvp().expect("MF-LOW physical BVP setup");
    let error = reactor
        .solver
        .solve_with_bvp_api_on_mesh(mesh, guess, correction_newton_settings())
        .expect_err("nonmonotone mesh must be rejected");
    assert!(
        error.to_string().contains("finite and strictly increasing"),
        "{error}"
    );
    assert!(reactor.solver.get_solution().is_none());
}

#[test]
fn mf1963_explicit_mesh_validates_domain_profile_shape_and_finiteness() {
    let case = cases()[1];
    let nodes = 31;
    let mesh = (0..nodes)
        .map(|node| node as f64 / (nodes - 1) as f64)
        .collect::<Vec<_>>();
    let guess = super::reference_seed(case, nodes);

    let mut wrong_domain = mesh.clone();
    wrong_domain[0] = 1e-3;
    let mut reactor = reactor_for(case);
    reactor.setup_bvp().expect("MF-HIGH physical BVP setup");
    let error = reactor
        .solver
        .solve_with_bvp_api_on_mesh(wrong_domain, guess.clone(), correction_newton_settings())
        .expect_err("mesh must span the configured solver domain");
    assert!(
        error
            .to_string()
            .contains("must span the configured domain")
    );

    let mut reactor = reactor_for(case);
    reactor.setup_bvp().expect("MF-HIGH physical BVP setup");
    let error = reactor
        .solver
        .solve_with_bvp_api_on_mesh(
            mesh.clone(),
            DMatrix::zeros(6, nodes - 1),
            correction_newton_settings(),
        )
        .expect_err("profile must contain one value per state and mesh node");
    assert!(error.to_string().contains("must have shape 6x31"));

    let mut nonfinite_guess = guess;
    nonfinite_guess[(1, 10)] = f64::NAN;
    let mut reactor = reactor_for(case);
    reactor.setup_bvp().expect("MF-HIGH physical BVP setup");
    let error = reactor
        .solver
        .solve_with_bvp_api_on_mesh(mesh, nonfinite_guess, correction_newton_settings())
        .expect_err("non-finite profile must be rejected");
    assert!(error.to_string().contains("non-finite value"));
}
