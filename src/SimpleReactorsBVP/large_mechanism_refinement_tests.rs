//! Capability probe for a connected acetaldehyde-oxidation subset of the
//! bundled GRI30 mechanism. Thermochemistry is resolved through public
//! `SubsData`; only one-way elementary reactions are included because the
//! current `FastElemReact` BVP model cannot represent reversible, falloff, or
//! third-body rows faithfully.

mod large_mechanism_fixture {
    use crate as root;
    include!("../../benches/support/large_mechanism_fixture.rs");
}

use std::collections::BTreeMap;

use super::SimpleReactorBVP::{NrbvpHandoffConfig, SimpleReactorTask};
use crate::Kinetics::molmass::calculate_molar_mass;
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
use large_mechanism_fixture::{
    REACTION_IDS, THERMO_TEMPERATURE_K, build_reactor, inlet_initial_guess, load_reactions,
    parse_equation, resolved_thermochemistry,
};

fn solve_damped(max_refinements: usize, nodes: usize) -> SimpleReactorTask {
    solve_damped_with_grid(max_refinements, nodes, GridRefinementMethod::DoublePoints)
}

fn solve_damped_with_grid(
    max_refinements: usize,
    nodes: usize,
    grid_method: GridRefinementMethod,
) -> SimpleReactorTask {
    let mut reactor = build_reactor();
    reactor.setup_bvp().expect("large-mechanism BVP setup");
    reactor.solver.arg_name = "z".to_string();
    reactor.solver.x_range = (0.0, 1.0);
    reactor.solver = reactor
        .solver
        .with_bvp_telemetry_mode(BvpApiTelemetryMode::Counters);
    let initial = inlet_initial_guess(&reactor, nodes);
    let strategy = SolverParams {
        max_jac: Some(1),
        max_damp_iter: Some(12),
        adaptive: (max_refinements > 0).then_some(AdaptiveGridConfig {
            version: 1,
            max_refinements,
            grid_method,
        }),
        ..SolverParams::default()
    };
    let config = NrbvpHandoffConfig::new(
        initial,
        0.0,
        1.0,
        nodes,
        "trapezoid".to_string(),
        "Damped".to_string(),
        Some(strategy),
        None,
        "Sparse".to_string(),
        1e-5,
        Some(
            reactor
                .solver
                .unknowns
                .iter()
                .map(|name| (name.clone(), 1e-5))
                .collect(),
        ),
        150,
        Some(
            reactor
                .solver
                .unknowns
                .iter()
                .map(|name| {
                    let bounds = if name == "Teta" {
                        (-1.0, 3.0)
                    } else if name == "q" || name.starts_with('J') {
                        (-1.0e4, 1.0e4)
                    } else {
                        (-0.1, 1.1)
                    };
                    (name.clone(), bounds)
                })
                .collect(),
        ),
        Some("off".to_string()),
    )
    .with_solver_backend_config(
        ReactorBvpSolverConfig::sparse_lambdify()
            .with_symbolic_backend(ReactorBvpSymbolicBackend::ExprLegacy),
    );
    reactor
        .solver
        .solve_with_bvp_api(config)
        .expect("large-mechanism Damped BVP solve");
    let telemetry = reactor
        .solver
        .last_telemetry
        .as_ref()
        .expect("Damped solve telemetry");
    let RustedSciThe::numerical::BVP_api::BvpSolverTelemetrySnapshot::Damped(_telemetry) =
        telemetry
    else {
        panic!("expected Damped telemetry snapshot");
    };
    reactor
        .postprocessing()
        .expect("publish dimensional large-mechanism Damped profile");
    reactor
}

fn solve_correction_newton(nodes: usize) -> Result<SimpleReactorTask, String> {
    let mut reactor = build_reactor();
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
            let bounds = if name == "Teta" {
                (-1.0, 3.0)
            } else if name == "q" || name.starts_with('J') {
                (-1.0e4, 1.0e4)
            } else {
                (0.0, 1.1)
            };
            (name.clone(), bounds)
        })
        .collect();
    let result = reactor.solver.solve_with_hybrid_settings(
        BvpSolverSettings::CorrectionNewton(options),
        inlet_initial_guess(&reactor, nodes),
        0.0,
        1.0,
        nodes,
    );
    if let Err(error) = result {
        return Err(format!(
            "CorrectionNewton solve: {error}; partial telemetry={:?}",
            reactor.solver.last_telemetry
        ));
    }
    reactor
        .postprocessing()
        .map_err(|error| format!("CorrectionNewton postprocessing: {error}"))?;
    Ok(reactor)
}

#[derive(Debug)]
struct PublishedMetrics {
    nodes: usize,
    min_temperature_k: f64,
    max_temperature_k: f64,
    min_mass_fraction: f64,
    max_mass_fraction: f64,
    max_mass_fraction_sum_error: f64,
    max_species_change: f64,
}

fn published_metrics(reactor: &SimpleReactorTask) -> PublishedMetrics {
    let solution = reactor
        .solver
        .get_solution()
        .expect("published solution matrix");
    let mesh = reactor.solver.x_mesh.as_ref().expect("published mesh");
    assert_eq!(solution.nrows(), mesh.len());
    assert_eq!(solution.ncols(), reactor.solver.unknowns.len());
    assert!(solution.iter().all(|value| value.is_finite()));
    let temperature_column = reactor
        .solver
        .unknowns
        .iter()
        .position(|name| name == "T")
        .or_else(|| {
            reactor
                .solver
                .unknowns
                .iter()
                .position(|name| name == "Teta")
        })
        .expect("temperature column");
    let temperatures = solution.column(temperature_column);
    let mut metrics = PublishedMetrics {
        nodes: solution.nrows(),
        min_temperature_k: temperatures.iter().copied().fold(f64::INFINITY, f64::min),
        max_temperature_k: temperatures
            .iter()
            .copied()
            .fold(f64::NEG_INFINITY, f64::max),
        min_mass_fraction: f64::INFINITY,
        max_mass_fraction: f64::NEG_INFINITY,
        max_mass_fraction_sum_error: 0.0,
        max_species_change: 0.0,
    };
    for node in 0..solution.nrows() {
        let mut fraction_sum = 0.0;
        for (species_index, name) in reactor.kindata.substances.iter().enumerate() {
            let fraction = solution[(node, 2 + 2 * species_index)];
            metrics.min_mass_fraction = metrics.min_mass_fraction.min(fraction);
            metrics.max_mass_fraction = metrics.max_mass_fraction.max(fraction);
            metrics.max_species_change = metrics
                .max_species_change
                .max((fraction - reactor.boundary_condition[name]).abs());
            fraction_sum += fraction;
        }
        metrics.max_mass_fraction_sum_error = metrics
            .max_mass_fraction_sum_error
            .max((fraction_sum - 1.0).abs());
    }
    metrics
}

fn assert_physical_profile(metrics: &PublishedMetrics, route: &str) {
    assert!(
        (200.0..=4000.0).contains(&metrics.min_temperature_k)
            && (200.0..=4000.0).contains(&metrics.max_temperature_k),
        "{route} temperature range [{}, {}] K is not physically admissible",
        metrics.min_temperature_k,
        metrics.max_temperature_k
    );
    assert!(
        metrics.min_mass_fraction >= -1e-8 && metrics.max_mass_fraction <= 1.0 + 1e-5,
        "{route} mass fractions [{}, {}] escape the physical envelope",
        metrics.min_mass_fraction,
        metrics.max_mass_fraction
    );
    assert!(
        metrics.max_mass_fraction_sum_error < 1e-3,
        "{route} mass-fraction closure error={}",
        metrics.max_mass_fraction_sum_error
    );
    assert!(
        metrics.max_species_change > 1e-8,
        "{route} did not exhibit reaction progress"
    );
}

fn solve_bvp_sci_reference() -> SimpleReactorTask {
    let mut reactor = build_reactor();
    reactor.setup_bvp().expect("large-mechanism BVP setup");
    reactor.solver.arg_name = "z".to_string();
    reactor.solver.x_range = (0.0, 1.0);
    reactor
        .solver
        .solve_BVPsci(inlet_initial_guess(&reactor, 33), 33, 256, 1e-5)
        .expect("BVP_sci reference solve");
    reactor
        .postprocessing()
        .expect("publish dimensional BVP_sci reference");
    reactor
}

fn interpolate_column(mesh: &[f64], values: &[f64], x: f64) -> f64 {
    let upper = mesh.partition_point(|point| *point < x).min(mesh.len() - 1);
    if upper == 0 {
        return values[0];
    }
    let lower = upper - 1;
    let span = mesh[upper] - mesh[lower];
    if span <= 0.0 {
        return values[upper];
    }
    let weight = (x - mesh[lower]) / span;
    values[lower] + weight * (values[upper] - values[lower])
}

fn profile_errors(candidate: &SimpleReactorTask, reference: &SimpleReactorTask) -> [f64; 4] {
    let candidate_solution = candidate.solver.get_solution().expect("candidate solution");
    let candidate_mesh = candidate
        .solver
        .x_mesh
        .as_ref()
        .expect("candidate mesh")
        .as_slice();
    let reference_solution = reference.solver.get_solution().expect("reference solution");
    let reference_mesh = reference
        .solver
        .x_mesh
        .as_ref()
        .expect("reference mesh")
        .as_slice();
    let reference_columns: Vec<Vec<f64>> = (0..reference_solution.ncols())
        .map(|column| reference_solution.column(column).iter().copied().collect())
        .collect();
    let mut errors = [0.0_f64; 4];
    for (node, x) in candidate_mesh.iter().copied().enumerate() {
        let expected_temperature = interpolate_column(reference_mesh, &reference_columns[0], x);
        let expected_heat_flux = interpolate_column(reference_mesh, &reference_columns[1], x);
        errors[0] = errors[0].max((candidate_solution[(node, 0)] - expected_temperature).abs());
        errors[3] = errors[3].max((candidate_solution[(node, 1)] - expected_heat_flux).abs());
        for species_index in 0..candidate.kindata.substances.len() {
            let mass_fraction_column = 2 + 2 * species_index;
            let flux_column = mass_fraction_column + 1;
            let expected_fraction =
                interpolate_column(reference_mesh, &reference_columns[mass_fraction_column], x);
            let expected_flux =
                interpolate_column(reference_mesh, &reference_columns[flux_column], x);
            errors[1] = errors[1]
                .max((candidate_solution[(node, mass_fraction_column)] - expected_fraction).abs());
            errors[2] =
                errors[2].max((candidate_solution[(node, flux_column)] - expected_flux).abs());
        }
    }
    errors
}

#[test]
fn larger_cantera_subset_has_resolved_thermochemistry_and_balanced_reactions() {
    let (reactions, species) = load_reactions();
    let (thermo, molar_masses, provenance) = resolved_thermochemistry(&species);
    for reaction in &reactions {
        let (reactants, products) = parse_equation(&reaction.equation);
        let mut elemental_balance = BTreeMap::<String, f64>::new();
        let mut delta_h = 0.0;
        for (name, coefficient) in reactants {
            let (_, elements) = calculate_molar_mass(name.clone(), None).unwrap();
            let enthalpy = thermo[&name].1;
            delta_h -= coefficient * enthalpy;
            for (element, count) in elements {
                *elemental_balance.entry(element).or_default() -= coefficient * count as f64;
            }
        }
        for (name, coefficient) in products {
            let (_, elements) = calculate_molar_mass(name.clone(), None).unwrap();
            delta_h += coefficient * thermo[&name].1;
            for (element, count) in elements {
                *elemental_balance.entry(element).or_default() += coefficient * count as f64;
            }
        }
        assert!(
            elemental_balance.values().all(|value| value.abs() < 1e-12),
            "reaction {} is not element-balanced: {elemental_balance:?}",
            reaction.equation
        );
        assert!(reaction.a_si.is_finite() && reaction.a_si > 0.0);
        assert!(reaction.activation_energy_j_per_mol.is_finite());
        assert!(delta_h.is_finite());
    }

    assert_eq!(reactions.len(), REACTION_IDS.len());
    assert!(
        species.len() >= 12,
        "fixture should exercise a connected multispecies state"
    );
    assert_eq!(molar_masses.len(), species.len());
    assert!(
        thermo
            .values()
            .all(|(cp, h)| cp.is_finite() && h.is_finite())
    );
    assert!(reactions.iter().all(|reaction| reaction.a_si > 0.0));
    assert!(REACTION_IDS.iter().all(|id| id.parse::<usize>().is_ok()));
    assert!(
        reactions
            .iter()
            .any(|reaction| reaction.temperature_exponent != 0.0)
    );
    println!(
        "large mechanism preparation: reactions={} species={} thermo_sources={:?} T={} K",
        reactions.len(),
        species.len(),
        provenance
            .values()
            .collect::<std::collections::BTreeSet<_>>(),
        THERMO_TEMPERATURE_K
    );
}

#[test]
fn larger_cantera_subset_damped_fixed_grid_converges_but_violates_species_bounds() {
    let reactor = solve_damped(0, 33);
    let telemetry = reactor
        .solver
        .last_telemetry
        .as_ref()
        .expect("Damped telemetry retained");
    let RustedSciThe::numerical::BVP_api::BvpSolverTelemetrySnapshot::Damped(telemetry) = telemetry
    else {
        panic!("expected Damped telemetry snapshot");
    };
    assert_eq!(
        telemetry.last_solve_status,
        RustedSciThe::numerical::BVP_Damp::telemetry::BvpSolveStatus::Converged,
        "RST 0.4.32 should converge on the former fixed-grid failure"
    );
    assert!(telemetry.counters.iterations > 0);
    assert!(telemetry.counters.residual_requests > 0);
    assert!(telemetry.counters.jacobian_requests > 0);
    let metrics = published_metrics(&reactor);
    assert!((200.0..=4000.0).contains(&metrics.min_temperature_k));
    assert!((200.0..=4000.0).contains(&metrics.max_temperature_k));
    assert!(metrics.max_mass_fraction_sum_error < 1e-3);
    assert!(metrics.max_species_change > 1e-8);
    assert!(
        metrics.min_mass_fraction < -1e-5,
        "this fixed-grid case should continue to expose the known physical bounds gap"
    );
    println!(
        "large-mechanism Damped fixed grid: status={:?} reactions={} species={} states={} nodes={} T=[{:.3},{:.3}]K min_Y={:.3e} max_Y={:.3e} closure={:.3e} max_dY={:.3e} iterations={} damping_trials={} linear_solves={}",
        telemetry.last_solve_status,
        REACTION_IDS.len(),
        reactor.kindata.substances.len(),
        reactor.solver.unknowns.len(),
        metrics.nodes,
        metrics.min_temperature_k,
        metrics.max_temperature_k,
        metrics.min_mass_fraction,
        metrics.max_mass_fraction,
        metrics.max_mass_fraction_sum_error,
        metrics.max_species_change,
        telemetry.counters.iterations,
        telemetry.counters.damping_trials,
        telemetry.counters.linear_solves
    );
}

#[test]
#[ignore = "large mechanism adaptive-grid strategy diagnostic; run explicitly"]
fn larger_cantera_subset_damped_grid_strategy_matrix_reports_physical_quality() {
    let methods = [
        ("DoublePoints", GridRefinementMethod::DoublePoints),
        ("Pearson-1e-2", GridRefinementMethod::Pearson(1e-2, 1.4)),
        (
            "GrcarSmooke-1e-2-5e-2",
            GridRefinementMethod::GrcarSmooke(1e-2, 5e-2, 1.4),
        ),
        (
            "TwoPoint-2e-2-C1.4",
            GridRefinementMethod::TwoPoint(2e-2, 2e-2, 1.4),
        ),
    ];
    let method_count = methods.len();
    let reference = solve_bvp_sci_reference();
    let fixed = solve_damped(0, 33);
    let fixed_metrics = published_metrics(&fixed);
    let fixed_errors = profile_errors(&fixed, &reference);
    println!(
        "large-mechanism Damped fixed baseline nodes={} min_Y={:.3e} dT={:.3e} dY={:.3e} dJ={:.3e} dq={:.3e}",
        fixed_metrics.nodes,
        fixed_metrics.min_mass_fraction,
        fixed_errors[0],
        fixed_errors[1],
        fixed_errors[2],
        fixed_errors[3],
    );
    let mut rows = Vec::new();
    for (name, method) in methods {
        let reactor = solve_damped_with_grid(2, 33, method);
        let metrics = published_metrics(&reactor);
        assert!((200.0..=4000.0).contains(&metrics.min_temperature_k));
        assert!((200.0..=4000.0).contains(&metrics.max_temperature_k));
        assert!(metrics.max_mass_fraction_sum_error < 1e-3);
        assert!(
            metrics.min_mass_fraction >= -1e-12,
            "strategy {name} published a negative mass fraction: {:.3e}",
            metrics.min_mass_fraction
        );
        let telemetry = reactor
            .solver
            .last_telemetry
            .as_ref()
            .expect("Damped strategy telemetry retained");
        let RustedSciThe::numerical::BVP_api::BvpSolverTelemetrySnapshot::Damped(telemetry) =
            telemetry
        else {
            panic!("expected Damped telemetry snapshot");
        };
        assert_eq!(
            telemetry.last_solve_status,
            RustedSciThe::numerical::BVP_Damp::telemetry::BvpSolveStatus::Converged,
            "strategy {name} did not converge"
        );
        let errors = profile_errors(&reactor, &reference);
        for (field, (error, baseline)) in errors.iter().zip(fixed_errors.iter()).enumerate() {
            assert!(
                error < baseline,
                "strategy {name} did not improve profile field {field}: {error:.3e} >= {baseline:.3e}"
            );
        }
        println!(
            "large-mechanism Damped grid={name} nodes={} T=[{:.3},{:.3}]K min_Y={:.3e} max_Y={:.3e} closure={:.3e} dT={:.3e} dY={:.3e} dJ={:.3e} dq={:.3e} refinements={} iterations={} linear_solves={}",
            metrics.nodes,
            metrics.min_temperature_k,
            metrics.max_temperature_k,
            metrics.min_mass_fraction,
            metrics.max_mass_fraction,
            metrics.max_mass_fraction_sum_error,
            errors[0],
            errors[1],
            errors[2],
            errors[3],
            telemetry.counters.grid_refinements,
            telemetry.counters.iterations,
            telemetry.counters.linear_solves,
        );
        rows.push((name, metrics.min_mass_fraction, metrics.nodes, errors));
    }
    assert_eq!(rows.len(), method_count);
}

#[test]
#[ignore = "large multispecies CorrectionNewton comparison; run explicitly"]
fn larger_cantera_subset_correction_newton_mesh_diagnostic() {
    let mut successes = 0;
    for nodes in [33, 65, 129] {
        match solve_correction_newton(nodes) {
            Ok(reactor) => {
                successes += 1;
                let metrics = published_metrics(&reactor);
                assert_physical_profile(&metrics, "CorrectionNewton");
                let telemetry = reactor
                    .solver
                    .last_telemetry
                    .as_ref()
                    .expect("CorrectionNewton telemetry retained");
                let BvpSolverTelemetrySnapshot::CorrectionNewton(Some(telemetry)) = telemetry
                else {
                    panic!("expected CorrectionNewton telemetry snapshot");
                };
                assert!(telemetry.residual_evaluations > 0);
                assert!(telemetry.jacobian_evaluations > 0);
                println!(
                    "large-mechanism CorrectionNewton nodes={nodes} outcome=PASS published_nodes={} T=[{:.3},{:.3}]K min_Y={:.3e} max_Y={:.3e} closure={:.3e} residual_calls={} jacobian_calls={} linear_solves={}",
                    metrics.nodes,
                    metrics.min_temperature_k,
                    metrics.max_temperature_k,
                    metrics.min_mass_fraction,
                    metrics.max_mass_fraction,
                    metrics.max_mass_fraction_sum_error,
                    telemetry.residual_evaluations,
                    telemetry.jacobian_evaluations,
                    telemetry.linear_solves,
                );
            }
            Err(error) => {
                println!("large-mechanism CorrectionNewton nodes={nodes} outcome=FAIL {error}")
            }
        }
    }
    assert!(
        successes > 0,
        "CorrectionNewton failed on every tested mesh"
    );
}

#[test]
#[ignore = "large multispecies BVP_sci correctness story; run explicitly"]
fn larger_cantera_subset_bvp_sci_solves_the_same_physical_case() {
    let mut reactor = build_reactor();
    reactor.setup_bvp().expect("large-mechanism BVP setup");
    reactor.solver.arg_name = "z".to_string();
    reactor.solver.x_range = (0.0, 1.0);
    reactor
        .solver
        .solve_BVPsci(inlet_initial_guess(&reactor, 33), 33, 256, 1e-5)
        .expect("BVP_sci should solve the same multispecies fixture");
    reactor
        .postprocessing()
        .expect("publish dimensional BVP_sci result");

    let solution = reactor
        .solver
        .get_solution()
        .expect("BVP_sci published solution");
    let mesh = reactor.solver.x_mesh.as_ref().expect("published mesh");
    assert_eq!(solution.nrows(), mesh.len());
    assert_eq!(solution.ncols(), reactor.solver.unknowns.len());
    assert!(solution.iter().all(|value| value.is_finite()));

    let temperature_column = reactor
        .solver
        .unknowns
        .iter()
        .position(|name| name == "Teta")
        .expect("temperature state");
    let temperatures = solution.column(temperature_column);
    let min_temperature = temperatures.iter().copied().fold(f64::INFINITY, f64::min);
    let max_temperature = temperatures
        .iter()
        .copied()
        .fold(f64::NEG_INFINITY, f64::max);
    assert!(
        (200.0..=4000.0).contains(&min_temperature) && (200.0..=4000.0).contains(&max_temperature),
        "BVP_sci temperature range [{min_temperature}, {max_temperature}] K is not physical"
    );

    let mut max_species_sum_error = 0.0_f64;
    for node in 0..solution.nrows() {
        let mut sum = 0.0;
        for (species_index, name) in reactor.kindata.substances.iter().enumerate() {
            let fraction = solution[(node, 2 + 2 * species_index)];
            assert!(
                (-1e-5..=1.0 + 1e-5).contains(&fraction),
                "BVP_sci mass fraction {name} at node {node} is outside [0,1]: {fraction}"
            );
            sum += fraction;
        }
        max_species_sum_error = max_species_sum_error.max((sum - 1.0).abs());
    }
    assert!(
        max_species_sum_error < 1e-3,
        "BVP_sci mass-fraction closure error={max_species_sum_error}"
    );
    let max_species_change = reactor
        .kindata
        .substances
        .iter()
        .enumerate()
        .map(|(species_index, name)| {
            solution
                .column(2 + 2 * species_index)
                .iter()
                .map(|value| (value - reactor.boundary_condition[name]).abs())
                .fold(0.0, f64::max)
        })
        .fold(0.0, f64::max);
    assert!(
        max_species_change > 1e-8,
        "fixture must exhibit nontrivial reaction progress"
    );

    let telemetry = reactor
        .solver
        .last_telemetry
        .as_ref()
        .expect("BVP_sci telemetry retained");
    let RustedSciThe::numerical::BVP_api::BvpSolverTelemetrySnapshot::Sci(_telemetry) = telemetry
    else {
        panic!("expected BVP_sci telemetry snapshot");
    };
    println!(
        "large-mechanism BVP_sci: reactions={} species={} nodes={} T=[{min_temperature:.3}, {max_temperature:.3}] K max_species_change={max_species_change:.3e} mass_fraction_sum_error={max_species_sum_error:.3e}",
        REACTION_IDS.len(),
        reactor.kindata.substances.len(),
        mesh.len(),
    );
}

#[test]
#[ignore = "large multispecies Frozen dense-grid correctness story; run explicitly"]
fn larger_cantera_subset_frozen_dense_65_node_diagnostic() {
    const FROZEN_NODES: usize = 65;
    let mut reactor = build_reactor();
    reactor.setup_bvp().expect("large-mechanism BVP setup");
    reactor.solver.arg_name = "z".to_string();
    reactor.solver.x_range = (0.0, 1.0);
    reactor
        .solver
        .solve_BVPfrozen(
            inlet_initial_guess(&reactor, FROZEN_NODES),
            FrozenSolverOptions::dense_naive()
                .with_tolerance(1e-5)
                .with_max_iterations(500),
        )
        .expect("Frozen should solve the same fixture on a dense fixed grid");
    reactor
        .postprocessing()
        .expect("publish dimensional Frozen result");

    let solution = reactor
        .solver
        .get_solution()
        .expect("Frozen published solution");
    let mesh = reactor.solver.x_mesh.as_ref().expect("published mesh");
    assert_eq!(mesh.len(), FROZEN_NODES);
    assert_eq!(solution.nrows(), mesh.len());
    assert_eq!(solution.ncols(), reactor.solver.unknowns.len());
    assert!(solution.iter().all(|value| value.is_finite()));

    let temperature_column = reactor
        .solver
        .unknowns
        .iter()
        .position(|name| name == "Teta")
        .expect("temperature state");
    let temperatures = solution.column(temperature_column);
    let min_temperature = temperatures.iter().copied().fold(f64::INFINITY, f64::min);
    let max_temperature = temperatures
        .iter()
        .copied()
        .fold(f64::NEG_INFINITY, f64::max);
    assert!(
        (200.0..=4000.0).contains(&min_temperature) && (200.0..=4000.0).contains(&max_temperature),
        "Frozen temperature range [{min_temperature}, {max_temperature}] K is not physical"
    );

    let mut max_species_sum_error = 0.0_f64;
    let mut max_species_change = 0.0_f64;
    let mut min_mass_fraction = f64::INFINITY;
    for node in 0..solution.nrows() {
        let mut sum = 0.0;
        for (species_index, name) in reactor.kindata.substances.iter().enumerate() {
            let fraction = solution[(node, 2 + 2 * species_index)];
            min_mass_fraction = min_mass_fraction.min(fraction);
            sum += fraction;
            max_species_change =
                max_species_change.max((fraction - reactor.boundary_condition[name]).abs());
        }
        max_species_sum_error = max_species_sum_error.max((sum - 1.0).abs());
    }
    assert!(max_species_sum_error < 1e-3);
    assert!(max_species_change > 1e-8);
    println!(
        "large-mechanism Frozen: reactions={} species={} nodes={} T=[{min_temperature:.3}, {max_temperature:.3}] K min_mass_fraction={min_mass_fraction:.3e} max_species_change={max_species_change:.3e} mass_fraction_sum_error={max_species_sum_error:.3e}",
        REACTION_IDS.len(),
        reactor.kindata.substances.len(),
        mesh.len()
    );
}

#[test]
#[ignore = "large multispecies Frozen banded-grid correctness story; run explicitly"]
fn larger_cantera_subset_frozen_banded_257_publishes_physical_profile() {
    const FROZEN_NODES: usize = 257;
    let mut reactor = build_reactor();
    reactor.setup_bvp().expect("large-mechanism BVP setup");
    reactor.solver.arg_name = "z".to_string();
    reactor.solver.x_range = (0.0, 1.0);
    reactor.solver = reactor
        .solver
        .with_bvp_telemetry_mode(BvpApiTelemetryMode::Counters);
    reactor
        .solver
        .solve_BVPfrozen(
            inlet_initial_guess(&reactor, FROZEN_NODES),
            FrozenSolverOptions::banded_frozen()
                .with_tolerance(1e-5)
                .with_max_iterations(500),
        )
        .expect("257-node banded Frozen should publish its converged profile");
    reactor
        .postprocessing()
        .expect("publish dimensional banded Frozen profile");

    let telemetry = reactor
        .solver
        .last_telemetry
        .as_ref()
        .expect("Frozen telemetry retained");
    let BvpSolverTelemetrySnapshot::Frozen(telemetry) = telemetry else {
        panic!("expected Frozen telemetry snapshot");
    };
    assert_eq!(
        telemetry.last_solve_status,
        RustedSciThe::numerical::BVP_Damp::telemetry::BvpSolveStatus::Converged
    );
    assert!(telemetry.counters.iterations > 0);

    let metrics = published_metrics(&reactor);
    assert_eq!(metrics.nodes, FROZEN_NODES);
    assert_physical_profile(&metrics, "257-node banded Frozen");
    println!(
        "large-mechanism Frozen banded: status={:?} reactions={} species={} nodes={} T=[{:.3},{:.3}]K min_Y={:.3e} max_Y={:.3e} closure={:.3e} max_dY={:.3e} iterations={}",
        telemetry.last_solve_status,
        REACTION_IDS.len(),
        reactor.kindata.substances.len(),
        metrics.nodes,
        metrics.min_temperature_k,
        metrics.max_temperature_k,
        metrics.min_mass_fraction,
        metrics.max_mass_fraction,
        metrics.max_mass_fraction_sum_error,
        metrics.max_species_change,
        telemetry.counters.iterations,
    );
}
