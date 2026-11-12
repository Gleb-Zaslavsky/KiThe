//! Damped adaptive-grid strategy matrix, isolated from branch and SciPy probes.

use super::*;

struct Strategy {
    name: &'static str,
    method: GridRefinementMethod,
}

fn strategies() -> [Strategy; 4] {
    [
        Strategy {
            name: "DoublePoints",
            method: GridRefinementMethod::DoublePoints,
        },
        Strategy {
            name: "Pearson(1e-2,1.4)",
            method: GridRefinementMethod::Pearson(1e-2, 1.4),
        },
        Strategy {
            name: "Grcar(1e-2,5e-2,1.4)",
            method: GridRefinementMethod::GrcarSmooke(1e-2, 5e-2, 1.4),
        },
        Strategy {
            name: "TwoPoint(2e-2,2e-2,1.4)",
            method: GridRefinementMethod::TwoPoint(2e-2, 2e-2, 1.4),
        },
    ]
}

#[test]
fn mf1963_damped_adaptive_strategy_matrix_reports_every_case() {
    let mut rows = 0;
    let mut outcomes = Vec::with_capacity(8);
    println!("MF-04A Damped adaptive matrix (diagnostic; failures do not short-circuit)");
    println!(
        "case | strategy | outcome | nodes | refinements | max_dT_K | max_dY | rel_dq | closure | defects"
    );

    for case in cases() {
        let initial_nodes = if case.name == "MF-LOW" { 2049 } else { 301 };
        for strategy in strategies() {
            rows += 1;
            let mut template = reactor_for(case);
            if let Err(error) = template.setup_bvp() {
                outcomes.push((case.name, strategy.name, "SETUP_FAIL"));
                println!(
                    "{} | {} | SETUP_FAIL | - | - | - | - | - | - | {error}",
                    case.name, strategy.name,
                );
                continue;
            }
            let seed = damped_reduced_seed_from_nondimensional_profile(
                &template,
                &reference_seed(case, initial_nodes),
                initial_nodes,
            );
            let (result, telemetry) = solve_mf_damped_with_seed_and_telemetry(
                case,
                seed,
                initial_nodes,
                2,
                strategy.method,
                50,
                true,
                false,
            );
            match result {
                Ok(reactor) => {
                    let solution = reactor.solver.get_solution().expect("published solution");
                    let errors = production_profile_errors(case, &reactor);
                    let defects = mf_discrete_trapezoid_defects(case, &reactor);
                    let refinements = match reactor.solver.last_telemetry.as_ref() {
                        Some(BvpSolverTelemetrySnapshot::Damped(snapshot)) => {
                            snapshot.counters.grid_refinements
                        }
                        other => panic!("expected Damped telemetry, got {other:?}"),
                    };
                    let physical = solution.iter().all(|value| value.is_finite())
                        && solution
                            .column(0)
                            .iter()
                            .all(|theta| T0 + T_SCALE * theta > 0.0)
                        && solution
                            .column(2)
                            .iter()
                            .chain(solution.column(4).iter())
                            .all(|fraction| (0.0..=1.0).contains(fraction));
                    let accepted = physical
                        && errors.0 < 0.5
                        && errors.1 < 2e-3
                        && errors.2 < 2e-3
                        && errors.3 < 2e-8
                        && defects.0 < 1e-7
                        && defects.1 < 1e-7
                        && defects.2 < 1e-8;
                    println!(
                        "{} | {} | {} | {} | {} | {:.3e} | {:.3e} | {:.3e} | {:.3e} | ({:.2e},{:.2e},{:.2e})",
                        case.name,
                        strategy.name,
                        if accepted { "PASS" } else { "ACCURACY_FAIL" },
                        solution.nrows(),
                        refinements,
                        errors.0,
                        errors.1,
                        errors.2,
                        errors.3,
                        defects.0,
                        defects.1,
                        defects.2,
                    );
                    outcomes.push((
                        case.name,
                        strategy.name,
                        if accepted { "PASS" } else { "ACCURACY_FAIL" },
                    ));
                }
                Err(error) => {
                    let reason = error.split("; telemetry=").next().unwrap_or(&error);
                    let diagnostic = match telemetry {
                        Some(BvpSolverTelemetrySnapshot::Damped(snapshot)) => format!(
                            "status={:?} iter={} residual_calls={} jacobian={} trials={} rejected={} refinements={} prep_ms={:.1}",
                            snapshot.last_solve_status,
                            snapshot.counters.iterations,
                            snapshot.counters.residual_calls,
                            snapshot.counters.jacobian_requests,
                            snapshot.counters.damping_trials,
                            snapshot.counters.damping_rejections,
                            snapshot.counters.grid_refinements,
                            snapshot
                                .generation
                                .as_ref()
                                .map_or(0.0, |g| g.total.as_secs_f64() * 1e3),
                        ),
                        other => format!("telemetry={other:?}"),
                    };
                    println!(
                        "{} | {} | SOLVE_FAIL | - | - | - | - | - | - | {reason}; {diagnostic}",
                        case.name, strategy.name,
                    );
                    outcomes.push((case.name, strategy.name, "SOLVE_FAIL"));
                }
            }
        }
    }
    assert_eq!(rows, 8, "every regime/strategy pair must be attempted");
    let expected = [
        ("MF-LOW", "DoublePoints", "SOLVE_FAIL"),
        ("MF-LOW", "Pearson(1e-2,1.4)", "SOLVE_FAIL"),
        ("MF-LOW", "Grcar(1e-2,5e-2,1.4)", "SOLVE_FAIL"),
        ("MF-LOW", "TwoPoint(2e-2,2e-2,1.4)", "SOLVE_FAIL"),
        ("MF-HIGH", "DoublePoints", "ACCURACY_FAIL"),
        ("MF-HIGH", "Pearson(1e-2,1.4)", "PASS"),
        ("MF-HIGH", "Grcar(1e-2,5e-2,1.4)", "PASS"),
        ("MF-HIGH", "TwoPoint(2e-2,2e-2,1.4)", "ACCURACY_FAIL"),
    ];
    assert_eq!(
        outcomes, expected,
        "refresh MF-04A findings when outcomes change"
    );
}
