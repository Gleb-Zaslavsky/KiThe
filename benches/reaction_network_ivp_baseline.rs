//! Criterion baseline for direct reaction-network IVP preparation and solves.
//!
//! The diagnostic table is a separate Detailed-telemetry pass. Criterion
//! samples use telemetry Off and split cold preparation from prepared reruns.

use std::time::Instant;

use KiThe::Kinetics::reaction_network_ivp::{
    NetworkIvpOptions, NetworkProblem, PreparedNetworkIvp, RateCoefficient, RateLaw, ReactionSpec,
    SpeciesSpec, StoichTerm, TemperatureProgram, solve_network,
};
use RustedSciThe::numerical::SolverTelemetryMode;
use RustedSciThe::numerical::StiffAPI::{
    StiffCallbackPolicy, StiffFrontend, StiffMatrixLayout, StiffSolverKind,
};
use criterion::{BenchmarkId, Criterion, black_box, criterion_group, criterion_main};
use tabled::settings::Style;
use tabled::{Table, Tabled};

const SIZES: [usize; 3] = [3, 15, 40];
const FAMILIES: [StiffSolverKind; 4] = [
    StiffSolverKind::Bdf,
    StiffSolverKind::Radau,
    StiffSolverKind::BackwardEuler,
    StiffSolverKind::Lsode2,
];
const FRONTENDS: [StiffFrontend; 2] = [StiffFrontend::ExprLegacy, StiffFrontend::AtomViewNative];
const LAYOUTS: [StiffMatrixLayout; 3] = [
    StiffMatrixLayout::Dense,
    StiffMatrixLayout::Sparse,
    StiffMatrixLayout::BandedAuto,
];

#[derive(Tabled)]
struct StageRow {
    states: usize,
    reactions: usize,
    family: String,
    frontend: String,
    layout: String,
    callback_policy: String,
    status: &'static str,
    prepared_wall_ms: String,
    assembly_ms: String,
    preparation_ms: String,
    integration_ms: String,
    output_ms: String,
    unattributed_ms: String,
}

fn chain_problem(size: usize) -> NetworkProblem {
    let species = (0..size)
        .map(|index| SpeciesSpec {
            id: format!("S{index}"),
            initial_concentration: if index == 0 { 1.0 } else { 0.0 },
        })
        .collect();
    let reactions = (0..size - 1)
        .map(|index| ReactionSpec {
            id: format!("R{index}"),
            reactants: vec![StoichTerm::new(format!("S{index}"), 1.0)],
            products: vec![StoichTerm::new(format!("S{}", index + 1), 1.0)],
            rate_law: RateLaw::MassAction {
                coefficient: RateCoefficient::Constant(0.3),
            },
        })
        .collect();
    NetworkProblem::new(
        species,
        reactions,
        TemperatureProgram::Constant(700.0),
        0.0,
        2.0,
    )
    .expect("linear chain fixture is valid")
}

fn options(
    family: StiffSolverKind,
    frontend: StiffFrontend,
    telemetry: SolverTelemetryMode,
    matrix_layout: StiffMatrixLayout,
    callback_policy: Option<StiffCallbackPolicy>,
) -> NetworkIvpOptions {
    NetworkIvpOptions {
        family,
        frontend,
        matrix_layout,
        callback_policy,
        telemetry,
        max_step: 0.02,
        ..NetworkIvpOptions::default()
    }
}

fn ms(value: Option<std::time::Duration>) -> f64 {
    value.map_or(0.0, |duration| duration.as_secs_f64() * 1_000.0)
}

fn print_stage_matrix() {
    let mut rows = Vec::new();
    for size in SIZES {
        let problem = chain_problem(size);
        for family in FAMILIES {
            for frontend in FRONTENDS {
                let diagnostic_options = options(
                    family,
                    frontend,
                    SolverTelemetryMode::Detailed,
                    StiffMatrixLayout::Dense,
                    None,
                );
                let started = Instant::now();
                let mut prepared = match PreparedNetworkIvp::new(&problem, &diagnostic_options) {
                    Ok(prepared) => prepared,
                    Err(error) => {
                        rows.push(StageRow {
                            states: size,
                            reactions: size - 1,
                            family: format!("{family:?}"),
                            frontend: format!("{frontend:?}"),
                            layout: "Dense".into(),
                            callback_policy: "default".into(),
                            status: "PREP_FAIL",
                            prepared_wall_ms: String::new(),
                            assembly_ms: String::new(),
                            preparation_ms: String::new(),
                            integration_ms: String::new(),
                            output_ms: String::new(),
                            unattributed_ms: String::new(),
                        });
                        eprintln!("stage probe {family:?}/{size}: {error}");
                        continue;
                    }
                };
                match prepared.solve() {
                    Ok(_) => {
                        let wall_ms = started.elapsed().as_secs_f64() * 1_000.0;
                        let stages = prepared
                            .stage_diagnostics()
                            .expect("Detailed mode records stage data");
                        let assembly_ms = ms(stages.equation_assembly.elapsed);
                        let preparation_ms = ms(stages.backend_preparation.elapsed);
                        let integration_ms = ms(stages.integration.elapsed);
                        let output_ms = ms(stages.output_postprocessing.elapsed);
                        let attributed = assembly_ms + preparation_ms + integration_ms + output_ms;
                        rows.push(StageRow {
                            states: size,
                            reactions: size - 1,
                            family: format!("{family:?}"),
                            frontend: format!("{frontend:?}"),
                            layout: "Dense".into(),
                            callback_policy: "default".into(),
                            status: "PASS",
                            prepared_wall_ms: format!("{wall_ms:.3}"),
                            assembly_ms: format!("{assembly_ms:.3}"),
                            preparation_ms: format!("{preparation_ms:.3}"),
                            integration_ms: format!("{integration_ms:.3}"),
                            output_ms: format!("{output_ms:.3}"),
                            unattributed_ms: format!("{:.3}", wall_ms - attributed),
                        });
                    }
                    Err(error) => {
                        rows.push(StageRow {
                            states: size,
                            reactions: size - 1,
                            family: format!("{family:?}"),
                            frontend: format!("{frontend:?}"),
                            layout: "Dense".into(),
                            callback_policy: "default".into(),
                            status: "SOLVE_FAIL",
                            prepared_wall_ms: format!(
                                "{:.3}",
                                started.elapsed().as_secs_f64() * 1_000.0
                            ),
                            assembly_ms: String::new(),
                            preparation_ms: String::new(),
                            integration_ms: String::new(),
                            output_ms: String::new(),
                            unattributed_ms: String::new(),
                        });
                        eprintln!("stage probe {family:?}/{size}: {error}");
                    }
                }
            }
        }
    }
    println!(
        "\nReaction-network IVP baseline stages (Detailed diagnostic pass; separate from Criterion)\n{}",
        Table::new(rows).with(Style::markdown())
    );
}

fn policies_for(family: StiffSolverKind) -> Vec<Option<StiffCallbackPolicy>> {
    if family == StiffSolverKind::BackwardEuler {
        return vec![None];
    }
    vec![
        Some(StiffCallbackPolicy::Sequential),
        Some(StiffCallbackPolicy::Parallel { min_work: 1 }),
        Some(StiffCallbackPolicy::Auto { min_work: 1 }),
    ]
}

fn print_large_configuration_matrix() {
    let size = *SIZES.last().unwrap();
    let problem = chain_problem(size);
    let mut rows = Vec::new();
    for family in FAMILIES {
        for frontend in FRONTENDS {
            for layout in LAYOUTS {
                for callback_policy in policies_for(family) {
                    let label = format!("{layout:?}/{callback_policy:?}");
                    let diagnostic_options = options(
                        family,
                        frontend,
                        SolverTelemetryMode::Detailed,
                        layout,
                        callback_policy,
                    );
                    let started = Instant::now();
                    let mut prepared = match PreparedNetworkIvp::new(&problem, &diagnostic_options)
                    {
                        Ok(prepared) => prepared,
                        Err(error) => {
                            rows.push(StageRow {
                                states: size,
                                reactions: size - 1,
                                family: format!("{family:?}"),
                                frontend: format!("{frontend:?}"),
                                layout: format!("{layout:?}"),
                                callback_policy: format!("{callback_policy:?}"),
                                status: "UNSUPPORTED",
                                prepared_wall_ms: String::new(),
                                assembly_ms: String::new(),
                                preparation_ms: String::new(),
                                integration_ms: String::new(),
                                output_ms: String::new(),
                                unattributed_ms: String::new(),
                            });
                            eprintln!(
                                "large-config preflight {label}/{family:?}/{frontend:?}: {error}"
                            );
                            continue;
                        }
                    };
                    match prepared.solve() {
                        Ok(_) => {
                            let wall_ms = started.elapsed().as_secs_f64() * 1_000.0;
                            let stages = prepared
                                .stage_diagnostics()
                                .expect("Detailed mode records stage data");
                            let assembly_ms = ms(stages.equation_assembly.elapsed);
                            let preparation_ms = ms(stages.backend_preparation.elapsed);
                            let integration_ms = ms(stages.integration.elapsed);
                            let output_ms = ms(stages.output_postprocessing.elapsed);
                            let attributed =
                                assembly_ms + preparation_ms + integration_ms + output_ms;
                            rows.push(StageRow {
                                states: size,
                                reactions: size - 1,
                                family: format!("{family:?}"),
                                frontend: format!("{frontend:?}"),
                                layout: format!("{layout:?}"),
                                callback_policy: format!("{callback_policy:?}"),
                                status: "PASS",
                                prepared_wall_ms: format!("{wall_ms:.3}"),
                                assembly_ms: format!("{assembly_ms:.3}"),
                                preparation_ms: format!("{preparation_ms:.3}"),
                                integration_ms: format!("{integration_ms:.3}"),
                                output_ms: format!("{output_ms:.3}"),
                                unattributed_ms: format!("{:.3}", wall_ms - attributed),
                            });
                        }
                        Err(error) => {
                            rows.push(StageRow {
                                states: size,
                                reactions: size - 1,
                                family: format!("{family:?}"),
                                frontend: format!("{frontend:?}"),
                                layout: format!("{layout:?}"),
                                callback_policy: format!("{callback_policy:?}"),
                                status: "SOLVE_FAIL",
                                prepared_wall_ms: format!(
                                    "{:.3}",
                                    started.elapsed().as_secs_f64() * 1_000.0
                                ),
                                assembly_ms: String::new(),
                                preparation_ms: String::new(),
                                integration_ms: String::new(),
                                output_ms: String::new(),
                                unattributed_ms: String::new(),
                            });
                            eprintln!(
                                "large-config solve {label}/{family:?}/{frontend:?}: {error}"
                            );
                        }
                    }
                }
            }
        }
    }
    println!(
        "\nReaction-network IVP 40-state matrix/layout/callback stages (unsupported native combinations are explicit)\n{}",
        Table::new(rows).with(Style::markdown())
    );
}

fn reaction_network_ivp_baseline(criterion: &mut Criterion) {
    print_stage_matrix();
    print_large_configuration_matrix();
    let mut group = criterion.benchmark_group("reaction_network_ivp");
    group.sample_size(20);
    for size in SIZES {
        let problem = chain_problem(size);
        for family in FAMILIES {
            for frontend in FRONTENDS {
                let solver_options = options(
                    family,
                    frontend,
                    SolverTelemetryMode::Off,
                    StiffMatrixLayout::Dense,
                    None,
                );
                let label = format!("{frontend:?}/{family:?}/{size}-states");
                let diagnostic = solve_network(
                    &problem,
                    &options(
                        family,
                        frontend,
                        SolverTelemetryMode::Off,
                        StiffMatrixLayout::Dense,
                        None,
                    ),
                );
                if let Err(error) = diagnostic {
                    eprintln!("Skipping Criterion {label}: {error}");
                    continue;
                }
                group.bench_with_input(
                    BenchmarkId::new("cold_prepare_and_solve", &label),
                    &label,
                    |bench, _| {
                        bench.iter(|| {
                            black_box(
                                solve_network(&problem, &solver_options)
                                    .expect("case passed preflight"),
                            )
                        })
                    },
                );

                let mut prepared = PreparedNetworkIvp::new(&problem, &solver_options)
                    .expect("case passed preflight");
                prepared.solve().expect("case passed preflight");
                group.bench_function(BenchmarkId::new("prepared_rerun", &label), |bench| {
                    bench.iter(|| {
                        let endpoint = prepared
                            .solve()
                            .expect("prepared rerun succeeds")
                            .concentrations
                            .last()
                            .expect("trajectory contains an endpoint")[0];
                        black_box(endpoint)
                    })
                });
            }
        }
    }
    let size = *SIZES.last().unwrap();
    let problem = chain_problem(size);
    for family in FAMILIES {
        for frontend in FRONTENDS {
            for layout in LAYOUTS {
                for callback_policy in policies_for(family) {
                    let solver_options = options(
                        family,
                        frontend,
                        SolverTelemetryMode::Off,
                        layout,
                        callback_policy,
                    );
                    let label = format!(
                        "{frontend:?}/{family:?}/{layout:?}/{callback_policy:?}/{size}-states"
                    );
                    if let Err(error) = solve_network(&problem, &solver_options) {
                        eprintln!("Skipping Criterion {label}: {error}");
                        continue;
                    }
                    group.bench_with_input(
                        BenchmarkId::new("large_matrix_callback_policy", &label),
                        &label,
                        |bench, _| {
                            bench.iter(|| {
                                black_box(
                                    solve_network(&problem, &solver_options)
                                        .expect("case passed preflight"),
                                )
                            })
                        },
                    );
                }
            }
        }
    }
    group.finish();
}

criterion_group!(benches, reaction_network_ivp_baseline);
criterion_main!(benches);
