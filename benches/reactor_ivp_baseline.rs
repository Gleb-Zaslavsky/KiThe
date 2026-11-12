//! PERF-IVP-01: matched family/frontend/layout baseline for bounded condensed IVPs.
//!
//! Criterion runs with telemetry disabled and measures task construction through
//! accepted solution publication. A separate detailed pass prints KiThe stages
//! and RST-native measurements; those instrumented values are not mixed into
//! Criterion's wall-time samples.
//!
//! Run with:
//! cargo bench --bench reactor_ivp_baseline --no-default-features -j 2 -- --noplot
//! Routes: LSODE2/BDF Dense; Radau Dense/Sparse/Banded, each ExprLegacy/AtomView.

use std::collections::HashMap;
use std::hint::black_box;
use std::time::Instant;

use criterion::{BenchmarkId, Criterion, criterion_group, criterion_main};
use nalgebra::DVector;
use tabled::settings::Style;
use tabled::{Table, Tabled};

use KiThe::Kinetics::mechfinder_api::ReactionData;
use KiThe::SimpleReactorsBVP::reactor_BVP_utils::ScalingConfig;
use KiThe::SimpleReactorsIVP::SimpleReactorIVP::{IvpSolveSnapshot, SimpleReactorTask};
use KiThe::SimpleReactorsIVP::solver_backend::{
    ReactorIvpCallbackPolicy, ReactorIvpExecutionBackend, ReactorIvpMatrixBackend,
    ReactorIvpMethod, ReactorIvpSolverConfig, ReactorIvpSolverFamily, ReactorIvpSymbolicBackend,
};
use RustedSciThe::numerical::LSODE2::{
    Lsode2AotProfile, Lsode2AotToolchain, Lsode2BackendConfig, solver::Lsode2Solver,
};
use RustedSciThe::numerical::SolverTelemetryMode;
use RustedSciThe::symbolic::ivp_telemetry::IvpTelemetry;

const DIAGNOSTIC_REPETITIONS: usize = 3;

#[derive(Clone, Copy, Debug)]
struct Route {
    family: ReactorIvpSolverFamily,
    frontend: ReactorIvpSymbolicBackend,
    matrix: ReactorIvpMatrixBackend,
}

#[derive(Clone, Copy, Debug)]
enum Fixture {
    Baseline,
    WarmInlet,
}

impl Fixture {
    fn label(self) -> &'static str {
        match self {
            Self::Baseline => "baseline-520K",
            Self::WarmInlet => "warm-inlet-900K",
        }
    }
}

impl Route {
    fn label(self) -> String {
        format!("{:?}/{:?}/{:?}", self.family, self.frontend, self.matrix)
    }

    fn config(self, telemetry: SolverTelemetryMode) -> ReactorIvpSolverConfig {
        ReactorIvpSolverConfig::default()
            .with_family(self.family)
            .with_symbolic_backend(self.frontend)
            .with_matrix_backend(self.matrix)
            .with_integration_domain(0.0, 0.02)
            .with_max_step(1.0e-3)
            .with_rtol(1.0e-6)
            .with_atol(1.0e-8)
            .with_telemetry(telemetry)
    }
}

const ROUTES: [Route; 10] = [
    Route {
        family: ReactorIvpSolverFamily::Lsode2,
        frontend: ReactorIvpSymbolicBackend::ExprLegacy,
        matrix: ReactorIvpMatrixBackend::Dense,
    },
    Route {
        family: ReactorIvpSolverFamily::Lsode2,
        frontend: ReactorIvpSymbolicBackend::AtomView,
        matrix: ReactorIvpMatrixBackend::Dense,
    },
    Route {
        family: ReactorIvpSolverFamily::Bdf,
        frontend: ReactorIvpSymbolicBackend::ExprLegacy,
        matrix: ReactorIvpMatrixBackend::Dense,
    },
    Route {
        family: ReactorIvpSolverFamily::Bdf,
        frontend: ReactorIvpSymbolicBackend::AtomView,
        matrix: ReactorIvpMatrixBackend::Dense,
    },
    Route {
        family: ReactorIvpSolverFamily::Radau,
        frontend: ReactorIvpSymbolicBackend::ExprLegacy,
        matrix: ReactorIvpMatrixBackend::Dense,
    },
    Route {
        family: ReactorIvpSolverFamily::Radau,
        frontend: ReactorIvpSymbolicBackend::AtomView,
        matrix: ReactorIvpMatrixBackend::Dense,
    },
    Route {
        family: ReactorIvpSolverFamily::Radau,
        frontend: ReactorIvpSymbolicBackend::ExprLegacy,
        matrix: ReactorIvpMatrixBackend::Sparse,
    },
    Route {
        family: ReactorIvpSolverFamily::Radau,
        frontend: ReactorIvpSymbolicBackend::AtomView,
        matrix: ReactorIvpMatrixBackend::Sparse,
    },
    Route {
        family: ReactorIvpSolverFamily::Radau,
        frontend: ReactorIvpSymbolicBackend::ExprLegacy,
        matrix: ReactorIvpMatrixBackend::Banded,
    },
    Route {
        family: ReactorIvpSolverFamily::Radau,
        frontend: ReactorIvpSymbolicBackend::AtomView,
        matrix: ReactorIvpMatrixBackend::Banded,
    },
];

#[derive(Tabled)]
struct MacroRow {
    route: String,
    status: String,
    states: usize,
    output_points: usize,
    solve_call_wall_ms: String,
    workflow_ms: String,
}

#[derive(Tabled)]
struct StageRow {
    route: String,
    stage: String,
    elapsed_ms: String,
    workflow_share_pct: String,
    calls: String,
}

#[derive(Tabled)]
struct NativeRow {
    route: String,
    metric: String,
    value: String,
}

#[derive(Tabled)]
struct LifecycleRow {
    route: String,
    lifecycle: String,
    preparation_ms: String,
    first_solve_ms: String,
    repeated_solve_median_ms: String,
    repeats: usize,
}

#[derive(Tabled)]
struct ExecutionRow {
    fixture: String,
    solver_mode: String,
    frontend: String,
    execution: String,
    setup_ms: String,
    problem_build_ms: String,
    rst_preparation_ms: String,
    integration_ms: String,
    total_ms: String,
    residual_calls: u64,
    jacobian_calls: u64,
}

#[derive(Clone, Copy, Debug)]
enum CallbackPolicy {
    Sequential,
    Auto,
    ForcedParallel,
}

impl CallbackPolicy {
    fn label(self) -> &'static str {
        match self {
            Self::Sequential => "Sequential",
            Self::Auto => "Auto",
            Self::ForcedParallel => "Forced Parallel",
        }
    }

    fn config(self) -> ReactorIvpCallbackPolicy {
        match self {
            Self::Sequential => ReactorIvpCallbackPolicy::Sequential,
            Self::Auto => ReactorIvpCallbackPolicy::Auto { min_work: 0 },
            Self::ForcedParallel => ReactorIvpCallbackPolicy::Parallel { min_work: 0 },
        }
    }
}

#[derive(Tabled)]
struct CallbackPolicyRow {
    solver: String,
    policy: String,
    median_workflow_ms: String,
    callback_ms: String,
    parallel_dispatches: String,
    sequential_dispatches: String,
    max_endpoint_delta: String,
}

fn build_task(route: Route, telemetry: SolverTelemetryMode) -> SimpleReactorTask {
    build_task_for(route, telemetry, Fixture::Baseline)
}

fn build_task_for(
    route: Route,
    telemetry: SolverTelemetryMode,
    fixture: Fixture,
) -> SimpleReactorTask {
    let mut task = SimpleReactorTask::new();
    task.kindata
        .set_reaction_data_directly(
            vec![ReactionData::new_elementary(
                "H2=>2H".to_string(),
                vec![1.0e-3, 0.0, 8.0e3],
                None,
            )],
            None,
        )
        .expect("benchmark reaction should parse");
    task.thermal_effects = vec![0.0];
    task.ro = 1350.0;
    task.Cp = 1100.0;
    task.Lambda = 0.22;
    task.m = 0.012;
    task.scaling = ScalingConfig::new(80.0, 0.04, 90.0);
    let (temperature, h2, h) = match fixture {
        Fixture::Baseline => (520.0, 0.85, 0.15),
        Fixture::WarmInlet => (900.0, 0.70, 0.30),
    };
    task.initial_conditions = HashMap::from([
        ("T".to_string(), temperature),
        ("q".to_string(), 0.0),
        ("H2".to_string(), h2),
        ("H".to_string(), h),
    ]);
    task.set_solver_backend_config(route.config(telemetry));
    task
}

fn solve_fixture_once(
    route: Route,
    telemetry: SolverTelemetryMode,
    fixture: Fixture,
) -> (IvpSolveSnapshot, f64) {
    solve_fixture_once_with_policy(route, telemetry, fixture, None)
}

fn solve_fixture_once_with_policy(
    route: Route,
    telemetry: SolverTelemetryMode,
    fixture: Fixture,
    policy: Option<CallbackPolicy>,
) -> (IvpSolveSnapshot, f64) {
    let mut task = build_task_for(route, telemetry, fixture);
    if let Some(policy) = policy {
        task.set_solver_backend_config(
            route
                .config(telemetry)
                .with_callback_policy(policy.config()),
        );
    }
    let started = Instant::now();
    let snapshot = task
        .solve()
        .unwrap_or_else(|error| panic!("{} / {} failed: {error}", route.label(), fixture.label()));
    let elapsed_ms = started.elapsed().as_secs_f64() * 1_000.0;
    assert!(
        snapshot.status().starts_with("finished"),
        "{}: {}",
        route.label(),
        snapshot.status()
    );
    assert!(snapshot.values.iter().all(|value| value.is_finite()));
    (snapshot, elapsed_ms)
}

fn report_callback_policy_axis() -> Vec<CallbackPolicyRow> {
    let mut rows = Vec::new();
    for family in [
        ReactorIvpSolverFamily::Lsode2,
        ReactorIvpSolverFamily::Bdf,
        ReactorIvpSolverFamily::Radau,
    ] {
        let route = Route {
            family,
            frontend: ReactorIvpSymbolicBackend::ExprLegacy,
            matrix: ReactorIvpMatrixBackend::Dense,
        };
        let mut policies = Vec::new();
        for policy in [
            CallbackPolicy::Sequential,
            CallbackPolicy::Auto,
            CallbackPolicy::ForcedParallel,
        ] {
            let samples = (0..DIAGNOSTIC_REPETITIONS)
                .map(|_| {
                    solve_fixture_once_with_policy(
                        route,
                        SolverTelemetryMode::Detailed,
                        Fixture::WarmInlet,
                        Some(policy),
                    )
                })
                .collect::<Vec<_>>();
            policies.push((policy, samples));
        }
        let sequential = &policies[0].1;
        let mut sequential_workflow_ms = sequential
            .iter()
            .map(|(snapshot, _)| workflow_elapsed_ms(snapshot))
            .collect::<Vec<_>>();
        let sequential_median_ms = median(&mut sequential_workflow_ms);
        let sequential_representative = &sequential
            .iter()
            .min_by(|(left, _), (right, _)| {
                (workflow_elapsed_ms(left) - sequential_median_ms)
                    .abs()
                    .total_cmp(&(workflow_elapsed_ms(right) - sequential_median_ms).abs())
            })
            .expect("sequential diagnostic samples should exist")
            .0;
        let reference = sequential_representative
            .values
            .row(sequential_representative.values.nrows() - 1)
            .iter()
            .copied()
            .collect::<Vec<_>>();
        for (policy, samples) in &policies {
            let mut workflow_samples = samples
                .iter()
                .map(|(snapshot, _)| workflow_elapsed_ms(snapshot))
                .collect::<Vec<_>>();
            let workflow_ms = median(&mut workflow_samples);
            let representative = &samples
                .iter()
                .min_by(|(left, _), (right, _)| {
                    (workflow_elapsed_ms(left) - workflow_ms)
                        .abs()
                        .total_cmp(&(workflow_elapsed_ms(right) - workflow_ms).abs())
                })
                .expect("callback policy diagnostic samples should exist")
                .0;
            let endpoint = representative
                .values
                .row(representative.values.nrows() - 1)
                .iter()
                .copied()
                .collect::<Vec<_>>();
            let max_endpoint_delta = reference
                .iter()
                .zip(&endpoint)
                .map(|(expected, actual)| (expected - actual).abs())
                .fold(0.0, f64::max);
            assert!(
                max_endpoint_delta <= 2e-6,
                "{} {} changed IVP endpoint by {max_endpoint_delta:e}",
                family_label(family),
                policy.label()
            );
            let native_rows = representative.diagnostics_report_rows();
            let counter = |name: &str| {
                native_rows
                    .iter()
                    .find(|row| row.field == format!("native_counter::{name}"))
                    .and_then(|row| row.value.parse::<u64>().ok())
            };
            let parallel = counter("parallel_dispatches");
            let sequential = counter("sequential_dispatches");
            if family != ReactorIvpSolverFamily::Bdf {
                assert!(parallel.unwrap_or_default() + sequential.unwrap_or_default() > 0);
                if matches!(policy, CallbackPolicy::ForcedParallel) {
                    assert!(parallel.unwrap_or_default() > 0);
                }
                if matches!(policy, CallbackPolicy::Sequential) {
                    assert_eq!(parallel, Some(0));
                }
            }
            let callback_ms = native_rows
                .iter()
                .filter(|row| row.field.starts_with("native_timing_ms::callback_"))
                .filter_map(|row| row.value.parse::<f64>().ok())
                .sum::<f64>();
            rows.push(CallbackPolicyRow {
                solver: family_label(family).into(),
                policy: policy.label().into(),
                median_workflow_ms: format!("{workflow_ms:.4}"),
                callback_ms: if parallel.unwrap_or_default() + sequential.unwrap_or_default() == 0 {
                    "not reported".into()
                } else {
                    format!("{callback_ms:.4}")
                },
                parallel_dispatches: parallel
                    .map_or_else(|| "not reported".into(), |n| n.to_string()),
                sequential_dispatches: sequential
                    .map_or_else(|| "not reported".into(), |n| n.to_string()),
                max_endpoint_delta: format!("{max_endpoint_delta:.3e}"),
            });
        }
    }
    rows
}

fn workflow_elapsed_ms(snapshot: &IvpSolveSnapshot) -> f64 {
    snapshot
        .stage_telemetry
        .stages
        .iter()
        .find(|stage| stage.name == "total_workflow")
        .and_then(|stage| stage.elapsed_ms)
        .expect("detailed telemetry should report total workflow time")
}

fn bench_callback_policy_axis(c: &mut Criterion) {
    let mut group = c.benchmark_group("reactor_ivp_callback_policy/warm-inlet-900K");
    group.sample_size(20);
    for family in [
        ReactorIvpSolverFamily::Lsode2,
        ReactorIvpSolverFamily::Bdf,
        ReactorIvpSolverFamily::Radau,
    ] {
        let route = Route {
            family,
            frontend: ReactorIvpSymbolicBackend::ExprLegacy,
            matrix: ReactorIvpMatrixBackend::Dense,
        };
        for policy in [
            CallbackPolicy::Sequential,
            CallbackPolicy::Auto,
            CallbackPolicy::ForcedParallel,
        ] {
            group.bench_function(
                BenchmarkId::new(family_label(family), policy.label()),
                |b| {
                    b.iter(|| {
                        let (snapshot, _) = solve_fixture_once_with_policy(
                            route,
                            SolverTelemetryMode::Off,
                            Fixture::WarmInlet,
                            Some(policy),
                        );
                        black_box(snapshot.values[(snapshot.time_points() - 1, 2)])
                    });
                },
            );
        }
    }
    group.finish();
}

fn family_label(family: ReactorIvpSolverFamily) -> &'static str {
    match family {
        ReactorIvpSolverFamily::Lsode2 => "LSODE2",
        ReactorIvpSolverFamily::Bdf => "BDF",
        ReactorIvpSolverFamily::Radau => "Radau",
    }
}

fn report_diagnostics() {
    let mut macros = Vec::new();
    let mut stages = Vec::new();
    let mut native = Vec::new();

    for fixture in [Fixture::Baseline, Fixture::WarmInlet] {
        for route in ROUTES {
            let label = format!("{} / {}", fixture.label(), route.label());
            let mut workflow_values = Vec::with_capacity(DIAGNOSTIC_REPETITIONS);
            let mut wall_values = Vec::with_capacity(DIAGNOSTIC_REPETITIONS);
            let mut diagnostic_samples = Vec::with_capacity(DIAGNOSTIC_REPETITIONS);
            let mut latest = None;
            for _ in 0..DIAGNOSTIC_REPETITIONS {
                let (snapshot, wall_ms) =
                    solve_fixture_once(route, SolverTelemetryMode::Detailed, fixture);
                wall_values.push(wall_ms);
                let report_rows = snapshot.diagnostics_report_rows();
                let workflow_ms = snapshot
                    .stage_telemetry
                    .stages
                    .iter()
                    .find(|stage| stage.name == "total_workflow")
                    .and_then(|stage| stage.elapsed_ms)
                    .expect("detailed telemetry should report total workflow time");
                workflow_values.push(workflow_ms);
                diagnostic_samples.push((workflow_ms, report_rows));
                latest = Some(snapshot);
            }
            let snapshot = latest.expect("diagnostic run should produce a snapshot");
            let workflow_ms = median(&mut workflow_values);
            let (representative_workflow_ms, representative_rows) = diagnostic_samples
                .into_iter()
                .min_by(|left, right| {
                    (left.0 - workflow_ms)
                        .abs()
                        .total_cmp(&(right.0 - workflow_ms).abs())
                })
                .expect("diagnostic samples should be present");
            let solve_call_wall_ms = median(&mut wall_values);
            macros.push(MacroRow {
                route: label.clone(),
                status: snapshot.status().to_string(),
                states: snapshot.variable_count(),
                output_points: snapshot.time_points(),
                solve_call_wall_ms: format!("{solve_call_wall_ms:.3}"),
                workflow_ms: format!("{workflow_ms:.3}"),
            });

            for row in &representative_rows {
                if !row.field.starts_with("ivp_stage::") {
                    continue;
                }
                let mut parts = row.field.split("::");
                let _ = parts.next();
                let Some(stage_name) = parts.next() else {
                    continue;
                };
                let Some(metric) = parts.next() else { continue };
                if stage_name == "total_workflow" || metric != "elapsed_ms" {
                    continue;
                }
                let Ok(value) = row.value.parse::<f64>() else {
                    continue;
                };
                let share = 100.0 * value / representative_workflow_ms;
                if share < 1.0 {
                    continue;
                }
                let calls = representative_rows
                    .iter()
                    .find(|candidate| candidate.field == format!("ivp_stage::{stage_name}::calls"))
                    .map(|candidate| candidate.value.clone())
                    .unwrap_or_else(|| "n/a".to_string());
                stages.push(StageRow {
                    route: label.clone(),
                    stage: stage_name.to_string(),
                    elapsed_ms: format!("{value:.4}"),
                    workflow_share_pct: format!("{share:.1}%"),
                    calls,
                });
            }
            native.extend(representative_rows.iter().filter_map(|row| {
                let metric = row.field.as_str();
                let Ok(value) = row.value.parse::<f64>() else {
                    return None;
                };
                let keep_counter = matches!(
                    metric,
                    "native_counter::accepted_steps"
                        | "native_counter::rejected_steps"
                        | "native_counter::residual_calls"
                        | "native_counter::jacobian_calls"
                        | "native_counter::linear_solve_attempts"
                        | "native_counter::residual_evaluations"
                        | "native_counter::jacobian_evaluations"
                        | "native_counter::linear_solves"
                        | "native_counter::newton_iterations"
                        | "native_counter::factorizations"
                );
                (value > 0.0 && (metric.starts_with("native_timing_ms::") || keep_counter)).then(
                    || NativeRow {
                        route: label.clone(),
                        metric: metric.to_string(),
                        value: format!("{value:.4}"),
                    },
                )
            }));
        }
    }

    stages.sort_by(|left: &StageRow, right: &StageRow| {
        left.route
            .cmp(&right.route)
            .then(left.stage.cmp(&right.stage))
    });
    native.sort_by(|left: &NativeRow, right: &NativeRow| {
        left.route
            .cmp(&right.route)
            .then(left.metric.cmp(&right.metric))
    });

    println!(
        "IVP PERF-IVP-01 | bounded synthetic H2=>2H fixtures: 520 K / 900 K, four states | x=[0,0.02] | rtol=1e-6 | atol=1e-8 | diagnostic repetitions={DIAGNOSTIC_REPETITIONS}"
    );
    println!(
        "\nMacro: instrumented wall time is separate from Criterion telemetry-Off measurements\n{}",
        Table::new(macros).with(Style::markdown())
    );
    println!(
        "\nKiThe stage detail: stage shares use KiThe total_workflow as denominator\n{}",
        Table::new(stages).with(Style::markdown())
    );
    println!(
        "\nRST native telemetry: solver-specific timings/counters; not additive to KiThe stages\n{}",
        Table::new(native).with(Style::markdown())
    );
    let callback_rows = report_callback_policy_axis();
    println!(
        "\nCallback policy axis: warm-inlet-900K, ExprLegacy/Dense; BDF callback dispatch telemetry is unavailable in the current RST snapshot\n{}",
        Table::new(callback_rows).with(Style::markdown())
    );
}

fn report_prepared_lifecycles() {
    let mut rows = Vec::new();
    for family in [ReactorIvpSolverFamily::Lsode2, ReactorIvpSolverFamily::Bdf] {
        let route = Route {
            family,
            frontend: ReactorIvpSymbolicBackend::ExprLegacy,
            matrix: ReactorIvpMatrixBackend::Dense,
        };
        for (lifecycle, initial_states) in [
            ("same initial state", vec![None; 5]),
            (
                "prepared initial-composition sweep",
                [0.84, 0.85, 0.86, 0.87, 0.88]
                    .into_iter()
                    .map(Some)
                    .collect(),
            ),
        ] {
            let mut task = build_task(route, SolverTelemetryMode::Off);
            task.setup_condensed_ivp()
                .expect("prepared lifecycle setup should pass");
            let base_state = task
                .build_initial_state_vector()
                .expect("fixture initial state should validate");
            let preparation_started = Instant::now();
            let mut prepared = task
                .build_condensed_stiff_solver()
                .expect("StiffAPI preparation should pass");
            let preparation_ms = preparation_started.elapsed().as_secs_f64() * 1_000.0;
            let first_started = Instant::now();
            prepared
                .try_solve()
                .expect("first prepared lifecycle solve should pass");
            let first_solve_ms = first_started.elapsed().as_secs_f64() * 1_000.0;

            let mut repeated_ms = Vec::with_capacity(initial_states.len());
            for initial_h2 in initial_states {
                let mut y0 = base_state.clone();
                if let Some(h2) = initial_h2 {
                    y0[2] = h2;
                    y0[3] = 1.0 - h2;
                }
                let started = Instant::now();
                prepared
                    .try_solve_with_parameter_values(y0, DVector::zeros(0))
                    .expect("prepared repeated lifecycle solve should pass");
                repeated_ms.push(started.elapsed().as_secs_f64() * 1_000.0);
            }
            rows.push(LifecycleRow {
                route: route.label(),
                lifecycle: lifecycle.to_string(),
                preparation_ms: format!("{preparation_ms:.4}"),
                first_solve_ms: format!("{first_solve_ms:.4}"),
                repeated_solve_median_ms: format!("{:.4}", median(&mut repeated_ms)),
                repeats: 5,
            });
        }
    }
    println!(
        "\nPrepared RST lifecycle: preparation is one-time; first and repeat solve timings exclude KiThe facade validation/publication\n{}",
        Table::new(rows).with(Style::markdown())
    );
}

fn tcc_available() -> bool {
    std::process::Command::new(if cfg!(windows) { "where" } else { "which" })
        .arg("tcc")
        .output()
        .map(|output| output.status.success())
        .unwrap_or(false)
}

fn report_lsode2_execution_axis() {
    if !matches!(
        std::env::var("KITHE_RUN_IVP_AOT_BENCH").ok().as_deref(),
        Some("1") | Some("true") | Some("TRUE") | Some("yes") | Some("YES")
    ) || !tcc_available()
    {
        println!(
            "\nLSODE2 AOT axis skipped: set KITHE_RUN_IVP_AOT_BENCH=1 and install tcc; AOT is not supported by standalone BDF or Radau."
        );
        return;
    }

    let mut rows = Vec::new();
    let mut frontend_endpoints: Vec<(String, Vec<f64>)> = Vec::new();
    for fixture in [Fixture::Baseline, Fixture::WarmInlet] {
        for method in [
            ReactorIvpMethod::Auto,
            ReactorIvpMethod::Bdf,
            ReactorIvpMethod::Adams,
        ] {
            for frontend in [
                ReactorIvpSymbolicBackend::ExprLegacy,
                ReactorIvpSymbolicBackend::AtomView,
            ] {
                for execution in [
                    ReactorIvpExecutionBackend::Lambdify,
                    ReactorIvpExecutionBackend::Aot {
                        toolchain: Lsode2AotToolchain::CTcc,
                        profile: Lsode2AotProfile::Release,
                    },
                ] {
                    let mut setup_ms = Vec::new();
                    let mut problem_build_ms = Vec::new();
                    let mut preparation_ms = Vec::new();
                    let mut integration_ms = Vec::new();
                    let mut total_ms = Vec::new();
                    let output_dir =
                        tempfile::tempdir().expect("AOT output directory should exist");
                    let mut last_report = None;

                    for _ in 0..DIAGNOSTIC_REPETITIONS {
                        let route = Route {
                            family: ReactorIvpSolverFamily::Lsode2,
                            frontend,
                            matrix: ReactorIvpMatrixBackend::Dense,
                        };
                        let config = route
                            .config(SolverTelemetryMode::Detailed)
                            .with_lsode2_controller(method)
                            .with_execution_backend(execution);
                        let mut task =
                            build_task_for(route, SolverTelemetryMode::Detailed, fixture);
                        task.set_solver_backend_config(config);
                        let total_started = Instant::now();
                        let started = Instant::now();
                        task.setup_condensed_ivp()
                            .expect("LSODE2 setup should pass");
                        setup_ms.push(started.elapsed().as_secs_f64() * 1_000.0);

                        let started = Instant::now();
                        let mut problem = task
                            .build_lsode2_problem_config()
                            .expect("LSODE2 problem should build");
                        if matches!(execution, ReactorIvpExecutionBackend::Aot { .. }) {
                            problem = problem.with_backend(Lsode2BackendConfig::dense_aot_c_tcc(
                                output_dir.path(),
                            ));
                        }
                        problem = problem.with_telemetry(IvpTelemetry::detailed());
                        problem_build_ms.push(started.elapsed().as_secs_f64() * 1_000.0);

                        let started = Instant::now();
                        let mut solver =
                            Lsode2Solver::new(problem).expect("LSODE2 solver should prepare");
                        solver.prepare().expect("LSODE2 callbacks should prepare");
                        preparation_ms.push(started.elapsed().as_secs_f64() * 1_000.0);
                        let started = Instant::now();
                        solver
                            .solve_with_summary()
                            .expect("LSODE2 integration should pass");
                        let (_, values) = solver.get_result();
                        assert!(values.iter().all(|value| value.is_finite()));
                        let final_row = values.row(values.nrows() - 1);
                        let final_temperature = final_row[0] * 90.0 + 80.0;
                        assert!((0.0..=3000.0).contains(&final_temperature));
                        assert!(
                            final_row
                                .iter()
                                .skip(2)
                                .all(|fraction| (0.0..=1.0).contains(fraction))
                        );
                        let endpoint = final_row.iter().copied().collect::<Vec<_>>();
                        let key = format!("{}/{method:?}/{frontend:?}", fixture.label());
                        if matches!(execution, ReactorIvpExecutionBackend::Lambdify) {
                            frontend_endpoints.push((key, endpoint));
                        } else {
                            let reference = frontend_endpoints
                                .iter()
                                .find(|(candidate, _)| candidate == &key)
                                .map(|(_, endpoint)| endpoint)
                                .expect("AOT route must have a matching Lambdify result");
                            for (index, (actual, expected)) in
                                endpoint.iter().zip(reference).enumerate()
                            {
                                assert!(
                                    (actual - expected).abs() <= 1.0e-8 * (1.0 + expected.abs()),
                                    "AOT/Lambdify mismatch at state[{index}]: {actual:e} vs {expected:e}"
                                );
                            }
                        }
                        integration_ms.push(started.elapsed().as_secs_f64() * 1_000.0);
                        total_ms.push(total_started.elapsed().as_secs_f64() * 1_000.0);
                        last_report = Some(solver.telemetry_snapshot());
                    }

                    let label = format!("{method:?}");
                    rows.push(ExecutionRow {
                        fixture: fixture.label().to_string(),
                        solver_mode: label,
                        frontend: format!("{frontend:?}"),
                        execution: format!("{execution:?}"),
                        setup_ms: format!("{:.4}", median(&mut setup_ms)),
                        problem_build_ms: format!("{:.4}", median(&mut problem_build_ms)),
                        rst_preparation_ms: format!("{:.4}", median(&mut preparation_ms)),
                        integration_ms: format!("{:.4}", median(&mut integration_ms)),
                        total_ms: format!("{:.4}", median(&mut total_ms)),
                        residual_calls: last_report
                            .as_ref()
                            .map_or(0, |report| report.residual_evaluations),
                        jacobian_calls: last_report
                            .as_ref()
                            .map_or(0, |report| report.jacobian_evaluations),
                    });
                    if matches!(execution, ReactorIvpExecutionBackend::Aot { .. }) {
                        if let Some(report) = last_report {
                            println!(
                                "\nRST detailed telemetry {method:?}/{frontend:?}/AOT:\n{}",
                                report.pretty_report()
                            );
                        }
                    }
                }
            }
        }
    }
    println!(
        "\nLSODE2 execution axis: same bounded fixture, Dense matrix layout, preparation and integration timed separately; AOT rows require tcc and are LSODE2-only\n{}",
        Table::new(rows).with(Style::markdown())
    );
}

fn median(values: &mut [f64]) -> f64 {
    values.sort_by(f64::total_cmp);
    values[values.len() / 2]
}

fn reactor_ivp_baseline(c: &mut Criterion) {
    report_diagnostics();
    report_prepared_lifecycles();
    report_lsode2_execution_axis();
    let mut group = c.benchmark_group("reactor_ivp_h2_dissociation");
    group.sample_size(20);
    for fixture in [Fixture::Baseline, Fixture::WarmInlet] {
        for route in ROUTES {
            group.bench_function(BenchmarkId::new(fixture.label(), route.label()), |b| {
                b.iter(|| {
                    let (snapshot, _) =
                        solve_fixture_once(route, SolverTelemetryMode::Off, fixture);
                    black_box(snapshot.values[(snapshot.time_points() - 1, 2)])
                });
            });
        }
    }
    group.finish();
    bench_callback_policy_axis(c);
}

criterion_group!(benches, reactor_ivp_baseline);
criterion_main!(benches);
