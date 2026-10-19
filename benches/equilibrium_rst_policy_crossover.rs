//! PERF-RST-06: Sequential versus Parallel Lambdify crossover.
//!
//! This is a RustedSciThe-side microbenchmark driven from KiThe's workspace.
//! It keeps the symbolic fixture and parameter binding identical and changes
//! only `LambdifyExecutionPolicy`. The diagnostic table measures both prepared
//! payload construction and repeated residual/Jacobian callback work.
//!
//! Run an optimized probe with:
//!
//! ```text
//! cargo bench --bench equilibrium_rst_policy_crossover --no-default-features -- --noplot
//! ```

use std::hint::black_box;
use std::time::{Duration, Instant};

use criterion::{BenchmarkId, Criterion, criterion_group, criterion_main};
use nalgebra::{DMatrix, DVector};
use tabled::{Table, Tabled};

use RustedSciThe::numerical::Nonlinear_systems::prelude::{
    JacobianProvider, LambdifyExecutionPolicy, NonlinearProblem, PreparedSymbolicNonlinearProblem,
    SymbolicProblemOptions,
};
use RustedSciThe::symbolic::symbolic_engine::Expr;

const SIZES: [usize; 11] = [3, 15, 18, 40, 64, 80, 100, 128, 256, 512, 1024];
const PREPARATION_REPETITIONS: usize = 3;
const CALLBACK_REPETITIONS: usize = 40;

#[derive(Clone)]
struct Fixture {
    equations: Vec<Expr>,
    variables: Vec<String>,
}

#[derive(Tabled)]
struct CrossoverRow {
    equations: usize,
    #[tabled(rename = "prep seq ms")]
    preparation_sequential_ms: String,
    #[tabled(rename = "prep par ms")]
    preparation_parallel_ms: String,
    #[tabled(rename = "prep winner")]
    preparation_winner: String,
    #[tabled(rename = "callbacks seq us")]
    callbacks_sequential_us: String,
    #[tabled(rename = "callbacks par us")]
    callbacks_parallel_us: String,
    #[tabled(rename = "callback winner")]
    callback_winner: String,
    #[tabled(rename = "parallel / sequential")]
    callback_ratio: String,
}

fn fixture(size: usize) -> Fixture {
    let variables = (0..size)
        .map(|index| format!("x{index}"))
        .collect::<Vec<_>>();
    let equations = variables
        .iter()
        .map(|variable| {
            Expr::parse_expression(&format!(
                "({variable} * {variable}) + (0.01 * {variable}) + (0.001 * {variable} * {variable} * {variable}) - scale"
            ))
        })
        .collect::<Vec<_>>();
    Fixture {
        equations,
        variables,
    }
}

fn options(fixture: &Fixture, policy: LambdifyExecutionPolicy) -> SymbolicProblemOptions {
    SymbolicProblemOptions::new()
        .with_variables(fixture.variables.clone())
        .with_equation_parameters(vec!["scale".to_string()])
        .with_equation_parameter_values(DVector::from_element(1, 1.0))
        .with_lambdify_execution_policy(policy)
        .with_lambdify_backend()
}

fn prepare(fixture: &Fixture, policy: LambdifyExecutionPolicy) -> PreparedSymbolicNonlinearProblem {
    PreparedSymbolicNonlinearProblem::from_expressions(
        fixture.equations.clone(),
        options(fixture, policy),
    )
    .expect("crossover prepared problem must build")
}

fn elapsed_preparation_ms(fixture: &Fixture, policy: LambdifyExecutionPolicy) -> f64 {
    let started = Instant::now();
    black_box(prepare(fixture, policy));
    started.elapsed().as_secs_f64() * 1.0e3
}

fn callback_total_us(
    prepared: &PreparedSymbolicNonlinearProblem,
    size: usize,
    repetitions: usize,
) -> f64 {
    let bound = prepared
        .bind_values(DVector::from_element(1, 1.0))
        .expect("crossover parameter binding must succeed");
    let point = DVector::from_element(size, 0.9);
    let mut residual = DVector::zeros(size);
    let mut jacobian = DMatrix::zeros(size, size);
    let started = Instant::now();
    for _ in 0..repetitions {
        bound
            .residual_into(&point, &mut residual)
            .expect("crossover residual evaluation must succeed");
        bound
            .jacobian_into(&point, &mut jacobian)
            .expect("crossover Jacobian evaluation must succeed");
        black_box((&residual, &jacobian));
    }
    started.elapsed().as_secs_f64() * 1.0e6 / repetitions as f64
}

fn mean(values: &[f64]) -> f64 {
    values.iter().sum::<f64>() / values.len() as f64
}

fn format_ms(value: f64) -> String {
    format!("{value:.3}")
}

fn format_us(value: f64) -> String {
    format!("{value:.1}")
}

fn winner(left: f64, right: f64, left_name: &str, right_name: &str) -> String {
    if left <= right {
        left_name.to_string()
    } else {
        right_name.to_string()
    }
}

fn print_probe() {
    let mut rows = Vec::with_capacity(SIZES.len());
    for size in SIZES {
        let fixture = fixture(size);
        let preparation_sequential = mean(
            &(0..PREPARATION_REPETITIONS)
                .map(|_| elapsed_preparation_ms(&fixture, LambdifyExecutionPolicy::Sequential))
                .collect::<Vec<_>>(),
        );
        let preparation_parallel = mean(
            &(0..PREPARATION_REPETITIONS)
                .map(|_| {
                    elapsed_preparation_ms(
                        &fixture,
                        LambdifyExecutionPolicy::Parallel { min_work: 1 },
                    )
                })
                .collect::<Vec<_>>(),
        );

        let sequential = prepare(&fixture, LambdifyExecutionPolicy::Sequential);
        let parallel = prepare(&fixture, LambdifyExecutionPolicy::Parallel { min_work: 1 });
        let callbacks_sequential = mean(
            &(0..PREPARATION_REPETITIONS)
                .map(|_| callback_total_us(&sequential, size, CALLBACK_REPETITIONS))
                .collect::<Vec<_>>(),
        );
        let callbacks_parallel = mean(
            &(0..PREPARATION_REPETITIONS)
                .map(|_| callback_total_us(&parallel, size, CALLBACK_REPETITIONS))
                .collect::<Vec<_>>(),
        );

        rows.push(CrossoverRow {
            equations: size,
            preparation_sequential_ms: format_ms(preparation_sequential),
            preparation_parallel_ms: format_ms(preparation_parallel),
            preparation_winner: winner(
                preparation_sequential,
                preparation_parallel,
                "Sequential",
                "Parallel",
            ),
            callbacks_sequential_us: format_us(callbacks_sequential),
            callbacks_parallel_us: format_us(callbacks_parallel),
            callback_winner: winner(
                callbacks_sequential,
                callbacks_parallel,
                "Sequential",
                "Parallel",
            ),
            callback_ratio: format!("{:.3}x", callbacks_parallel / callbacks_sequential),
        });
    }

    println!("\nPERF-RST-06 Sequential versus Parallel crossover");
    println!("fixture=independent_parameterized_equations callback=repeated residual+Jacobian");
    println!("{}", Table::new(rows));
    println!(
        "Interpretation: callback ratio < 1.0 means Parallel is faster; preparation is reported separately."
    );
}

fn criterion_benchmark(c: &mut Criterion) {
    print_probe();
    if std::env::var_os("RST_CROSSOVER_PROBE_ONLY").is_some() {
        return;
    }
    let mut group = c.benchmark_group("equilibrium_rst_policy_crossover");
    group.sample_size(10);
    group.warm_up_time(Duration::from_secs(1));
    group.measurement_time(Duration::from_secs(2));

    for size in SIZES {
        let fixture = fixture(size);
        for (label, policy) in [
            ("sequential", LambdifyExecutionPolicy::Sequential),
            (
                "parallel",
                LambdifyExecutionPolicy::Parallel { min_work: 1 },
            ),
        ] {
            let prepared = prepare(&fixture, policy);
            group.bench_function(
                BenchmarkId::new("residual_plus_jacobian", format!("{size}_{label}")),
                |bencher| {
                    bencher.iter(|| {
                        black_box(callback_total_us(&prepared, size, 1));
                    });
                },
            );
        }
    }
    group.finish();
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
