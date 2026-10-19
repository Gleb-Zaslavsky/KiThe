//! PERF-LEGACY-01: isolated callback timing for retained legacy backends.
//!
//! Every case below selects exactly one `Legacy(LM|NR|TR)` backend through the
//! public equilibrium facade. The diagnostic tables are intentionally separate:
//! accepted callback work is not averaged together with failed or rejected
//! attempts, while the ledger keeps every started attempt visible.
//!
//! The accepted table also carries the PERF-NL-01 measurement: adjacent
//! residual/Jacobian callbacks whose input iterate is identical by floating
//! point bits. This is evidence for or against a shared numeric workspace;
//! it is not itself an optimization.
//!
//! Run in an optimized profile with:
//!
//! ```text
//! cargo bench --bench equilibrium_legacy_callbacks --no-default-features \
//!   --features equilibrium-instrumentation -- --noplot
//! ```

use std::time::{Duration, Instant};

use criterion::{BenchmarkId, Criterion, black_box, criterion_group, criterion_main};
use tabled::{Table, Tabled};

use KiThe::Thermodynamics::ChemEquilibrium::prelude::*;

mod support;

use support::{
    PRESSURE_PA, REFERENCE_PRESSURE_PA, REFERENCE_TEMPERATURE_K, RP1311_40, STANJAN_15,
    TP1907_GAS_17, component, tp1907_element_seed,
};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Universe {
    Small3,
    Stanjan15,
    Tp1907_18,
    Rp1311_40,
}

impl Universe {
    fn label(self) -> &'static str {
        match self {
            Self::Small3 => "small3",
            Self::Stanjan15 => "stanjan15",
            Self::Tp1907_18 => "tp1907_18",
            Self::Rp1311_40 => "rp1311_40",
        }
    }

    fn temperature(self) -> f64 {
        match self {
            Self::Tp1907_18 => 700.0,
            _ => REFERENCE_TEMPERATURE_K,
        }
    }

    fn reference_pressure(self) -> f64 {
        match self {
            Self::Tp1907_18 => PRESSURE_PA,
            _ => REFERENCE_PRESSURE_PA,
        }
    }

    fn bounds(self) -> TemperatureBounds {
        let (lower, upper) = match self {
            Self::Tp1907_18 => (699.0, 700.0),
            _ => (
                (self.temperature() - 400.0).max(300.0),
                self.temperature() + 400.0,
            ),
        };
        TemperatureBounds::new(lower, upper).expect("benchmark bounds must be valid")
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Task {
    PtMolecular,
    PtElementInventory,
    PhMolecular,
    PhElementInventory,
}

impl Task {
    fn label(self) -> &'static str {
        match self {
            Self::PtMolecular => "pt_molecular",
            Self::PtElementInventory => "pt_element_inventory",
            Self::PhMolecular => "ph_molecular",
            Self::PhElementInventory => "ph_element_inventory",
        }
    }

    fn is_elemental(self) -> bool {
        matches!(self, Self::PtElementInventory | Self::PhElementInventory)
    }

    fn is_ph(self) -> bool {
        matches!(self, Self::PhMolecular | Self::PhElementInventory)
    }
}

#[derive(Clone, Copy, Debug)]
struct Case {
    universe: Universe,
    task: Task,
    backend: LegacyEquilibriumSolver,
}

impl Case {
    fn label(self) -> String {
        format!(
            "{}_{}_legacy_{:?}",
            self.universe.label(),
            self.task.label(),
            self.backend
        )
    }

    fn route(self) -> &'static str {
        if !self.task.is_ph() {
            return "P,T";
        }
        if matches!(self.universe, Universe::Tp1907_18 | Universe::Rp1311_40) {
            "P,H nested"
        } else {
            "P,H monolithic"
        }
    }
}

fn cases() -> Vec<Case> {
    let tasks = [
        (Universe::Small3, Task::PtMolecular),
        (Universe::Small3, Task::PtElementInventory),
        (Universe::Stanjan15, Task::PtMolecular),
        (Universe::Stanjan15, Task::PtElementInventory),
        (Universe::Tp1907_18, Task::PtElementInventory),
        (Universe::Rp1311_40, Task::PtElementInventory),
        (Universe::Small3, Task::PhMolecular),
        (Universe::Small3, Task::PhElementInventory),
        (Universe::Stanjan15, Task::PhMolecular),
        (Universe::Stanjan15, Task::PhElementInventory),
        (Universe::Tp1907_18, Task::PhElementInventory),
        (Universe::Rp1311_40, Task::PhElementInventory),
    ];
    [
        LegacyEquilibriumSolver::LM,
        LegacyEquilibriumSolver::NR,
        LegacyEquilibriumSolver::TR,
    ]
    .into_iter()
    .flat_map(|backend| {
        tasks.into_iter().map(move |(universe, task)| Case {
            universe,
            task,
            backend,
        })
    })
    .collect()
}

fn inventory(universe: Universe) -> ElementInventory {
    let values = match universe {
        Universe::Small3 => vec![("H", 4.0), ("O", 2.0)],
        Universe::Stanjan15 => vec![("C", 6.0), ("H", 16.0), ("O", 20.0), ("N", 75.2)],
        Universe::Tp1907_18 => vec![
            ("Ar", 0.05364816),
            ("C", 1.00182736),
            ("H", 2.0),
            ("N", 8.9462088),
            ("O", 2.40365472),
        ],
        Universe::Rp1311_40 => vec![("Ar", 0.1), ("C", 0.1), ("H", 2.2), ("O", 4.3), ("N", 6.0)],
    };
    ElementInventory::from_amounts(values).expect("benchmark inventory must validate")
}

fn base_builder(universe: Universe) -> EquilibriumCalculatorBuilder {
    match universe {
        Universe::Small3 => EquilibriumCalculator::builder()
            .single_ideal_gas(["H2", "O2", "H2O"])
            .expect("small phase"),
        Universe::Stanjan15 => EquilibriumCalculator::builder()
            .ideal_gas_phase("gas", STANJAN_15)
            .expect("STANJAN phase"),
        Universe::Tp1907_18 => EquilibriumCalculator::builder()
            .ideal_gas_phase("gas", TP1907_GAS_17)
            .expect("TP-1907 gas phase")
            .pure_condensed_phase("graphite", "C(gr)", PhysicalState::Solid)
            .expect("graphite phase"),
        Universe::Rp1311_40 => EquilibriumCalculator::builder()
            .single_ideal_gas(RP1311_40)
            .expect("RP-1311 phase"),
    }
}

fn molecular_builder(universe: Universe) -> EquilibriumCalculatorBuilder {
    match universe {
        Universe::Small3 => base_builder(universe).initial_moles([0.1, 0.05, 1.9]),
        Universe::Stanjan15 => base_builder(universe).initial_phase_moles([
            (component(Some("gas"), "CH4"), 4.0),
            (component(Some("gas"), "CO2"), 2.0),
            (component(Some("gas"), "O2"), 8.0),
            (component(Some("gas"), "N2"), 37.6),
        ]),
        _ => unreachable!("molecular benchmark is defined only for gas controls"),
    }
}

fn builder_for(case: Case, target_enthalpy: f64) -> EquilibriumCalculatorBuilder {
    let builder = if case.task.is_elemental() {
        base_builder(case.universe).element_inventory(inventory(case.universe))
    } else {
        molecular_builder(case.universe)
    };
    let builder = if matches!(
        case,
        Case {
            universe: Universe::Tp1907_18,
            task: Task::PhElementInventory,
            ..
        }
    ) {
        builder.element_numerical_seed(tp1907_element_seed())
    } else {
        builder
    };
    let builder = builder
        .prefer_libraries(["NASA_gas", "NASA_cond"])
        .offline_only()
        .pressure_pa(PRESSURE_PA)
        .reference_pressure_pa(case.universe.reference_pressure())
        .solver_policy(SolverPolicy::Single(SolverBackend::Legacy(case.backend)))
        .expect("single legacy policy must validate")
        .preset(EquilibriumCalculatorPreset::Timed)
        .instrumentation(EquilibriumInstrumentationMode::Enabled);
    if case.task.is_ph() {
        let builder = if matches!(case.universe, Universe::Tp1907_18 | Universe::Rp1311_40) {
            builder.ph_solve_mode(PhSolveMode::NestedTemperature)
        } else {
            builder.ph_solve_mode(PhSolveMode::Monolithic)
        };
        builder.at_total_enthalpy(
            TotalEnthalpyJoules::new(target_enthalpy).expect("target must be finite"),
            case.universe.temperature(),
            case.universe.bounds(),
        )
    } else {
        builder.at_temperature(case.universe.temperature())
    }
}

fn reference_enthalpy(universe: Universe, elemental: bool) -> f64 {
    if matches!(universe, Universe::Tp1907_18) {
        return -420_039.2761;
    }
    let builder = if elemental {
        base_builder(universe).element_inventory(inventory(universe))
    } else {
        molecular_builder(universe)
    }
    .at_temperature(universe.temperature());
    let outcome = builder
        .prefer_libraries(["NASA_gas", "NASA_cond"])
        .offline_only()
        .pressure_pa(PRESSURE_PA)
        .reference_pressure_pa(universe.reference_pressure())
        .production_cascade()
        .preset(EquilibriumCalculatorPreset::Silent)
        .solve()
        .expect("reference P,T solve");
    let EquilibriumCalculatorOutcome::PtPoint(point) = outcome else {
        panic!("reference must be a P,T point")
    };
    let thermochemistry = ResolvedThermochemistry::from_resolved_system(point.resolved())
        .expect("reference thermochemistry");
    thermochemistry
        .enthalpy_model()
        .evaluate_total(point.solution().component_moles(), universe.temperature())
        .expect("reference enthalpy")
}

fn micros(duration: Duration) -> f64 {
    duration.as_secs_f64() * 1.0e6
}

fn milliseconds(duration: Duration) -> String {
    format!("{:.3}", duration.as_secs_f64() * 1_000.0)
}

fn percent(part: f64, total: f64) -> String {
    if total <= 0.0 {
        return "n/a".to_owned();
    }
    format!("{:.1}%", 100.0 * part / total)
}

fn pair_share(pairs: usize, residual_calls: usize, jacobian_calls: usize) -> String {
    let possible = residual_calls.min(jacobian_calls);
    percent(pairs as f64, possible as f64)
}

#[derive(Tabled)]
struct AcceptedCallbackRow {
    case: String,
    attempt: usize,
    route: &'static str,
    backend: String,
    backend_wall_ms: String,
    iterations: usize,
    residual_calls: usize,
    residual_ms: String,
    residual_share: String,
    jacobian_calls: usize,
    jacobian_ms: String,
    jacobian_share: String,
    same_iterate_pairs: String,
    same_pair_share: String,
    linear_solves: usize,
    linear_ms: String,
    linear_share: String,
    solver_overhead_ms: String,
    overhead_share: String,
}

#[derive(Tabled)]
struct AttemptLedgerRow {
    case: String,
    attempt: usize,
    route: &'static str,
    backend: String,
    outcome: String,
    wall_ms: String,
    metrics: &'static str,
    iterations: String,
    residual: String,
    jacobian: String,
    same_iterate_pairs: String,
    linear: String,
    overhead: String,
}

fn outcome_label(outcome: &SolverAttemptOutcome) -> String {
    match outcome {
        SolverAttemptOutcome::Accepted => "accepted".to_owned(),
        SolverAttemptOutcome::Failed { kind, .. } => format!("failed:{kind:?}"),
        SolverAttemptOutcome::RejectedCandidate { .. } => "rejected_candidate".to_owned(),
        SolverAttemptOutcome::Skipped { .. } => "skipped".to_owned(),
    }
}

fn attempt_wall(attempt: &SolverAttemptReport) -> Option<Duration> {
    if let Some(metrics) = attempt.metrics.as_ref() {
        return Some(Duration::from_micros(metrics.elapsed_micros as u64));
    }
    match attempt.outcome {
        SolverAttemptOutcome::Failed { elapsed_micros, .. } if elapsed_micros > 0 => {
            Some(Duration::from_micros(elapsed_micros as u64))
        }
        _ => None,
    }
}

fn add_attempt_rows(
    case: Case,
    attempts: &[SolverAttemptReport],
    accepted: &mut Vec<AcceptedCallbackRow>,
    ledger: &mut Vec<AttemptLedgerRow>,
) {
    for (attempt_index, attempt) in attempts.iter().enumerate() {
        let attempt_number = attempt_index + 1;
        let backend = format!("{:?}", attempt.backend);
        let wall = attempt_wall(attempt);
        if let (SolverAttemptOutcome::Accepted, Some(metrics), Some(timing), Some(wall)) = (
            &attempt.outcome,
            attempt.metrics.as_ref(),
            attempt
                .metrics
                .as_ref()
                .and_then(|metrics| metrics.evaluation_timing.as_ref()),
            wall,
        ) {
            let total_us = micros(wall);
            let residual_us = timing.residual_evaluation_micros as f64;
            let jacobian_us = timing.jacobian_evaluation_micros as f64;
            let linear_us = timing.linear_system_micros.unwrap_or_default() as f64;
            let overhead_us = timing.solver_overhead_micros as f64;
            accepted.push(AcceptedCallbackRow {
                case: case.label(),
                attempt: attempt_number,
                route: case.route(),
                backend: backend.clone(),
                backend_wall_ms: milliseconds(wall),
                iterations: metrics.iterations,
                residual_calls: metrics.residual_evaluations,
                residual_ms: format!("{:.3}", residual_us / 1_000.0),
                residual_share: percent(residual_us, total_us),
                jacobian_calls: metrics.jacobian_evaluations,
                jacobian_ms: format!("{:.3}", jacobian_us / 1_000.0),
                jacobian_share: percent(jacobian_us, total_us),
                same_iterate_pairs: timing
                    .legacy_same_iterate_callback_pairs
                    .map_or_else(|| "n/a".to_owned(), |pairs| pairs.to_string()),
                same_pair_share: timing.legacy_same_iterate_callback_pairs.map_or_else(
                    || "n/a".to_owned(),
                    |pairs| {
                        pair_share(
                            pairs,
                            metrics.residual_evaluations,
                            metrics.jacobian_evaluations,
                        )
                    },
                ),
                linear_solves: metrics.linear_solves,
                linear_ms: format!("{:.3}", linear_us / 1_000.0),
                linear_share: percent(linear_us, total_us),
                solver_overhead_ms: format!("{:.3}", overhead_us / 1_000.0),
                overhead_share: percent(overhead_us, total_us),
            });
        }

        let timing = attempt
            .metrics
            .as_ref()
            .and_then(|metrics| metrics.evaluation_timing.as_ref());
        ledger.push(AttemptLedgerRow {
            case: case.label(),
            attempt: attempt_number,
            route: case.route(),
            backend,
            outcome: outcome_label(&attempt.outcome),
            wall_ms: wall.map_or_else(|| "n/a".to_owned(), milliseconds),
            metrics: if timing.is_some() {
                "reported"
            } else {
                "wall-only"
            },
            iterations: attempt
                .metrics
                .as_ref()
                .map_or_else(|| "n/a".to_owned(), |m| m.iterations.to_string()),
            residual: timing.map_or_else(
                || "n/a".to_owned(),
                |t| format!("{:.3} ms", t.residual_evaluation_micros as f64 / 1_000.0),
            ),
            jacobian: timing.map_or_else(
                || "n/a".to_owned(),
                |t| format!("{:.3} ms", t.jacobian_evaluation_micros as f64 / 1_000.0),
            ),
            same_iterate_pairs: timing
                .and_then(|t| t.legacy_same_iterate_callback_pairs)
                .map_or_else(|| "n/a".to_owned(), |pairs| pairs.to_string()),
            linear: timing.and_then(|t| t.linear_system_micros).map_or_else(
                || "n/a".to_owned(),
                |us| format!("{:.3} ms", us as f64 / 1_000.0),
            ),
            overhead: timing.map_or_else(
                || "n/a".to_owned(),
                |t| format!("{:.3} ms", t.solver_overhead_micros as f64 / 1_000.0),
            ),
        });
    }
}

fn attempts_from_ph_report(report: &PhTemperatureSolveReport) -> Vec<SolverAttemptReport> {
    if let Some(evidence) = report.monolithic_evidence() {
        return evidence.solve_report().attempts.clone();
    }
    report
        .trials()
        .iter()
        .filter_map(|trial| trial.inner_evidence())
        .flat_map(|evidence| evidence.solve_report().attempts.iter().cloned())
        .collect()
}

fn attempts_from_error(error: &EquilibriumCalculatorError) -> Vec<SolverAttemptReport> {
    match error {
        EquilibriumCalculatorError::Preparation(cause)
        | EquilibriumCalculatorError::Solve(cause) => cause.backend_attempt_reports(),
        _ => Vec::new(),
    }
}

fn collect_case(
    case: Case,
    target: f64,
    accepted: &mut Vec<AcceptedCallbackRow>,
    ledger: &mut Vec<AttemptLedgerRow>,
) {
    let started = Instant::now();
    let result = builder_for(case, target).solve();
    let elapsed = started.elapsed();
    let attempts = match &result {
        Ok(EquilibriumCalculatorOutcome::PtPoint(point)) => {
            point.solution().solve_report().attempts.clone()
        }
        Ok(EquilibriumCalculatorOutcome::PhPoint(point)) => {
            attempts_from_ph_report(point.solution().report())
        }
        Ok(_) => Vec::new(),
        Err(error) => attempts_from_error(error),
    };
    add_attempt_rows(case, &attempts, accepted, ledger);
    let status = if result.is_ok() { "accepted" } else { "failed" };
    println!(
        "PERF-LEGACY-01 {} status={} wall_ms={} attempts={}",
        case.label(),
        status,
        milliseconds(elapsed),
        attempts.len()
    );
}

fn print_dashboard() {
    let mut accepted = Vec::new();
    let mut ledger = Vec::new();
    let target_cache = [
        (
            Universe::Small3,
            false,
            reference_enthalpy(Universe::Small3, false),
        ),
        (
            Universe::Small3,
            true,
            reference_enthalpy(Universe::Small3, true),
        ),
        (
            Universe::Stanjan15,
            false,
            reference_enthalpy(Universe::Stanjan15, false),
        ),
        (
            Universe::Stanjan15,
            true,
            reference_enthalpy(Universe::Stanjan15, true),
        ),
        (
            Universe::Tp1907_18,
            true,
            reference_enthalpy(Universe::Tp1907_18, true),
        ),
        (
            Universe::Rp1311_40,
            true,
            reference_enthalpy(Universe::Rp1311_40, true),
        ),
    ];
    for case in cases() {
        let elemental = case.task.is_elemental();
        let target = target_cache
            .iter()
            .find(|(universe, is_elemental, _)| {
                *universe == case.universe && *is_elemental == elemental
            })
            .expect("target cache must cover the case")
            .2;
        collect_case(case, target, &mut accepted, &mut ledger);
    }
    println!("\nPERF-LEGACY-01 accepted callback breakdown");
    println!(
        "Only accepted attempts are shown. Percentages use that attempt's backend wall time; this table is the callback optimization input."
    );
    println!("{}", Table::new(accepted));
    println!("\nPERF-LEGACY-01 complete attempt ledger");
    println!(
        "Every started attempt is retained, including rejected candidates and failed calls. `wall-only` means this backend has no structured callback split."
    );
    println!("{}", Table::new(ledger));
}

fn criterion_benchmark(c: &mut Criterion) {
    print_dashboard();
    let mut group = c.benchmark_group("equilibrium_legacy_callbacks");
    group.sample_size(10);
    group.measurement_time(Duration::from_secs(2));
    for case in cases() {
        let target = reference_enthalpy(case.universe, case.task.is_elemental());
        group.bench_function(BenchmarkId::from_parameter(case.label()), |bencher| {
            bencher.iter(|| black_box(builder_for(case, target).solve()));
        });
    }
    group.finish();
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
