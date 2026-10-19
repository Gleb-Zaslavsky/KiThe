//! PERF-LEGACY-02: explain the retained Legacy solver wall time.
//!
//! This is intentionally a focused diagnostic benchmark, not a second broad
//! equilibrium dashboard. Each case runs through the public facade with one
//! selected Legacy backend. The first table compares the complete facade call
//! with the accepted and all-attempt backend wall time. The second table keeps
//! every attempt, including failed attempts and P,H recovery attempts, and
//! expands the Legacy solver remainder into measurable internal operations.
//!
//! Run in release mode with:
//!
//! ```text
//! cargo bench --bench equilibrium_legacy_overhead --no-default-features \
//!   --features equilibrium-instrumentation -- --noplot
//! ```

use std::time::{Duration, Instant};

use criterion::{Criterion, black_box, criterion_group, criterion_main};
use tabled::{Table, Tabled};

use KiThe::Thermodynamics::ChemEquilibrium::equilibrium_solver_policy::SolverEvaluationTiming;
use KiThe::Thermodynamics::ChemEquilibrium::prelude::*;

mod support;

use support::{
    PRESSURE_PA, REFERENCE_PRESSURE_PA, REFERENCE_TEMPERATURE_K, RP1311_40, STANJAN_15,
    TP1907_GAS_17, component, tp1907_element_seed,
};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum CaseKind {
    SmallPtMolecular,
    StanjanPhMolecular,
    Tp1907PtElemental,
    Tp1907PhElemental,
    Rp1311PtElemental,
    Rp1311PhElemental,
}

impl CaseKind {
    fn label(self) -> &'static str {
        match self {
            Self::SmallPtMolecular => "small3_pt_molecular",
            Self::StanjanPhMolecular => "stanjan15_ph_molecular_monolithic",
            Self::Tp1907PtElemental => "tp1907_18_pt_element_inventory",
            Self::Tp1907PhElemental => "tp1907_18_ph_element_inventory_nested",
            Self::Rp1311PtElemental => "rp1311_40_pt_element_inventory",
            Self::Rp1311PhElemental => "rp1311_40_ph_element_inventory_nested",
        }
    }

    fn temperature(self) -> f64 {
        match self {
            Self::Tp1907PtElemental | Self::Tp1907PhElemental => 700.0,
            _ => REFERENCE_TEMPERATURE_K,
        }
    }

    fn pressure(self) -> f64 {
        PRESSURE_PA
    }

    fn reference_pressure(self) -> f64 {
        match self {
            Self::Tp1907PtElemental | Self::Tp1907PhElemental => PRESSURE_PA,
            _ => REFERENCE_PRESSURE_PA,
        }
    }

    fn is_ph(self) -> bool {
        matches!(
            self,
            Self::StanjanPhMolecular | Self::Tp1907PhElemental | Self::Rp1311PhElemental
        )
    }

    fn is_elemental(self) -> bool {
        !matches!(self, Self::SmallPtMolecular | Self::StanjanPhMolecular)
    }

    fn route(self) -> &'static str {
        match self {
            Self::StanjanPhMolecular => "P,H monolithic",
            Self::Tp1907PhElemental | Self::Rp1311PhElemental => "P,H nested",
            _ => "P,T",
        }
    }

    fn bounds(self) -> TemperatureBounds {
        let (lower, upper) = match self {
            Self::Tp1907PhElemental => (699.0, 700.0),
            _ => (
                (self.temperature() - 400.0).max(300.0),
                self.temperature() + 400.0,
            ),
        };
        TemperatureBounds::new(lower, upper).expect("benchmark bounds must be valid")
    }
}

#[derive(Clone, Copy, Debug)]
struct Case {
    kind: CaseKind,
    backend: LegacyEquilibriumSolver,
}

impl Case {
    fn label(self) -> String {
        format!("{}_legacy_{:?}", self.kind.label(), self.backend)
    }
}

fn cases() -> Vec<Case> {
    [
        CaseKind::SmallPtMolecular,
        CaseKind::StanjanPhMolecular,
        CaseKind::Tp1907PtElemental,
        CaseKind::Tp1907PhElemental,
        CaseKind::Rp1311PtElemental,
        CaseKind::Rp1311PhElemental,
    ]
    .into_iter()
    .flat_map(|kind| {
        [
            LegacyEquilibriumSolver::LM,
            LegacyEquilibriumSolver::NR,
            LegacyEquilibriumSolver::TR,
        ]
        .into_iter()
        .map(move |backend| Case { kind, backend })
    })
    .collect()
}

fn inventory(kind: CaseKind) -> ElementInventory {
    let values = match kind {
        CaseKind::Tp1907PtElemental | CaseKind::Tp1907PhElemental => vec![
            ("Ar", 0.05364816),
            ("C", 1.00182736),
            ("H", 2.0),
            ("N", 8.9462088),
            ("O", 2.40365472),
        ],
        CaseKind::Rp1311PtElemental | CaseKind::Rp1311PhElemental => {
            vec![("Ar", 0.1), ("C", 0.1), ("H", 2.2), ("O", 4.3), ("N", 6.0)]
        }
        _ => unreachable!("molecular cases do not use an elemental inventory"),
    };
    ElementInventory::from_amounts(values).expect("benchmark inventory must validate")
}

fn base_builder(kind: CaseKind) -> EquilibriumCalculatorBuilder {
    match kind {
        CaseKind::SmallPtMolecular => EquilibriumCalculator::builder()
            .single_ideal_gas(["H2", "O2", "H2O"])
            .expect("small phase"),
        CaseKind::StanjanPhMolecular => EquilibriumCalculator::builder()
            .ideal_gas_phase("gas", STANJAN_15)
            .expect("STANJAN phase"),
        CaseKind::Tp1907PtElemental | CaseKind::Tp1907PhElemental => {
            EquilibriumCalculator::builder()
                .ideal_gas_phase("gas", TP1907_GAS_17)
                .expect("TP-1907 gas phase")
                .pure_condensed_phase("graphite", "C(gr)", PhysicalState::Solid)
                .expect("graphite phase")
        }
        CaseKind::Rp1311PtElemental | CaseKind::Rp1311PhElemental => {
            EquilibriumCalculator::builder()
                .single_ideal_gas(RP1311_40)
                .expect("RP-1311 phase")
        }
    }
}

fn molecular_builder(kind: CaseKind) -> EquilibriumCalculatorBuilder {
    match kind {
        CaseKind::SmallPtMolecular => base_builder(kind).initial_moles([0.1, 0.05, 1.9]),
        CaseKind::StanjanPhMolecular => base_builder(kind).initial_phase_moles([
            (component(Some("gas"), "CH4"), 4.0),
            (component(Some("gas"), "CO2"), 2.0),
            (component(Some("gas"), "O2"), 8.0),
            (component(Some("gas"), "N2"), 37.6),
        ]),
        _ => unreachable!("molecular builder is defined only for molecular cases"),
    }
}

fn reference_enthalpy(kind: CaseKind) -> f64 {
    let builder = if kind.is_elemental() {
        base_builder(kind).element_inventory(inventory(kind))
    } else {
        molecular_builder(kind)
    };
    let outcome = builder
        .prefer_libraries(["NASA_gas", "NASA_cond"])
        .offline_only()
        .pressure_pa(kind.pressure())
        .reference_pressure_pa(kind.reference_pressure())
        .production_cascade()
        .preset(EquilibriumCalculatorPreset::Silent)
        .at_temperature(kind.temperature())
        .solve()
        .expect("reference P,T solve");
    let EquilibriumCalculatorOutcome::PtPoint(point) = outcome else {
        panic!("reference must be a P,T point")
    };
    let thermochemistry = ResolvedThermochemistry::from_resolved_system(point.resolved())
        .expect("reference thermochemistry");
    thermochemistry
        .enthalpy_model()
        .evaluate_total(point.solution().component_moles(), kind.temperature())
        .expect("reference enthalpy")
}

fn builder_for(case: Case, target_enthalpy: f64) -> EquilibriumCalculatorBuilder {
    let builder = if case.kind.is_elemental() {
        base_builder(case.kind).element_inventory(inventory(case.kind))
    } else {
        molecular_builder(case.kind)
    };
    let builder = if matches!(case.kind, CaseKind::Tp1907PhElemental) {
        builder.element_numerical_seed(tp1907_element_seed())
    } else {
        builder
    };
    let builder = builder
        .prefer_libraries(["NASA_gas", "NASA_cond"])
        .offline_only()
        .pressure_pa(case.kind.pressure())
        .reference_pressure_pa(case.kind.reference_pressure())
        .solver_policy(SolverPolicy::Single(SolverBackend::Legacy(case.backend)))
        .expect("single legacy policy must validate")
        .preset(EquilibriumCalculatorPreset::Timed)
        .instrumentation(EquilibriumInstrumentationMode::Enabled);
    if case.kind.is_ph() {
        let builder = if matches!(
            case.kind,
            CaseKind::Tp1907PhElemental | CaseKind::Rp1311PhElemental
        ) {
            builder.ph_solve_mode(PhSolveMode::NestedTemperature)
        } else {
            builder.ph_solve_mode(PhSolveMode::Monolithic)
        };
        builder.at_total_enthalpy(
            TotalEnthalpyJoules::new(target_enthalpy).expect("target must be finite"),
            case.kind.temperature(),
            case.kind.bounds(),
        )
    } else {
        builder.at_temperature(case.kind.temperature())
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

fn ms(duration: Duration) -> String {
    format!("{:.3}", duration.as_secs_f64() * 1_000.0)
}

fn share(part: u128, total: u128) -> String {
    if total == 0 {
        return "n/a".to_owned();
    }
    format!("{:.1}%", 100.0 * part as f64 / total as f64)
}

fn stage_ms(micros: u128) -> String {
    format!("{:.3}", micros as f64 / 1_000.0)
}

fn outcome_label(outcome: &SolverAttemptOutcome) -> &'static str {
    match outcome {
        SolverAttemptOutcome::Accepted => "accepted",
        SolverAttemptOutcome::Failed { .. } => "failed",
        SolverAttemptOutcome::RejectedCandidate { .. } => "rejected",
        SolverAttemptOutcome::Skipped { .. } => "skipped",
    }
}

#[derive(Tabled)]
struct SummaryRow {
    case: String,
    backend: String,
    route: &'static str,
    facade_ms: String,
    accepted_ms: String,
    accepted_facade_share: String,
    all_attempts_ms: String,
    all_attempts_facade_share: String,
    accepted_overhead: String,
    accepted_detail_ms: String,
    accepted_unattributed_ms: String,
}

#[derive(Tabled)]
struct DetailRow {
    case: String,
    attempt: usize,
    backend: String,
    outcome: &'static str,
    wall_ms: String,
    callbacks_ms: String,
    linear_ms: String,
    matrix_ms: String,
    norms_ms: String,
    trial_ms: String,
    feasible_ms: String,
    global_ms: String,
    unattributed_ms: String,
}

fn timing_parts(timing: &SolverEvaluationTiming) -> [u128; 8] {
    [
        timing.residual_evaluation_micros,
        timing.jacobian_evaluation_micros,
        timing.linear_system_micros.unwrap_or_default(),
        timing.matrix_operations_micros,
        timing.norm_micros,
        timing.trial_vector_micros,
        timing.feasibility_micros,
        timing.globalization_micros,
    ]
}

fn detail_total(parts: [u128; 8]) -> u128 {
    parts.into_iter().sum()
}

fn add_attempt_rows(
    case: Case,
    facade_wall: Duration,
    attempts: &[SolverAttemptReport],
    summaries: &mut Vec<SummaryRow>,
    details: &mut Vec<DetailRow>,
) {
    let mut accepted_row = None;
    let mut all_attempts_us = 0_u128;
    for (index, attempt) in attempts.iter().enumerate() {
        let Some(metrics) = attempt.metrics.as_ref() else {
            continue;
        };
        let Some(timing) = metrics.evaluation_timing.as_ref() else {
            continue;
        };
        let wall_us = metrics.elapsed_micros;
        all_attempts_us = all_attempts_us.saturating_add(wall_us);
        let parts = timing_parts(timing);
        let detail_us = detail_total(parts);
        let unattributed_us = timing.solver_overhead_micros.saturating_sub(
            detail_us
                .saturating_sub(parts[0])
                .saturating_sub(parts[1])
                .saturating_sub(parts[2]),
        );
        details.push(DetailRow {
            case: case.kind.label().to_owned(),
            attempt: index + 1,
            backend: format!("{:?}", case.backend),
            outcome: outcome_label(&attempt.outcome),
            wall_ms: format!(
                "{} ({})",
                stage_ms(wall_us),
                share(wall_us, facade_wall.as_micros())
            ),
            callbacks_ms: format!(
                "{} ({})",
                stage_ms(parts[0].saturating_add(parts[1])),
                share(parts[0].saturating_add(parts[1]), wall_us)
            ),
            linear_ms: format!("{} ({})", stage_ms(parts[2]), share(parts[2], wall_us)),
            matrix_ms: format!("{} ({})", stage_ms(parts[3]), share(parts[3], wall_us)),
            norms_ms: format!("{} ({})", stage_ms(parts[4]), share(parts[4], wall_us)),
            trial_ms: format!("{} ({})", stage_ms(parts[5]), share(parts[5], wall_us)),
            feasible_ms: format!("{} ({})", stage_ms(parts[6]), share(parts[6], wall_us)),
            global_ms: format!("{} ({})", stage_ms(parts[7]), share(parts[7], wall_us)),
            unattributed_ms: format!(
                "{} ({})",
                stage_ms(unattributed_us),
                share(unattributed_us, wall_us)
            ),
        });
        if matches!(attempt.outcome, SolverAttemptOutcome::Accepted) && accepted_row.is_none() {
            accepted_row = Some((wall_us, timing, unattributed_us));
        }
    }

    let zero_timing = SolverEvaluationTiming {
        residual_evaluation_micros: 0,
        jacobian_evaluation_micros: 0,
        legacy_same_iterate_callback_pairs: Some(0),
        linear_system_micros: Some(0),
        solver_overhead_micros: 0,
        matrix_operations_micros: 0,
        norm_micros: 0,
        trial_vector_micros: 0,
        feasibility_micros: 0,
        globalization_micros: 0,
        state_residual_evaluations: None,
        trial_residual_evaluations: None,
        state_jacobian_evaluations: None,
        trial_jacobian_evaluations: None,
        linear_factorization_micros: None,
        linear_system_solve_micros: None,
        linear_factorizations: None,
        accepted_steps: None,
        rejected_steps: None,
        reusable_trial_points: None,
        iteration_telemetry: None,
    };
    let (accepted_us, accepted_timing, accepted_unattributed_us) =
        accepted_row.unwrap_or((0, &zero_timing, 0));
    let accepted_detail_us = timing_parts(accepted_timing).into_iter().sum::<u128>();
    summaries.push(SummaryRow {
        case: case.kind.label().to_owned(),
        backend: format!("{:?}", case.backend),
        route: case.kind.route(),
        facade_ms: ms(facade_wall),
        accepted_ms: stage_ms(accepted_us),
        accepted_facade_share: if accepted_us == 0 {
            "n/a".to_owned()
        } else {
            share(accepted_us, facade_wall.as_micros())
        },
        all_attempts_ms: stage_ms(all_attempts_us),
        all_attempts_facade_share: share(all_attempts_us, facade_wall.as_micros()),
        accepted_overhead: share(accepted_timing.solver_overhead_micros, accepted_us),
        accepted_detail_ms: stage_ms(accepted_detail_us),
        accepted_unattributed_ms: stage_ms(accepted_unattributed_us),
    });
}

fn collect(
    case: Case,
    target_enthalpy: f64,
    summaries: &mut Vec<SummaryRow>,
    details: &mut Vec<DetailRow>,
) {
    let started = Instant::now();
    let result = builder_for(case, target_enthalpy).solve();
    let facade_wall = started.elapsed();
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
    add_attempt_rows(case, facade_wall, &attempts, summaries, details);
}

fn print_dashboard() {
    let mut summaries = Vec::new();
    let mut details = Vec::new();
    for case in cases() {
        collect(
            case,
            reference_enthalpy(case.kind),
            &mut summaries,
            &mut details,
        );
    }
    println!("PERF-LEGACY-04 macro facade/backend comparison");
    println!("All percentages in the detail table use the individual attempt wall time.");
    println!("{}", Table::new(summaries));
    println!("\nPERF-LEGACY-04 internal Legacy attempt detail");
    println!(
        "`unattributed` is the remaining Legacy wall after callbacks, linear systems, and the listed coarse operations."
    );
    println!("{}", Table::new(details));
}

fn criterion_benchmark(c: &mut Criterion) {
    print_dashboard();
    let mut group = c.benchmark_group("equilibrium_legacy_overhead");
    group.sample_size(10);
    group.measurement_time(Duration::from_secs(2));
    for case in cases() {
        let target = reference_enthalpy(case.kind);
        group.bench_function(case.label(), |bencher| {
            bencher.iter(|| black_box(builder_for(case, target).solve()));
        });
    }
    group.finish();
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
