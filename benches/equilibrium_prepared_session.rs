//! PERF-PREP-02: cold construction versus a reusable fixed-P,T session.
//!
//! This benchmark measures the first prepared point separately from a later
//! point on the same public-facade session. It is intentionally fixed-phase:
//! P,H and phase-control routes have different ownership and rollback rules.
//!
//! Run in an optimized profile with:
//!
//! ```text
//! cargo bench --bench equilibrium_prepared_session --features equilibrium-instrumentation -- --noplot
//! ```

use std::time::Instant;

use criterion::{BatchSize, BenchmarkId, Criterion, black_box, criterion_group, criterion_main};
use tabled::{Table, Tabled};

use KiThe::Thermodynamics::ChemEquilibrium::prelude::*;

mod support;

use support::{
    PRESSURE_PA, REFERENCE_PRESSURE_PA, RP1311_40, STANJAN_15, TP1907_GAS_17, tp1907_element_seed,
};

#[derive(Clone, Copy, Debug)]
enum Case {
    SmallMolecular,
    Stanjan15Molecular,
    Tp1907Elemental,
    Rp1311Molecular,
    SmallElemental,
    Stanjan15Elemental,
}

impl Case {
    fn label(self) -> &'static str {
        match self {
            Self::SmallMolecular => "small3_molecular_pt",
            Self::Stanjan15Molecular => "stanjan15_molecular_pt",
            Self::Tp1907Elemental => "tp1907_18_elemental_pt",
            Self::Rp1311Molecular => "rp1311_40_molecular_pt",
            Self::SmallElemental => "small3_elemental_pt",
            Self::Stanjan15Elemental => "stanjan15_elemental_pt",
        }
    }

    fn temperature(self) -> f64 {
        match self {
            Self::Tp1907Elemental => 700.0,
            _ => 2_500.0,
        }
    }

    fn reference_pressure(self) -> f64 {
        match self {
            Self::Tp1907Elemental => PRESSURE_PA,
            _ => REFERENCE_PRESSURE_PA,
        }
    }

    fn components(self) -> usize {
        match self {
            Self::SmallMolecular | Self::SmallElemental => 3,
            Self::Stanjan15Molecular | Self::Stanjan15Elemental => 15,
            Self::Tp1907Elemental => 18,
            Self::Rp1311Molecular => 40,
        }
    }
}

fn base_builder(case: Case) -> EquilibriumCalculatorBuilder {
    match case {
        Case::SmallMolecular | Case::SmallElemental => EquilibriumCalculator::builder()
            .single_ideal_gas(["H2", "O2", "H2O"])
            .expect("small phase"),
        Case::Stanjan15Molecular | Case::Stanjan15Elemental => EquilibriumCalculator::builder()
            .ideal_gas_phase("gas", STANJAN_15)
            .expect("STANJAN phase"),
        Case::Tp1907Elemental => EquilibriumCalculator::builder()
            .ideal_gas_phase("gas", TP1907_GAS_17)
            .expect("TP-1907 gas phase")
            .pure_condensed_phase("graphite", "C(gr)", PhysicalState::Solid)
            .expect("graphite phase"),
        Case::Rp1311Molecular => EquilibriumCalculator::builder()
            .single_ideal_gas(RP1311_40)
            .expect("RP-1311 phase"),
    }
}

fn builder_for(case: Case) -> EquilibriumCalculatorBuilder {
    let builder = match case {
        Case::SmallMolecular => base_builder(case).initial_moles([0.1, 0.05, 1.9]),
        Case::Stanjan15Molecular => base_builder(case).initial_phase_moles([
            (
                PhaseComponentId::new(PhaseId::new(Some("gas".to_owned())), "CH4"),
                4.0,
            ),
            (
                PhaseComponentId::new(PhaseId::new(Some("gas".to_owned())), "CO2"),
                2.0,
            ),
            (
                PhaseComponentId::new(PhaseId::new(Some("gas".to_owned())), "O2"),
                8.0,
            ),
            (
                PhaseComponentId::new(PhaseId::new(Some("gas".to_owned())), "N2"),
                37.6,
            ),
        ]),
        Case::Tp1907Elemental => base_builder(case)
            .element_inventory(
                ElementInventory::from_amounts([
                    ("Ar", 0.05364816),
                    ("C", 1.00182736),
                    ("H", 2.0),
                    ("N", 8.9462088),
                    ("O", 2.40365472),
                ])
                .unwrap(),
            )
            .element_numerical_seed(tp1907_element_seed()),
        Case::Rp1311Molecular => base_builder(case).initial_phase_moles([
            (PhaseComponentId::new(PhaseId::new(None), "Ar"), 0.1),
            (PhaseComponentId::new(PhaseId::new(None), "CO2"), 0.1),
            (PhaseComponentId::new(PhaseId::new(None), "H2"), 1.0),
            (PhaseComponentId::new(PhaseId::new(None), "H2O"), 0.1),
            (PhaseComponentId::new(PhaseId::new(None), "N2"), 3.0),
            (PhaseComponentId::new(PhaseId::new(None), "O2"), 2.0),
        ]),
        Case::SmallElemental => base_builder(case)
            .element_inventory(ElementInventory::from_amounts([("H", 4.0), ("O", 2.0)]).unwrap()),
        Case::Stanjan15Elemental => base_builder(case).element_inventory(
            ElementInventory::from_amounts([("C", 6.0), ("H", 16.0), ("O", 20.0), ("N", 75.2)])
                .unwrap(),
        ),
    };
    builder
        .prefer_libraries(["NASA_gas", "NASA_cond"])
        .offline_only()
        .pressure_pa(PRESSURE_PA)
        .reference_pressure_pa(case.reference_pressure())
        .production_cascade()
        .preset(EquilibriumCalculatorPreset::Timed)
        .instrumentation(EquilibriumInstrumentationMode::Enabled)
        .at_temperature(case.temperature())
}

fn prepare(case: Case) -> PreparedEquilibriumPtSession {
    builder_for(case)
        .prepare_pt_session()
        .unwrap_or_else(|error| panic!("{} preparation failed: {error:?}", case.label()))
}

fn solve(case: Case, session: &mut PreparedEquilibriumPtSession) -> EquilibriumCalculatorPoint {
    session
        .solve_at_temperature(case.temperature())
        .unwrap_or_else(|error| panic!("{} solve failed: {error:?}", case.label()))
}

#[derive(Tabled)]
struct MacroRow {
    case: &'static str,
    components: usize,
    preparation_ms: String,
    cold_operation_ms: String,
    hot_point_ms: String,
    saved_vs_cold: String,
    rst_graph: &'static str,
}

#[derive(Tabled)]
struct StageRow {
    case: &'static str,
    scope: &'static str,
    total_ms: String,
    thermochemistry_ms: String,
    numeric_problem_ms: String,
    symbolic_ms: String,
    rst_preparation_ms: String,
    nonlinear_ms: String,
    validation_ms: String,
}

#[derive(Tabled)]
struct WorkflowRow {
    case: &'static str,
    one_shot_loop_ms: String,
    prepared_session_ms: String,
    builtin_range_ms: String,
    session_vs_loop: String,
    range_vs_loop: String,
}

fn ms(duration: std::time::Duration) -> String {
    format!("{:.3}", duration.as_secs_f64() * 1_000.0)
}

fn stage_row(case: Case, scope: &'static str, report: EquilibriumTimingReport) -> StageRow {
    StageRow {
        case: case.label(),
        scope,
        total_ms: ms(report.total()),
        thermochemistry_ms: ms(report.thermochemistry_preparation()),
        numeric_problem_ms: ms(report.numerical_problem_preparation()),
        symbolic_ms: ms(report.symbolic_construction()),
        rst_preparation_ms: ms(report.rst_preparation()),
        nonlinear_ms: ms(report.nonlinear_solve()),
        validation_ms: ms(report.validation()),
    }
}

fn capture(case: Case) -> (MacroRow, Vec<StageRow>) {
    let started = Instant::now();
    let mut session = prepare(case);
    let preparation = session.preparation_timing();
    let cold = solve(case, &mut session);
    let cold_wall = started.elapsed();
    let hot_temperature = case.temperature()
        + if case.temperature() < 1_000.0 {
            0.1
        } else {
            100.0
        };
    let hot = session
        .solve_at_temperature(hot_temperature)
        .unwrap_or_else(|error| panic!("{} hot solve failed: {error:?}", case.label()));
    let hot_report = *hot.solution().timing_report();
    let cold_report = *cold.solution().timing_report();
    let cold_ms = cold_wall.as_secs_f64() * 1_000.0;
    let hot_ms = hot_report.total().as_secs_f64() * 1_000.0;
    let saved = (1.0 - hot_ms / cold_ms) * 100.0;
    (
        MacroRow {
            case: case.label(),
            components: case.components(),
            preparation_ms: ms(preparation.total()),
            cold_operation_ms: format!("{cold_ms:.3}"),
            hot_point_ms: format!("{hot_ms:.3}"),
            saved_vs_cold: format!("{saved:.1}%"),
            rst_graph: if session.rst_graph_prepared() {
                "yes"
            } else {
                "no"
            },
        },
        vec![
            stage_row(case, "prepared construction", preparation),
            stage_row(case, "first point solve", cold_report),
            stage_row(case, "hot point", hot_report),
        ],
    )
}

fn temperatures(case: Case) -> [f64; 2] {
    [
        case.temperature(),
        case.temperature()
            + if case.temperature() < 1_000.0 {
                0.1
            } else {
                100.0
            },
    ]
}

fn one_shot_loop_wall(case: Case) -> std::time::Duration {
    let started = Instant::now();
    for temperature in temperatures(case) {
        let outcome = builder_for(case)
            .at_temperature(temperature)
            .solve()
            .unwrap_or_else(|error| panic!("{} one-shot failed: {error:?}", case.label()));
        black_box(outcome);
    }
    started.elapsed()
}

fn prepared_session_wall(case: Case) -> std::time::Duration {
    let started = Instant::now();
    let mut session = prepare(case);
    for temperature in temperatures(case) {
        let point = session
            .solve_at_temperature(temperature)
            .unwrap_or_else(|error| panic!("{} session failed: {error:?}", case.label()));
        black_box(point);
    }
    started.elapsed()
}

fn builtin_range_wall(case: Case) -> std::time::Duration {
    let started = Instant::now();
    let range = builder_for(case)
        .over_temperature_range(temperatures(case).to_vec())
        .unwrap_or_else(|error| panic!("{} range setup failed: {error:?}", case.label()))
        .solve()
        .unwrap_or_else(|error| panic!("{} range failed: {error:?}", case.label()));
    black_box(range);
    started.elapsed()
}

fn workflow_row(case: Case) -> WorkflowRow {
    let one_shot = one_shot_loop_wall(case);
    let prepared = prepared_session_wall(case);
    let range = builtin_range_wall(case);
    let one_shot_ms = one_shot.as_secs_f64() * 1_000.0;
    let prepared_ms = prepared.as_secs_f64() * 1_000.0;
    let range_ms = range.as_secs_f64() * 1_000.0;
    WorkflowRow {
        case: case.label(),
        one_shot_loop_ms: format!("{one_shot_ms:.3}"),
        prepared_session_ms: format!("{prepared_ms:.3}"),
        builtin_range_ms: format!("{range_ms:.3}"),
        session_vs_loop: format!("{:.1}%", (1.0 - prepared_ms / one_shot_ms) * 100.0),
        range_vs_loop: format!("{:.1}%", (1.0 - range_ms / one_shot_ms) * 100.0),
    }
}

fn cases() -> [Case; 6] {
    [
        Case::SmallMolecular,
        Case::Stanjan15Molecular,
        Case::Tp1907Elemental,
        Case::Rp1311Molecular,
        Case::SmallElemental,
        Case::Stanjan15Elemental,
    ]
}

fn print_capture() {
    let mut macro_rows = Vec::new();
    let mut stage_rows = Vec::new();
    let mut workflow_rows = Vec::new();
    for case in cases() {
        let (macro_row, mut rows) = capture(case);
        macro_rows.push(macro_row);
        stage_rows.append(&mut rows);
        workflow_rows.push(workflow_row(case));
    }
    println!("\nPERF-PREP-02 prepared P,T session macro table");
    println!("{}", Table::new(macro_rows));
    println!("\nPERF-PREP-02 concrete stage table (inclusive stages; do not sum nested columns)");
    println!("{}", Table::new(stage_rows));
    println!("\nPERF-PREP-02 equal two-point workflow comparison");
    println!("{}", Table::new(workflow_rows));
}

fn criterion_benchmark(c: &mut Criterion) {
    print_capture();
    let mut group = c.benchmark_group("prepared_pt_session");
    group.sample_size(10);
    for case in cases() {
        group.bench_function(
            BenchmarkId::new("cold_prepare_and_solve", case.label()),
            |bencher| {
                bencher.iter(|| {
                    let mut session = prepare(case);
                    black_box(solve(case, &mut session))
                });
            },
        );
        group.bench_function(
            BenchmarkId::new("one_shot_two_point_loop", case.label()),
            |bencher| bencher.iter(|| black_box(one_shot_loop_wall(case))),
        );
        group.bench_function(
            BenchmarkId::new("prepared_two_point_session", case.label()),
            |bencher| bencher.iter(|| black_box(prepared_session_wall(case))),
        );
        group.bench_function(
            BenchmarkId::new("builtin_two_point_range", case.label()),
            |bencher| bencher.iter(|| black_box(builtin_range_wall(case))),
        );
        group.bench_function(
            BenchmarkId::new("hot_solve_only", case.label()),
            |bencher| {
                bencher.iter_batched(
                    || prepare(case),
                    |mut session| black_box(solve(case, &mut session)),
                    BatchSize::SmallInput,
                );
            },
        );
    }
    group.finish();
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
