//! PERF-PH-14: matched benchmark for native Gibbs coefficient-interval reuse.
//!
//! It contains two matched measurements: an isolated Gibbs-vector sequence
//! and a full changing-temperature P,H facade range. The first isolates the
//! closure-construction effect; the second checks whether that effect survives
//! alongside orchestration, nonlinear solve, and validation. Both use the
//! same public facade setup and a control feature that forces refreshes.
//!
//! Current path:
//!
//! ```text
//! cargo bench --bench equilibrium_gibbs_interval \
//!   --features "equilibrium-instrumentation" -- --noplot
//! ```
//!
//! Matched refresh-every-time control:
//!
//! ```text
//! cargo bench --bench equilibrium_gibbs_interval \
//!   --features "equilibrium-instrumentation equilibrium-benchmark-no-gibbs-interval-reuse" -- --noplot
//! ```

use std::time::Instant;

use criterion::{Criterion, black_box, criterion_group, criterion_main};
use tabled::{Table, Tabled};

use KiThe::Thermodynamics::ChemEquilibrium::prelude::*;

mod support;

use support::{PRESSURE_PA, STANJAN_15, TP1907_GAS_17, component};

const TEMPERATURES_K: [f64; 3] = [700.0, 700.5, 699.5];
const REPETITIONS: usize = 128;
const PH_TARGET_OFFSETS_J: [f64; 3] = [-100.0, 0.0, 100.0];

#[derive(Tabled)]
struct ProbeRow {
    mode: &'static str,
    vectors: usize,
    wall_ms: String,
    per_vector_us: String,
    gibbs_refreshes: usize,
    interval_reuses: usize,
}

#[derive(Tabled)]
struct PhProbeRow {
    case: &'static str,
    mode: &'static str,
    points: usize,
    wall_ms: String,
    mean_point_ms: String,
    gibbs_refreshes: usize,
    interval_reuses: usize,
    continuation_points: usize,
}

fn mode() -> &'static str {
    if cfg!(feature = "equilibrium-benchmark-no-gibbs-interval-reuse") {
        "refresh-every-temperature"
    } else {
        "native-interval-reuse"
    }
}

fn resolved_thermochemistry() -> ResolvedThermochemistry {
    let outcome = EquilibriumCalculator::builder()
        .ideal_gas_phase("gas", TP1907_GAS_17)
        .expect("TP-1907 gas phase")
        .pure_condensed_phase("graphite", "C(gr)", PhysicalState::Solid)
        .expect("TP-1907 graphite phase")
        // Use the frozen reference's gas-only molecular feed. The benchmark
        // measures the resolved NASA closure below, so it intentionally keeps
        // basin discovery out of the timed operation and out of its setup.
        .initial_phase_moles([
            (component(Some("gas"), "CH4"), 0.5),
            (component(Some("gas"), "CO2"), 0.50182736),
            (component(Some("gas"), "O2"), 0.2),
            (component(Some("gas"), "N2"), 4.4731044),
            (component(Some("gas"), "Ar"), 0.05364816),
        ])
        .prefer_libraries(["NASA_gas", "NASA_cond"])
        .offline_only()
        .pressure_pa(PRESSURE_PA)
        .reference_pressure_pa(100_000.0)
        .production_cascade()
        .preset(EquilibriumCalculatorPreset::Timed)
        .instrumentation(EquilibriumInstrumentationMode::Enabled)
        .at_temperature(700.0)
        .solve()
        .expect("TP-1907 P,T preparation");
    let EquilibriumCalculatorOutcome::PtPoint(point) = outcome else {
        panic!("TP-1907 benchmark must produce a P,T point")
    };
    ResolvedThermochemistry::from_resolved_system(point.resolved())
        .expect("TP-1907 thermochemistry bundle")
}

fn reference_enthalpy() -> f64 {
    let outcome = EquilibriumCalculator::builder()
        .ideal_gas_phase("gas", TP1907_GAS_17)
        .expect("TP-1907 gas phase")
        .pure_condensed_phase("graphite", "C(gr)", PhysicalState::Solid)
        .expect("TP-1907 graphite phase")
        .initial_phase_moles([
            (component(Some("gas"), "CH4"), 0.5),
            (component(Some("gas"), "CO2"), 0.50182736),
            (component(Some("gas"), "O2"), 0.2),
            (component(Some("gas"), "N2"), 4.4731044),
            (component(Some("gas"), "Ar"), 0.05364816),
        ])
        .prefer_libraries(["NASA_gas", "NASA_cond"])
        .offline_only()
        .pressure_pa(PRESSURE_PA)
        .reference_pressure_pa(100_000.0)
        .production_cascade()
        .preset(EquilibriumCalculatorPreset::Timed)
        .instrumentation(EquilibriumInstrumentationMode::Enabled)
        .at_temperature(700.0)
        .solve()
        .expect("TP-1907 P,T reference");
    let EquilibriumCalculatorOutcome::PtPoint(point) = outcome else {
        panic!("TP-1907 reference must produce a P,T point")
    };
    let thermochemistry = ResolvedThermochemistry::from_resolved_system(point.resolved())
        .expect("TP-1907 reference thermochemistry");
    thermochemistry
        .enthalpy_model()
        .evaluate_total(point.solution().component_moles(), 700.0)
        .expect("TP-1907 reference enthalpy")
}

fn stanjan15_reference_enthalpy() -> f64 {
    let outcome = EquilibriumCalculator::builder()
        .ideal_gas_phase("gas", STANJAN_15)
        .expect("STANJAN-15 gas phase")
        .initial_phase_moles([
            (component(Some("gas"), "CH4"), 4.0),
            (component(Some("gas"), "CO2"), 2.0),
            (component(Some("gas"), "O2"), 8.0),
            (component(Some("gas"), "N2"), 37.6),
        ])
        .prefer_libraries(["NASA_gas", "NASA_cond"])
        .offline_only()
        .pressure_pa(PRESSURE_PA)
        .reference_pressure_pa(100_000.0)
        .production_cascade()
        .preset(EquilibriumCalculatorPreset::Timed)
        .instrumentation(EquilibriumInstrumentationMode::Enabled)
        .at_temperature(2_500.0)
        .solve()
        .expect("STANJAN-15 P,T reference");
    let EquilibriumCalculatorOutcome::PtPoint(point) = outcome else {
        panic!("STANJAN-15 reference must produce a P,T point")
    };
    let thermochemistry = ResolvedThermochemistry::from_resolved_system(point.resolved())
        .expect("STANJAN-15 reference thermochemistry");
    thermochemistry
        .enthalpy_model()
        .evaluate_total(point.solution().component_moles(), 2_500.0)
        .expect("STANJAN-15 reference enthalpy")
}

fn solve_changing_temperature_ph(target: f64) -> EquilibriumCalculatorPhRange {
    let targets = PH_TARGET_OFFSETS_J
        .into_iter()
        .map(|offset| target + offset)
        .collect::<Vec<_>>();
    let outcome = EquilibriumCalculator::builder()
        .ideal_gas_phase("gas", TP1907_GAS_17)
        .expect("TP-1907 gas phase")
        .pure_condensed_phase("graphite", "C(gr)", PhysicalState::Solid)
        .expect("TP-1907 graphite phase")
        .initial_phase_moles([
            (component(Some("gas"), "CH4"), 0.5),
            (component(Some("gas"), "CO2"), 0.50182736),
            (component(Some("gas"), "O2"), 0.2),
            (component(Some("gas"), "N2"), 4.4731044),
            (component(Some("gas"), "Ar"), 0.05364816),
        ])
        .prefer_libraries(["NASA_gas", "NASA_cond"])
        .offline_only()
        .pressure_pa(PRESSURE_PA)
        .reference_pressure_pa(100_000.0)
        .production_cascade()
        .preset(EquilibriumCalculatorPreset::Timed)
        .instrumentation(EquilibriumInstrumentationMode::Enabled)
        .ph_solve_mode(PhSolveMode::NestedTemperature)
        .over_enthalpy_range(
            targets,
            700.0,
            TemperatureBounds::new(698.0, 702.0).expect("TP-1907 P,H bounds"),
        )
        .expect("TP-1907 P,H range request")
        .solve()
        .expect("TP-1907 changing-temperature P,H range");
    let EquilibriumCalculatorOutcome::PhRange(range) = outcome else {
        panic!("TP-1907 benchmark must produce a P,H range")
    };
    range
}

fn solve_changing_temperature_stanjan15_ph(target: f64) -> EquilibriumCalculatorPhRange {
    let targets = [-1_000.0, 0.0, 1_000.0]
        .into_iter()
        .map(|offset| target + offset)
        .collect::<Vec<_>>();
    let outcome = EquilibriumCalculator::builder()
        .ideal_gas_phase("gas", STANJAN_15)
        .expect("STANJAN-15 gas phase")
        .initial_phase_moles([
            (component(Some("gas"), "CH4"), 4.0),
            (component(Some("gas"), "CO2"), 2.0),
            (component(Some("gas"), "O2"), 8.0),
            (component(Some("gas"), "N2"), 37.6),
        ])
        .prefer_libraries(["NASA_gas", "NASA_cond"])
        .offline_only()
        .pressure_pa(PRESSURE_PA)
        .reference_pressure_pa(100_000.0)
        .production_cascade()
        .preset(EquilibriumCalculatorPreset::Timed)
        .instrumentation(EquilibriumInstrumentationMode::Enabled)
        .ph_solve_mode(PhSolveMode::Monolithic)
        .over_enthalpy_range(
            targets,
            2_500.0,
            TemperatureBounds::new(2_100.0, 2_900.0).expect("STANJAN-15 P,H bounds"),
        )
        .expect("STANJAN-15 P,H range request")
        .solve()
        .expect("STANJAN-15 changing-temperature P,H range");
    let EquilibriumCalculatorOutcome::PhRange(range) = outcome else {
        panic!("STANJAN-15 benchmark must produce a P,H range")
    };
    range
}

fn ph_probe_row(
    case: &'static str,
    range: EquilibriumCalculatorPhRange,
    started: Instant,
) -> PhProbeRow {
    let wall = started.elapsed();
    let report = range.solution().report();
    let counters = range.instrumentation_report().counters();
    PhProbeRow {
        case,
        mode: mode(),
        points: report.point_count(),
        wall_ms: format!("{:.3}", wall.as_secs_f64() * 1_000.0),
        mean_point_ms: format!(
            "{:.3}",
            wall.as_secs_f64() * 1_000.0 / report.point_count() as f64
        ),
        gibbs_refreshes: counters.gibbs_coefficient_refreshes(),
        interval_reuses: counters.gibbs_interval_reuses(),
        continuation_points: report.continuation_points(),
    }
}

fn evaluate_sequence(thermochemistry: &ResolvedThermochemistry) -> f64 {
    let mut checksum = 0.0;
    for _ in 0..REPETITIONS {
        for temperature in TEMPERATURES_K {
            let values = thermochemistry
                .evaluate_gibbs(temperature)
                .expect("Gibbs vector must remain evaluable");
            checksum += values.iter().copied().sum::<f64>();
            black_box(&values);
        }
    }
    black_box(checksum)
}

fn probe_row(thermochemistry: &ResolvedThermochemistry) -> ProbeRow {
    let started = Instant::now();
    let (_, report) = capture_equilibrium_instrumentation(|| evaluate_sequence(thermochemistry));
    let wall = started.elapsed();
    let vectors = REPETITIONS * TEMPERATURES_K.len();
    let per_vector_us = wall.as_secs_f64() * 1_000_000.0 / vectors as f64;
    let counters = report.counters();
    ProbeRow {
        mode: mode(),
        vectors,
        wall_ms: format!("{:.3}", wall.as_secs_f64() * 1_000.0),
        per_vector_us: format!("{per_vector_us:.2}"),
        gibbs_refreshes: counters.gibbs_coefficient_refreshes(),
        interval_reuses: counters.gibbs_interval_reuses(),
    }
}

fn criterion_benchmark(c: &mut Criterion) {
    let thermochemistry = resolved_thermochemistry();
    let tp1907_target_enthalpy = reference_enthalpy();
    let stanjan15_target_enthalpy = stanjan15_reference_enthalpy();
    println!("\nPERF-PH-14 matched Gibbs interval benchmark");
    println!("temperature_sequence={TEMPERATURES_K:?} repetitions={REPETITIONS}");
    println!("{}", Table::new([probe_row(&thermochemistry)]));
    println!("\nPERF-PH-14 changing-temperature P,H facade benchmark");
    println!("TP-1907 nested target_offsets_j={PH_TARGET_OFFSETS_J:?} initial_temperature_k=700.0");
    let tp1907_started = Instant::now();
    let tp1907_range = solve_changing_temperature_ph(tp1907_target_enthalpy);
    let stanjan15_started = Instant::now();
    let stanjan15_range = solve_changing_temperature_stanjan15_ph(stanjan15_target_enthalpy);
    println!(
        "{}",
        Table::new([
            ph_probe_row("tp1907_nested", tp1907_range, tp1907_started),
            ph_probe_row("stanjan15_monolithic", stanjan15_range, stanjan15_started),
        ])
    );

    c.bench_function("tp1907_gibbs_temperature_sequence", |bencher| {
        bencher.iter(|| evaluate_sequence(black_box(&thermochemistry)));
    });
    c.bench_function("tp1907_changing_temperature_ph_range", |bencher| {
        bencher.iter(|| solve_changing_temperature_ph(black_box(tp1907_target_enthalpy)));
    });
    c.bench_function("stanjan15_changing_temperature_ph_range", |bencher| {
        bencher
            .iter(|| solve_changing_temperature_stanjan15_ph(black_box(stanjan15_target_enthalpy)));
    });
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
