//! PERF-PH-01: focused P,H preparation and continuation baseline.
//!
//! The existing facade baseline contains a broad matrix. This benchmark is
//! deliberately narrower and prints the P,H-specific owners of wall time:
//! monolithic formulation/RST preparation, nested trial/template preparation,
//! inner execution, bracket orchestration, and every accepted/rejected inner
//! backend attempt. It uses only the public equilibrium facade.
//!
//! Run in an optimized profile with:
//!
//! ```text
//! cargo bench --bench equilibrium_ph_preparation --features equilibrium-instrumentation -- --noplot
//! ```

use std::time::{Duration, Instant};

#[cfg(feature = "equilibrium-benchmark-allocation")]
use std::alloc::{GlobalAlloc, Layout, System};
#[cfg(feature = "equilibrium-benchmark-allocation")]
use std::sync::atomic::{AtomicU64, Ordering};

use criterion::{BenchmarkId, Criterion, black_box, criterion_group, criterion_main};
use tabled::{Table, Tabled};

use KiThe::Thermodynamics::ChemEquilibrium::prelude::*;

mod support;

use support::{
    PRESSURE_PA, REFERENCE_PRESSURE_PA, REFERENCE_TEMPERATURE_K, RP1311_40, STANJAN_15,
    TP1907_GAS_17, component, tp1907_element_seed,
};

#[cfg(feature = "equilibrium-benchmark-allocation")]
static ALLOCATION_CALLS: AtomicU64 = AtomicU64::new(0);
#[cfg(feature = "equilibrium-benchmark-allocation")]
static ALLOCATION_BYTES: AtomicU64 = AtomicU64::new(0);

#[cfg(feature = "equilibrium-benchmark-allocation")]
struct CountingAllocator;

#[cfg(feature = "equilibrium-benchmark-allocation")]
unsafe impl GlobalAlloc for CountingAllocator {
    unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
        ALLOCATION_CALLS.fetch_add(1, Ordering::Relaxed);
        ALLOCATION_BYTES.fetch_add(layout.size() as u64, Ordering::Relaxed);
        unsafe { System.alloc(layout) }
    }

    unsafe fn alloc_zeroed(&self, layout: Layout) -> *mut u8 {
        ALLOCATION_CALLS.fetch_add(1, Ordering::Relaxed);
        ALLOCATION_BYTES.fetch_add(layout.size() as u64, Ordering::Relaxed);
        unsafe { System.alloc_zeroed(layout) }
    }

    unsafe fn dealloc(&self, ptr: *mut u8, layout: Layout) {
        unsafe { System.dealloc(ptr, layout) };
    }

    unsafe fn realloc(&self, ptr: *mut u8, layout: Layout, new_size: usize) -> *mut u8 {
        ALLOCATION_CALLS.fetch_add(1, Ordering::Relaxed);
        ALLOCATION_BYTES.fetch_add(new_size as u64, Ordering::Relaxed);
        unsafe { System.realloc(ptr, layout, new_size) }
    }
}

#[cfg(feature = "equilibrium-benchmark-allocation")]
#[global_allocator]
static GLOBAL_ALLOCATOR: CountingAllocator = CountingAllocator;

#[derive(Clone, Copy)]
struct AllocationSnapshot {
    calls: u64,
    bytes: u64,
}

fn reset_allocations() {
    #[cfg(feature = "equilibrium-benchmark-allocation")]
    {
        ALLOCATION_CALLS.store(0, Ordering::Relaxed);
        ALLOCATION_BYTES.store(0, Ordering::Relaxed);
    }
}

fn allocation_snapshot() -> AllocationSnapshot {
    #[cfg(feature = "equilibrium-benchmark-allocation")]
    {
        return AllocationSnapshot {
            calls: ALLOCATION_CALLS.load(Ordering::Relaxed),
            bytes: ALLOCATION_BYTES.load(Ordering::Relaxed),
        };
    }
    #[cfg(not(feature = "equilibrium-benchmark-allocation"))]
    AllocationSnapshot { calls: 0, bytes: 0 }
}

fn allocation_mode() -> &'static str {
    if cfg!(feature = "equilibrium-benchmark-allocation") {
        "enabled"
    } else {
        "disabled"
    }
}

fn trial_cache_mode() -> &'static str {
    if cfg!(feature = "equilibrium-benchmark-no-trial-cache") {
        "disabled-baseline"
    } else {
        "enabled"
    }
}

#[derive(Clone, Copy, Debug)]
enum Case {
    SmallMolecular,
    Stanjan15Molecular,
    Stanjan15Elemental,
    Tp1907Elemental,
    Rp1311Elemental,
}

impl Case {
    fn label(self) -> &'static str {
        match self {
            Self::SmallMolecular => "small3_molecular",
            Self::Stanjan15Molecular => "stanjan15_molecular",
            Self::Stanjan15Elemental => "stanjan15_elemental",
            Self::Tp1907Elemental => "tp1907_18_elemental",
            Self::Rp1311Elemental => "rp1311_40_elemental",
        }
    }

    fn temperature(self) -> f64 {
        match self {
            Self::Tp1907Elemental => 700.0,
            _ => REFERENCE_TEMPERATURE_K,
        }
    }

    fn reference_pressure(self) -> f64 {
        match self {
            Self::Tp1907Elemental => PRESSURE_PA,
            _ => REFERENCE_PRESSURE_PA,
        }
    }

    fn bounds(self) -> TemperatureBounds {
        if matches!(self, Self::Tp1907Elemental) {
            return TemperatureBounds::new(699.0, 700.0).expect("benchmark bounds must validate");
        }
        let lower = (self.temperature() - 400.0).max(300.0);
        let upper = self.temperature() + 400.0;
        TemperatureBounds::new(lower, upper).expect("benchmark bounds must validate")
    }

    fn elemental(self) -> bool {
        !matches!(self, Self::SmallMolecular | Self::Stanjan15Molecular)
    }

    fn nested(self) -> bool {
        matches!(self, Self::Tp1907Elemental | Self::Rp1311Elemental)
    }

    fn prepared_target_delta(self) -> f64 {
        if matches!(self, Self::Tp1907Elemental) {
            -100.0
        } else {
            1_000.0
        }
    }
}

fn inventory(case: Case) -> ElementInventory {
    let values = match case {
        Case::SmallMolecular => unreachable!("molecular case has no inventory"),
        Case::Stanjan15Molecular => unreachable!("molecular case has no inventory"),
        Case::Stanjan15Elemental => vec![("C", 6.0), ("H", 16.0), ("O", 20.0), ("N", 75.2)],
        Case::Tp1907Elemental => vec![
            ("Ar", 0.05364816),
            ("C", 1.00182736),
            ("H", 2.0),
            ("N", 8.9462088),
            ("O", 2.40365472),
        ],
        Case::Rp1311Elemental => vec![("Ar", 0.1), ("C", 0.1), ("H", 2.2), ("O", 4.3), ("N", 6.0)],
    };
    ElementInventory::from_amounts(values).expect("benchmark inventory must validate")
}

fn base_builder(case: Case) -> EquilibriumCalculatorBuilder {
    match case {
        Case::SmallMolecular => EquilibriumCalculator::builder()
            .single_ideal_gas(["H2", "O2", "H2O"])
            .expect("small gas phase"),
        Case::Stanjan15Molecular | Case::Stanjan15Elemental => EquilibriumCalculator::builder()
            .ideal_gas_phase("gas", STANJAN_15)
            .expect("STANJAN phase"),
        Case::Tp1907Elemental => EquilibriumCalculator::builder()
            .ideal_gas_phase("gas", TP1907_GAS_17)
            .expect("TP-1907 gas phase")
            .pure_condensed_phase("graphite", "C(gr)", PhysicalState::Solid)
            .expect("graphite phase"),
        Case::Rp1311Elemental => EquilibriumCalculator::builder()
            .single_ideal_gas(RP1311_40)
            .expect("RP-1311 phase"),
    }
}

fn molecular_builder(case: Case) -> EquilibriumCalculatorBuilder {
    match case {
        Case::SmallMolecular => base_builder(case).initial_moles([0.1, 0.05, 1.9]),
        Case::Stanjan15Molecular => base_builder(case).initial_phase_moles([
            (component(Some("gas"), "CH4"), 4.0),
            (component(Some("gas"), "CO2"), 2.0),
            (component(Some("gas"), "O2"), 8.0),
            (component(Some("gas"), "N2"), 37.6),
        ]),
        _ => unreachable!("elemental case has no molecular seed"),
    }
}

fn pt_builder(case: Case) -> EquilibriumCalculatorBuilder {
    let builder = if case.elemental() {
        base_builder(case).element_inventory(inventory(case))
    } else {
        molecular_builder(case)
    };
    let builder = if matches!(case, Case::Tp1907Elemental) {
        builder.element_numerical_seed(tp1907_element_seed())
    } else {
        builder
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

fn reference_enthalpy(case: Case) -> f64 {
    if matches!(case, Case::Tp1907Elemental) {
        return -420_039.2761;
    }
    let outcome = pt_builder(case).solve().expect("reference P,T must solve");
    let EquilibriumCalculatorOutcome::PtPoint(point) = outcome else {
        panic!("reference calculation must be a P,T point")
    };
    let thermo = ResolvedThermochemistry::from_resolved_system(point.resolved())
        .expect("reference thermochemistry");
    thermo
        .enthalpy_model()
        .evaluate_total(point.solution().component_moles(), case.temperature())
        .expect("reference enthalpy must evaluate")
}

fn ph_builder(case: Case, target: f64) -> EquilibriumCalculatorBuilder {
    let builder = pt_builder(case);
    let builder = if case.nested() {
        builder.ph_solve_mode(PhSolveMode::NestedTemperature)
    } else {
        builder.ph_solve_mode(PhSolveMode::Monolithic)
    };
    builder.at_total_enthalpy(
        TotalEnthalpyJoules::new(target).expect("target must be finite"),
        case.temperature(),
        case.bounds(),
    )
}

fn prepared_ph_builder(case: Case, target: f64) -> EquilibriumCalculatorBuilder {
    ph_builder(case, target).ph_solve_mode(PhSolveMode::Monolithic)
}

fn production_policy_without_backtracking() -> SolverPolicy {
    let mut backends = RustedSciTheSolver::recommended_cascade()
        .into_iter()
        .filter(|backend| !matches!(backend, RustedSciTheSolver::BacktrackingLevenbergMarquardt))
        .map(SolverBackend::RustedSciThe)
        .collect::<Vec<_>>();
    backends.extend([
        SolverBackend::Legacy(LegacyEquilibriumSolver::LM),
        SolverBackend::Legacy(LegacyEquilibriumSolver::NR),
        SolverBackend::Legacy(LegacyEquilibriumSolver::TR),
    ]);
    SolverPolicy::Cascade(backends)
}

fn prepared_ph_builder_with_policy(
    case: Case,
    target: f64,
    include_backtracking: bool,
) -> EquilibriumCalculatorBuilder {
    if include_backtracking {
        prepared_ph_builder(case, target)
    } else {
        prepared_ph_builder(case, target)
            .solver_policy(production_policy_without_backtracking())
            .expect("comparison solver policy must validate")
    }
}

fn solve_point(case: Case, target: f64) -> FixedPressureEnthalpySolution {
    let outcome = ph_builder(case, target)
        .solve()
        .unwrap_or_else(|error| panic!("{} P,H solve failed: {error:?}", case.label()));
    let EquilibriumCalculatorOutcome::PhPoint(point) = outcome else {
        panic!("{} must return a P,H point", case.label())
    };
    point.solution().clone()
}

fn solve_prepared_route_point(case: Case, target: f64) -> FixedPressureEnthalpySolution {
    let outcome = prepared_ph_builder(case, target)
        .solve()
        .unwrap_or_else(|error| {
            panic!(
                "{} prepared-route P,H solve failed: {error:?}",
                case.label()
            )
        });
    let EquilibriumCalculatorOutcome::PhPoint(point) = outcome else {
        panic!("{} must return a P,H point", case.label())
    };
    point.solution().clone()
}

fn backend_name(report: &EquilibriumSolveReport) -> String {
    format!("{:?}", report.accepted_backend)
}

#[derive(Tabled)]
struct StageRow {
    case: &'static str,
    route: &'static str,
    wall_ms: String,
    ph_total_ms: String,
    trial_preparation_ms: String,
    trial_template_ms: String,
    monolithic_preparation_ms: String,
    monolithic_rst_ms: String,
    monolithic_backend_ms: String,
    inner_execution_ms: String,
    bracket_ms: String,
    builds: usize,
    reuses: usize,
    inner_attempts: usize,
    accepted_backend: String,
}

#[derive(Tabled)]
struct RecoveryCloneRow {
    case: &'static str,
    route: &'static str,
    wall_ms: String,
    ph_total_ms: String,
    recovery_clone_ms: String,
    clone_of_ph_pct: String,
    clone_of_wall_pct: String,
}

#[derive(Tabled)]
struct AttemptRow {
    case: &'static str,
    scope: String,
    index: usize,
    backend: String,
    outcome: String,
    wall_ms: String,
    iterations: usize,
    residual_calls: usize,
    jacobian_calls: usize,
}

#[derive(Tabled)]
struct RangeRow {
    case: &'static str,
    targets: usize,
    total_ms: String,
    mean_point_ms: String,
    builds: usize,
    reuses: usize,
    continuation_points: usize,
}

#[derive(Tabled)]
struct PreparedSessionRow {
    case: &'static str,
    allocation_mode: &'static str,
    one_shot_two_points_ms: String,
    preparation_ms: String,
    first_point_ms: String,
    second_point_ms: String,
    first_request_validation_ms: String,
    second_request_validation_ms: String,
    first_initial_enthalpy_ms: String,
    second_initial_enthalpy_ms: String,
    first_initial_enthalpy_reused: bool,
    second_initial_enthalpy_reused: bool,
    first_retarget_ms: String,
    second_retarget_ms: String,
    first_gibbs_ms: String,
    second_gibbs_ms: String,
    first_gibbs_reused: bool,
    second_gibbs_reused: bool,
    first_candidate_snapshots: usize,
    second_candidate_snapshots: usize,
    first_alloc_calls: u64,
    first_alloc_bytes: u64,
    second_alloc_calls: u64,
    second_alloc_bytes: u64,
    first_builds: usize,
    first_reuses: usize,
    second_builds: usize,
    second_reuses: usize,
}

#[derive(Tabled)]
struct PreparedRouteAvailabilityRow {
    case: &'static str,
    requested_route: &'static str,
    status: &'static str,
    reason: &'static str,
}

#[derive(Tabled)]
struct CascadeAttributionRow {
    case: &'static str,
    sample: usize,
    policy: &'static str,
    preparation_ms: String,
    first_point_ms: String,
    second_point_ms: String,
    total_ms: String,
    first_backend: String,
    second_backend: String,
    status: &'static str,
}

#[derive(Tabled)]
struct PreparationBreakdownRow {
    case: &'static str,
    facade_setup_ms: String,
    state_total_ms: String,
    thermochemistry_ms: String,
    normalization_ms: String,
    equation_ms: String,
    symbolic_ms: String,
    ph_formulation_ms: String,
    rst_ms: String,
    state_unattributed_ms: String,
}

#[derive(Tabled)]
struct NestedTrialRow {
    case: &'static str,
    trial: usize,
    temperature_k: String,
    step: String,
    preparation: String,
    template_build: bool,
    trial_total_ms: String,
    inner_execution_ms: String,
    inner_equilibrium_ms: String,
    enthalpy_ms: String,
    postprocessing_ms: String,
    equation_ms: String,
    rst_ms: String,
    nonlinear_ms: String,
    validation_ms: String,
}

#[derive(Tabled)]
struct NestedTemplatePreparationRow {
    case: &'static str,
    trial: usize,
    preparation: String,
    template_ms: String,
    facade_ms: String,
    request_ms: String,
    thermochemistry_ms: String,
    equation_ms: String,
    symbolic_ms: String,
    closure_ms: String,
    problem_ms: String,
    rst_ms: String,
    active_set_ms: String,
    projection_ms: String,
    unattributed_ms: String,
}

fn millis(duration: Duration) -> String {
    format!("{:.3}", duration.as_secs_f64() * 1_000.0)
}

fn percent_of(part: Duration, whole: Duration) -> String {
    let whole = whole.as_secs_f64();
    if whole <= 0.0 {
        return "n/a".to_owned();
    }
    format!("{:.1}%", part.as_secs_f64() * 100.0 / whole)
}

fn point_rows(
    case: Case,
    target: f64,
    stage_rows: &mut Vec<StageRow>,
    attempt_rows: &mut Vec<AttemptRow>,
    recovery_clone_rows: &mut Vec<RecoveryCloneRow>,
) {
    let started = Instant::now();
    let solution = solve_point(case, target);
    let wall = started.elapsed();
    let report = solution.report();
    let timing = report.timing();
    let route = if case.nested() {
        "nested"
    } else {
        "monolithic"
    };
    let accepted_backend = if let Some(evidence) = report.monolithic_evidence() {
        for (index, attempt) in evidence.solve_report().attempts.iter().enumerate() {
            attempt_rows.push(attempt_row(case, "monolithic", index, attempt));
        }
        backend_name(evidence.solve_report())
    } else {
        for (trial_index, trial) in report.trials().iter().enumerate() {
            if let Some(evidence) = trial.inner_evidence() {
                for (index, attempt) in evidence.solve_report().attempts.iter().enumerate() {
                    attempt_rows.push(attempt_row(
                        case,
                        &format!("trial-{trial_index}"),
                        index,
                        attempt,
                    ));
                }
            }
        }
        format!(
            "nested trials={}, attempts={}",
            report.trials().len(),
            report.inner_backend_attempts()
        )
    };
    stage_rows.push(StageRow {
        case: case.label(),
        route,
        wall_ms: millis(wall),
        ph_total_ms: millis(timing.total()),
        trial_preparation_ms: millis(timing.trial_preparation()),
        trial_template_ms: millis(timing.trial_template_preparation()),
        monolithic_preparation_ms: millis(timing.monolithic_preparation()),
        monolithic_rst_ms: millis(timing.monolithic_rst_preparation()),
        monolithic_backend_ms: millis(timing.monolithic_backend_execution()),
        inner_execution_ms: millis(timing.trial_inner_execution()),
        bracket_ms: millis(timing.bracket_orchestration()),
        builds: report.fixed_formulation_builds(),
        reuses: report.fixed_formulation_reuses(),
        inner_attempts: report.inner_backend_attempts(),
        accepted_backend,
    });
    recovery_clone_rows.push(RecoveryCloneRow {
        case: case.label(),
        route,
        wall_ms: millis(wall),
        ph_total_ms: millis(timing.total()),
        recovery_clone_ms: millis(timing.recovery_request_clone()),
        clone_of_ph_pct: percent_of(timing.recovery_request_clone(), timing.total()),
        clone_of_wall_pct: percent_of(timing.recovery_request_clone(), wall),
    });
}

fn attempt_row(case: Case, scope: &str, index: usize, attempt: &SolverAttemptReport) -> AttemptRow {
    let metrics = attempt.metrics.as_ref();
    AttemptRow {
        case: case.label(),
        scope: scope.to_owned(),
        index,
        backend: format!("{:?}", attempt.backend),
        outcome: format!("{:?}", attempt.outcome),
        wall_ms: metrics.map_or_else(
            || "n/a".to_owned(),
            |item| millis(Duration::from_micros(item.elapsed_micros as u64)),
        ),
        iterations: metrics.map_or(0, |item| item.iterations),
        residual_calls: metrics.map_or(0, |item| item.residual_evaluations),
        jacobian_calls: metrics.map_or(0, |item| item.jacobian_evaluations),
    }
}

fn range_row(case: Case, target: f64) -> RangeRow {
    let targets = if case.nested() {
        vec![target - 100.0, target]
    } else {
        vec![target, target + 1_000.0, target + 2_000.0]
    };
    let started = Instant::now();
    let outcome = ph_builder(case, target)
        .over_enthalpy_range(targets, case.temperature(), case.bounds())
        .unwrap_or_else(|error| panic!("{} range setup failed: {error:?}", case.label()))
        .solve()
        .unwrap_or_else(|error| panic!("{} range failed: {error:?}", case.label()));
    let EquilibriumCalculatorOutcome::PhRange(range) = outcome else {
        panic!("{} must return a P,H range", case.label())
    };
    let report = range.solution().report();
    RangeRow {
        case: case.label(),
        targets: report.point_count(),
        total_ms: millis(started.elapsed()),
        mean_point_ms: millis(report.point_timing().mean()),
        builds: report.formulation_builds(),
        reuses: report.formulation_reuses(),
        continuation_points: report.continuation_points(),
    }
}

fn prepared_session_row(case: Case, target: f64) -> PreparedSessionRow {
    let target_delta = case.prepared_target_delta();
    let one_shot_started = Instant::now();
    let _ = solve_prepared_route_point(case, target);
    let _ = solve_prepared_route_point(case, target + target_delta);
    let one_shot_two_points_ms = millis(one_shot_started.elapsed());
    let preparation_started = Instant::now();
    let mut session = prepared_ph_builder(case, target)
        .prepare_ph_session()
        .unwrap_or_else(|error| panic!("{} session setup failed: {error:?}", case.label()));
    let preparation_ms = millis(preparation_started.elapsed());

    let first_started = Instant::now();
    reset_allocations();
    let first = session
        .solve_at_total_enthalpy(TotalEnthalpyJoules::new(target).expect("target must be finite"))
        .unwrap_or_else(|error| panic!("{} first session point failed: {error:?}", case.label()));
    let first_point_ms = millis(first_started.elapsed());
    let first_allocations = allocation_snapshot();

    let second_started = Instant::now();
    reset_allocations();
    let second = session
        .solve_at_total_enthalpy(
            TotalEnthalpyJoules::new(target + target_delta).expect("target must be finite"),
        )
        .unwrap_or_else(|error| panic!("{} second session point failed: {error:?}", case.label()));
    let second_point_ms = millis(second_started.elapsed());
    let second_allocations = allocation_snapshot();
    let first_timing = first.solution().report().timing();
    let second_timing = second.solution().report().timing();
    PreparedSessionRow {
        case: case.label(),
        allocation_mode: allocation_mode(),
        one_shot_two_points_ms,
        preparation_ms,
        first_point_ms,
        second_point_ms,
        first_request_validation_ms: millis(first_timing.request_validation()),
        second_request_validation_ms: millis(second_timing.request_validation()),
        first_initial_enthalpy_ms: millis(first_timing.initial_enthalpy()),
        second_initial_enthalpy_ms: millis(second_timing.initial_enthalpy()),
        first_initial_enthalpy_reused: first_timing.initial_enthalpy_reused(),
        second_initial_enthalpy_reused: second_timing.initial_enthalpy_reused(),
        first_retarget_ms: millis(first_timing.prepared_problem_retarget()),
        second_retarget_ms: millis(second_timing.prepared_problem_retarget()),
        first_gibbs_ms: millis(first_timing.gibbs_snapshot()),
        second_gibbs_ms: millis(second_timing.gibbs_snapshot()),
        first_gibbs_reused: first_timing.gibbs_snapshot_reused(),
        second_gibbs_reused: second_timing.gibbs_snapshot_reused(),
        first_candidate_snapshots: first
            .instrumentation_report()
            .counters()
            .candidate_snapshot_evaluations(),
        second_candidate_snapshots: second
            .instrumentation_report()
            .counters()
            .candidate_snapshot_evaluations(),
        first_alloc_calls: first_allocations.calls,
        first_alloc_bytes: first_allocations.bytes,
        second_alloc_calls: second_allocations.calls,
        second_alloc_bytes: second_allocations.bytes,
        first_builds: first.solution().report().fixed_formulation_builds(),
        first_reuses: first.solution().report().fixed_formulation_reuses(),
        second_builds: second.solution().report().fixed_formulation_builds(),
        second_reuses: second.solution().report().fixed_formulation_reuses(),
    }
}

fn cascade_attribution_row(
    case: Case,
    target: f64,
    include_backtracking: bool,
    sample: usize,
) -> CascadeAttributionRow {
    let policy = if include_backtracking {
        "production + Backtracking LM"
    } else {
        "production without Backtracking LM"
    };
    let preparation_started = Instant::now();
    let mut session = prepared_ph_builder_with_policy(case, target, include_backtracking)
        .prepare_ph_session()
        .unwrap_or_else(|error| panic!("{} comparison setup failed: {error:?}", case.label()));
    let preparation = preparation_started.elapsed();

    let first_started = Instant::now();
    let first = session
        .solve_at_total_enthalpy(TotalEnthalpyJoules::new(target).expect("target must be finite"));
    let first_point = first_started.elapsed();
    let first_backend = first
        .as_ref()
        .ok()
        .and_then(|outcome| outcome.solution().report().monolithic_evidence())
        .map(|evidence| backend_name(evidence.solve_report()))
        .unwrap_or_else(|| "error".to_owned());

    let second_started = Instant::now();
    let second = if first.is_ok() {
        Some(
            session.solve_at_total_enthalpy(
                TotalEnthalpyJoules::new(target + case.prepared_target_delta())
                    .expect("target must be finite"),
            ),
        )
    } else {
        None
    };
    let second_point = second_started.elapsed();
    let second_backend = second
        .as_ref()
        .and_then(|result| result.as_ref().ok())
        .and_then(|outcome| outcome.solution().report().monolithic_evidence())
        .map(|evidence| backend_name(evidence.solve_report()))
        .unwrap_or_else(|| "error".to_owned());
    let status = if first.is_ok() && second.as_ref().is_some_and(Result::is_ok) {
        "accepted"
    } else {
        "failed"
    };
    CascadeAttributionRow {
        case: case.label(),
        sample,
        policy,
        preparation_ms: millis(preparation),
        first_point_ms: millis(first_point),
        second_point_ms: millis(second_point),
        total_ms: millis(preparation + first_point + second_point),
        first_backend,
        second_backend,
        status,
    }
}

fn preparation_breakdown_row(case: Case, target: f64) -> PreparationBreakdownRow {
    let started = Instant::now();
    let session = prepared_ph_builder(case, target)
        .prepare_ph_session()
        .unwrap_or_else(|error| panic!("{} session setup failed: {error:?}", case.label()));
    let wall = started.elapsed();
    let timing = session.preparation_timing();
    let known = [
        timing.request_preparation(),
        timing.normalization(),
        timing.thermochemistry_preparation(),
        timing.equation_construction(),
        timing.ph_formulation_preparation(),
        timing.rst_preparation(),
        timing.numeric_closure_construction(),
        timing.symbolic_construction(),
        timing.numerical_problem_preparation(),
        timing.active_set_formulation_build(),
        timing.projection_build(),
    ]
    .into_iter()
    .fold(Duration::ZERO, |total, duration| total + duration);
    PreparationBreakdownRow {
        case: case.label(),
        facade_setup_ms: millis(wall.saturating_sub(timing.total())),
        state_total_ms: millis(timing.total()),
        thermochemistry_ms: millis(timing.thermochemistry_preparation()),
        normalization_ms: millis(timing.normalization()),
        equation_ms: millis(timing.equation_construction()),
        symbolic_ms: millis(timing.symbolic_construction()),
        ph_formulation_ms: millis(timing.ph_formulation_preparation()),
        rst_ms: millis(timing.rst_preparation()),
        state_unattributed_ms: millis(timing.total().saturating_sub(known)),
    }
}

fn nested_trial_rows(
    case: Case,
    target: f64,
    rows: &mut Vec<NestedTrialRow>,
    preparation_rows: &mut Vec<NestedTemplatePreparationRow>,
) {
    let solution = solve_point(case, target);
    for (index, trial) in solution.report().trials().iter().enumerate() {
        let timing = trial.timing();
        let inner = trial.inner_timing();
        rows.push(NestedTrialRow {
            case: case.label(),
            trial: index,
            temperature_k: format!("{:.3}", trial.temperature()),
            step: format!("{:?}", trial.step_kind()),
            preparation: format!("{:?}", trial.preparation()),
            template_build: matches!(
                trial.preparation(),
                PhTrialPreparation::FixedFormulationInitial
            ),
            trial_total_ms: millis(timing.total()),
            inner_execution_ms: millis(timing.inner_execution()),
            inner_equilibrium_ms: millis(timing.inner_equilibrium()),
            enthalpy_ms: millis(timing.enthalpy_evaluation()),
            postprocessing_ms: millis(timing.postprocessing()),
            equation_ms: millis(inner.equation_construction()),
            rst_ms: millis(inner.rst_preparation()),
            nonlinear_ms: millis(inner.nonlinear_solve()),
            validation_ms: millis(inner.validation()),
        });
        let template = timing.template_preparation();
        let known = [
            template.facade_resolution(),
            template.request_preparation(),
            template.thermochemistry_preparation(),
            template.ph_formulation_preparation(),
            template.numeric_closure_construction(),
            template.symbolic_construction(),
            template.equation_construction(),
            template.numerical_problem_preparation(),
            template.rst_preparation(),
            template.active_set_formulation_build(),
            template.projection_build(),
            template.validation(),
            template.postprocessing(),
        ]
        .into_iter()
        .fold(Duration::ZERO, |total, duration| total + duration);
        preparation_rows.push(NestedTemplatePreparationRow {
            case: case.label(),
            trial: index,
            preparation: format!("{:?}", trial.preparation()),
            template_ms: millis(template.total()),
            facade_ms: millis(template.facade_resolution()),
            request_ms: millis(template.request_preparation()),
            thermochemistry_ms: millis(template.thermochemistry_preparation()),
            equation_ms: millis(template.equation_construction()),
            symbolic_ms: millis(template.symbolic_construction()),
            closure_ms: millis(template.numeric_closure_construction()),
            problem_ms: millis(template.numerical_problem_preparation()),
            rst_ms: millis(template.rst_preparation()),
            active_set_ms: millis(template.active_set_formulation_build()),
            projection_ms: millis(template.projection_build()),
            unattributed_ms: millis(template.total().saturating_sub(known)),
        });
    }
}

fn cases() -> [Case; 5] {
    [
        Case::SmallMolecular,
        Case::Stanjan15Molecular,
        Case::Stanjan15Elemental,
        Case::Tp1907Elemental,
        Case::Rp1311Elemental,
    ]
}

fn print_capture() {
    let mut stages = Vec::new();
    let mut attempts = Vec::new();
    let mut recovery_clone_rows = Vec::new();
    for case in cases() {
        point_rows(
            case,
            reference_enthalpy(case),
            &mut stages,
            &mut attempts,
            &mut recovery_clone_rows,
        );
    }
    let ranges = [
        Case::SmallMolecular,
        Case::Stanjan15Elemental,
        Case::Rp1311Elemental,
    ]
    .into_iter()
    .map(|case| range_row(case, reference_enthalpy(case)))
    .collect::<Vec<_>>();
    let prepared_sessions = [
        Case::SmallMolecular,
        Case::Stanjan15Molecular,
        Case::Stanjan15Elemental,
        Case::Rp1311Elemental,
    ]
    .into_iter()
    .map(|case| prepared_session_row(case, reference_enthalpy(case)))
    .collect::<Vec<_>>();
    let preparation_breakdown = [
        Case::SmallMolecular,
        Case::Stanjan15Molecular,
        Case::Stanjan15Elemental,
        Case::Rp1311Elemental,
    ]
    .into_iter()
    .map(|case| preparation_breakdown_row(case, reference_enthalpy(case)))
    .collect::<Vec<_>>();
    let prepared_route_availability = [
        PreparedRouteAvailabilityRow {
            case: Case::SmallMolecular.label(),
            requested_route: "monolithic",
            status: "measured",
            reason: "prepared session is a production route",
        },
        PreparedRouteAvailabilityRow {
            case: Case::Stanjan15Molecular.label(),
            requested_route: "monolithic",
            status: "measured",
            reason: "prepared session is a production route",
        },
        PreparedRouteAvailabilityRow {
            case: Case::Stanjan15Elemental.label(),
            requested_route: "monolithic",
            status: "measured",
            reason: "prepared session is a production route",
        },
        PreparedRouteAvailabilityRow {
            case: Case::Tp1907Elemental.label(),
            requested_route: "monolithic",
            status: "N/A",
            reason: "current production route is nested; fixed route rejects this case",
        },
        PreparedRouteAvailabilityRow {
            case: Case::Rp1311Elemental.label(),
            requested_route: "monolithic",
            status: "measured",
            reason: "prepared session is a production route",
        },
    ];
    let cascade_attribution = [
        Case::SmallMolecular,
        Case::Stanjan15Molecular,
        Case::Stanjan15Elemental,
        Case::Rp1311Elemental,
    ]
    .into_iter()
    .flat_map(|case| {
        let target = reference_enthalpy(case);
        (1..=3).flat_map(move |sample| {
            [
                cascade_attribution_row(case, target, true, sample),
                cascade_attribution_row(case, target, false, sample),
            ]
        })
    })
    .collect::<Vec<_>>();
    let mut nested_trials = Vec::new();
    let mut nested_template_preparation = Vec::new();
    for case in [Case::Tp1907Elemental, Case::Rp1311Elemental] {
        nested_trial_rows(
            case,
            reference_enthalpy(case),
            &mut nested_trials,
            &mut nested_template_preparation,
        );
    }
    println!("\nPERF-PH-01 P,H point stage ownership (inclusive fields)");
    println!("{}", Table::new(stages));
    println!("\nPERF-PH-01 all inner backend attempts");
    println!("{}", Table::new(attempts));
    println!("\nPERF-PH-17 recovery request clone attribution");
    println!("{}", Table::new(recovery_clone_rows));
    println!("\nPERF-PH-01 P,H range preparation reuse");
    println!("{}", Table::new(ranges));
    println!("\nPERF-PH-09/PERF-PH-11 prepared monolithic P,H retarget detail");
    println!("{}", Table::new(&prepared_sessions));
    println!("\nPERF-PH-09 prepared-route availability");
    println!("{}", Table::new(prepared_route_availability));
    println!("\nPERF-PH-14 production cascade attribution");
    println!("{}", Table::new(cascade_attribution));
    println!("\nPERF-PH-03 prepared monolithic cold preparation breakdown");
    println!("{}", Table::new(preparation_breakdown));
    println!(
        "\nPERF-PH-03 nested P,H trial detail (template cost is aggregate in the stage table)"
    );
    println!("{}", Table::new(nested_trials));
    println!("\nPERF-PH-18 nested fixed-template preparation detail");
    println!("{}", Table::new(nested_template_preparation));
    println!("\nPERF-PH-12/PERF-PH-13 compact changing-target evidence");
    println!("trial_cache_mode={}", trial_cache_mode());
    for row in &prepared_sessions {
        println!(
            "case={} preparation_ms={} first_point_ms={} second_point_ms={} first_snapshots={} second_snapshots={} first_alloc_bytes={} second_alloc_bytes={}",
            row.case,
            row.preparation_ms,
            row.first_point_ms,
            row.second_point_ms,
            row.first_candidate_snapshots,
            row.second_candidate_snapshots,
            row.first_alloc_bytes,
            row.second_alloc_bytes,
        );
    }
}

fn criterion_benchmark(c: &mut Criterion) {
    print_capture();
    let mut group = c.benchmark_group("equilibrium_ph_preparation");
    group.sample_size(10);
    for case in cases() {
        let target = reference_enthalpy(case);
        group.bench_function(BenchmarkId::new("point", case.label()), |bencher| {
            bencher.iter(|| black_box(solve_point(case, target)));
        });
    }
    for case in [
        Case::SmallMolecular,
        Case::Stanjan15Elemental,
        Case::Rp1311Elemental,
    ] {
        let target = reference_enthalpy(case);
        group.bench_function(BenchmarkId::new("range", case.label()), |bencher| {
            bencher.iter(|| black_box(range_row(case, target)));
        });
    }
    for case in [
        Case::SmallMolecular,
        Case::Stanjan15Molecular,
        Case::Stanjan15Elemental,
        Case::Rp1311Elemental,
    ] {
        let target = reference_enthalpy(case);
        let target_delta = case.prepared_target_delta();
        group.bench_function(
            BenchmarkId::new("one_shot_two_points", case.label()),
            |bencher| {
                bencher.iter(|| {
                    black_box(solve_prepared_route_point(case, target));
                    black_box(solve_prepared_route_point(case, target + target_delta));
                });
            },
        );
        group.bench_function(
            BenchmarkId::new("prepared_session_two_points", case.label()),
            |bencher| {
                bencher.iter(|| {
                    let mut session = prepared_ph_builder(case, target)
                        .prepare_ph_session()
                        .expect("prepared P,H session must build");
                    black_box(
                        session
                            .solve_at_total_enthalpy(TotalEnthalpyJoules::new(target).unwrap())
                            .expect("first prepared P,H point must solve"),
                    );
                    black_box(
                        session
                            .solve_at_total_enthalpy(
                                TotalEnthalpyJoules::new(target + target_delta).unwrap(),
                            )
                            .expect("second prepared P,H point must solve"),
                    );
                });
            },
        );
        group.bench_function(
            BenchmarkId::new("prepared_session_hot_point", case.label()),
            |bencher| {
                let target_delta = case.prepared_target_delta();
                let mut session = prepared_ph_builder(case, target)
                    .prepare_ph_session()
                    .expect("prepared P,H session must build");
                session
                    .solve_at_total_enthalpy(TotalEnthalpyJoules::new(target).unwrap())
                    .expect("session seed point must solve");
                bencher.iter(|| {
                    black_box(
                        session
                            .solve_at_total_enthalpy(
                                TotalEnthalpyJoules::new(target + target_delta).unwrap(),
                            )
                            .expect("hot prepared P,H point must solve"),
                    );
                });
            },
        );
    }
    group.finish();
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
