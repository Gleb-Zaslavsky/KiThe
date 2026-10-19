//! Opt-in evaluation counters for equilibrium characterization.
//!
//! Stage durations remain owned by [`super::equilibrium_timing`]. This module
//! owns fine-grained evaluation counts that are useful when deciding whether
//! a formulation, callback, or cache does unnecessary work. The normal library
//! build compiles these hooks out; enable the `equilibrium-instrumentation`
//! feature, or run tests, when a benchmark needs the counters.

use RustedSciThe::numerical::SolverTelemetryMode;
#[cfg(any(test, feature = "equilibrium-instrumentation"))]
use std::cell::RefCell;
#[cfg(any(test, feature = "equilibrium-instrumentation"))]
use std::rc::Rc;
#[cfg(any(test, feature = "equilibrium-instrumentation"))]
use std::sync::{Arc, Mutex};
use std::time::Duration;

/// Fine-grained thermochemistry and monolithic P,H evaluation counts.
///
/// The fields are intentionally private so new counters can be added without
/// making the diagnostic payload part of the numerical API. Use the accessors
/// when rendering a benchmark report.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct EquilibriumEvaluationCounters {
    context_evaluations: usize,
    residual_context_evaluations: usize,
    jacobian_context_evaluations: usize,
    residual_evaluations: usize,
    jacobian_evaluations: usize,
    candidate_snapshot_evaluations: usize,
    gibbs_vector_evaluations: usize,
    enthalpy_vector_evaluations: usize,
    heat_capacity_vector_evaluations: usize,
    gibbs_component_calls: usize,
    enthalpy_component_calls: usize,
    heat_capacity_component_calls: usize,
    gibbs_coefficient_refreshes: usize,
    gibbs_interval_reuses: usize,
    property_coefficient_selections: usize,
    property_calculations: usize,
    phase_equilibrium_builds: usize,
}

/// Controls automatic fine-grained instrumentation at the facade boundary.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum EquilibriumInstrumentationMode {
    /// Keep counter hooks out of the ordinary solve path.
    #[default]
    Disabled,
    /// Capture evaluation counts for the enclosing facade operation.
    Enabled,
}

/// Resolves calculator diagnostics into the matching RST collection level.
/// RST counters follow the active facade instrumentation capture, while its
/// clock reads follow the independently configured stage-timing mode.
pub(crate) fn rst_solver_telemetry_mode(timing_enabled: bool) -> SolverTelemetryMode {
    let counters_enabled = instrumentation_capture_active();
    match (counters_enabled, timing_enabled) {
        (false, false) => SolverTelemetryMode::Off,
        (true, false) => SolverTelemetryMode::Counters,
        (false, true) => SolverTelemetryMode::Timings,
        (true, true) => SolverTelemetryMode::Detailed,
    }
}

fn instrumentation_capture_active() -> bool {
    #[cfg(any(test, feature = "equilibrium-instrumentation"))]
    {
        return ACTIVE_COUNTERS.with(|active| active.borrow().is_some());
    }
    #[cfg(not(any(test, feature = "equilibrium-instrumentation")))]
    {
        false
    }
}

impl EquilibriumEvaluationCounters {
    pub fn context_evaluations(&self) -> usize {
        self.context_evaluations
    }

    pub fn residual_context_evaluations(&self) -> usize {
        self.residual_context_evaluations
    }

    pub fn jacobian_context_evaluations(&self) -> usize {
        self.jacobian_context_evaluations
    }

    pub fn residual_evaluations(&self) -> usize {
        self.residual_evaluations
    }

    pub fn jacobian_evaluations(&self) -> usize {
        self.jacobian_evaluations
    }

    pub fn candidate_snapshot_evaluations(&self) -> usize {
        self.candidate_snapshot_evaluations
    }

    pub fn gibbs_vector_evaluations(&self) -> usize {
        self.gibbs_vector_evaluations
    }

    pub fn enthalpy_vector_evaluations(&self) -> usize {
        self.enthalpy_vector_evaluations
    }

    pub fn heat_capacity_vector_evaluations(&self) -> usize {
        self.heat_capacity_vector_evaluations
    }

    pub fn gibbs_component_calls(&self) -> usize {
        self.gibbs_component_calls
    }

    pub fn enthalpy_component_calls(&self) -> usize {
        self.enthalpy_component_calls
    }

    pub fn heat_capacity_component_calls(&self) -> usize {
        self.heat_capacity_component_calls
    }

    pub fn gibbs_coefficient_refreshes(&self) -> usize {
        self.gibbs_coefficient_refreshes
    }

    pub fn gibbs_interval_reuses(&self) -> usize {
        self.gibbs_interval_reuses
    }

    pub fn property_coefficient_selections(&self) -> usize {
        self.property_coefficient_selections
    }

    pub fn property_calculations(&self) -> usize {
        self.property_calculations
    }

    /// Number of canonical bridge formulation builds in the captured scope.
    pub fn phase_equilibrium_builds(&self) -> usize {
        self.phase_equilibrium_builds
    }
}

/// Snapshot returned by one explicit instrumentation capture.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct EquilibriumInstrumentationReport {
    enabled: bool,
    elapsed: Duration,
    counters: EquilibriumEvaluationCounters,
}

/// Per-attempt diagnostics produced by the retained legacy nonlinear
/// algorithms. The adapter copies this snapshot into `SolverAttemptMetrics`,
/// so the accepted and rejected backend trace remains the source of truth.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub(crate) struct LegacySolverDiagnosticsSnapshot {
    pub residual_calls: usize,
    pub jacobian_calls: usize,
    pub same_iterate_callback_pairs: usize,
    pub linear_system_solves: usize,
    pub residual_micros: u128,
    pub jacobian_micros: u128,
    pub linear_system_micros: u128,
    pub matrix_operations_micros: u128,
    pub norm_micros: u128,
    pub trial_vector_micros: u128,
    pub feasibility_micros: u128,
    pub globalization_micros: u128,
}

/// Coarse operations inside a retained Legacy nonlinear attempt.
///
/// These timings are deliberately solver-facing rather than thermochemistry-
/// facing: residual and Jacobian callback time is reported separately, while
/// this enum explains the remainder of the Legacy backend wall time.
#[derive(Debug, Clone, Copy)]
pub(crate) enum LegacySolverStage {
    MatrixOperations,
    Norm,
    TrialVector,
    Feasibility,
    Globalization,
}

#[cfg(any(test, feature = "equilibrium-instrumentation"))]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum LegacyCallbackKind {
    Residual,
    Jacobian,
}

#[cfg(any(test, feature = "equilibrium-instrumentation"))]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct IterateSignature {
    len: usize,
    first: u64,
    second: u64,
}

#[cfg(any(test, feature = "equilibrium-instrumentation"))]
#[derive(Debug, Clone, Copy)]
struct PreviousLegacyCallback {
    kind: LegacyCallbackKind,
    signature: IterateSignature,
}

#[cfg(any(test, feature = "equilibrium-instrumentation"))]
#[derive(Debug, Default)]
pub(crate) struct LegacySolverDiagnosticsState {
    snapshot: LegacySolverDiagnosticsSnapshot,
    previous_callback: Option<PreviousLegacyCallback>,
}

#[cfg(any(test, feature = "equilibrium-instrumentation"))]
thread_local! {
    static ACTIVE_LEGACY_SOLVER_DIAGNOSTICS:
        RefCell<Option<Rc<RefCell<LegacySolverDiagnosticsState>>>> =
        const { RefCell::new(None) };
}

/// Restores a nested legacy diagnostic scope when the backend attempt ends.
#[cfg(any(test, feature = "equilibrium-instrumentation"))]
#[derive(Debug, Default)]
pub(crate) struct LegacySolverDiagnosticsGuard {
    #[cfg(any(test, feature = "equilibrium-instrumentation"))]
    previous: Option<Rc<RefCell<LegacySolverDiagnosticsState>>>,
}

/// Opaque per-attempt session used by the legacy backend adapter.
pub(crate) struct LegacySolverDiagnosticsSession {
    #[cfg(any(test, feature = "equilibrium-instrumentation"))]
    guard: LegacySolverDiagnosticsGuard,
    #[cfg(any(test, feature = "equilibrium-instrumentation"))]
    state: Rc<RefCell<LegacySolverDiagnosticsState>>,
}

#[cfg(any(test, feature = "equilibrium-instrumentation"))]
impl Drop for LegacySolverDiagnosticsGuard {
    fn drop(&mut self) {
        ACTIVE_LEGACY_SOLVER_DIAGNOSTICS.with(|active| {
            *active.borrow_mut() = self.previous.take();
        });
    }
}

#[cfg(any(test, feature = "equilibrium-instrumentation"))]
fn legacy_solver_diagnostics_enabled() -> bool {
    #[cfg(any(test, feature = "equilibrium-instrumentation"))]
    {
        return ACTIVE_LEGACY_SOLVER_DIAGNOSTICS.with(|active| active.borrow().is_some());
    }

    #[cfg(not(any(test, feature = "equilibrium-instrumentation")))]
    {
        false
    }
}

#[cfg(any(test, feature = "equilibrium-instrumentation"))]
pub(crate) fn legacy_solver_timing_enabled() -> bool {
    legacy_solver_diagnostics_enabled()
}

#[cfg(any(test, feature = "equilibrium-instrumentation"))]
pub(crate) fn record_legacy_residual(duration: Duration) {
    #[cfg(any(test, feature = "equilibrium-instrumentation"))]
    ACTIVE_LEGACY_SOLVER_DIAGNOSTICS.with(|active| {
        if let Some(state) = active.borrow().as_ref() {
            let mut state = state.borrow_mut();
            state.snapshot.residual_calls += 1;
            state.snapshot.residual_micros += duration.as_micros();
        }
    });
}

#[cfg(any(test, feature = "equilibrium-instrumentation"))]
fn iterate_signature(values: &[f64]) -> IterateSignature {
    let mut first = 0xcbf29ce484222325_u64;
    let mut second = 0x9e3779b97f4a7c15_u64;
    for (index, value) in values.iter().enumerate() {
        let bits = value.to_bits().wrapping_add(index as u64);
        first ^= bits;
        first = first.wrapping_mul(0x100000001b3_u64);
        second ^= bits.rotate_left((index & 63) as u32);
        second = second.rotate_left(27).wrapping_mul(0x94d049bb133111eb_u64);
    }
    IterateSignature {
        len: values.len(),
        first,
        second,
    }
}

/// Records adjacent Legacy residual/Jacobian callbacks that received the
/// same exact floating-point iterate. The diagnostic signature avoids cloning
/// every callback input while retaining a collision-resistant exact-bit key.
#[cfg(any(test, feature = "equilibrium-instrumentation"))]
pub(crate) fn record_legacy_callback(kind: LegacyCallbackKind, values: &[f64]) {
    let signature = iterate_signature(values);
    ACTIVE_LEGACY_SOLVER_DIAGNOSTICS.with(|active| {
        if let Some(state) = active.borrow().as_ref() {
            let mut state = state.borrow_mut();
            if let Some(previous) = state.previous_callback {
                if previous.kind != kind && previous.signature == signature {
                    state.snapshot.same_iterate_callback_pairs += 1;
                }
            }
            state.previous_callback = Some(PreviousLegacyCallback { kind, signature });
        }
    });
}

#[cfg(any(test, feature = "equilibrium-instrumentation"))]
pub(crate) fn record_legacy_jacobian(duration: Duration) {
    #[cfg(any(test, feature = "equilibrium-instrumentation"))]
    ACTIVE_LEGACY_SOLVER_DIAGNOSTICS.with(|active| {
        if let Some(state) = active.borrow().as_ref() {
            let mut state = state.borrow_mut();
            state.snapshot.jacobian_calls += 1;
            state.snapshot.jacobian_micros += duration.as_micros();
        }
    });
}

#[cfg(any(test, feature = "equilibrium-instrumentation"))]
pub(crate) fn record_legacy_linear_system(duration: Duration) {
    #[cfg(any(test, feature = "equilibrium-instrumentation"))]
    ACTIVE_LEGACY_SOLVER_DIAGNOSTICS.with(|active| {
        if let Some(state) = active.borrow().as_ref() {
            let mut state = state.borrow_mut();
            state.snapshot.linear_system_solves += 1;
            state.snapshot.linear_system_micros += duration.as_micros();
        }
    });
}

#[cfg(any(test, feature = "equilibrium-instrumentation"))]
pub(crate) fn record_legacy_stage(stage: LegacySolverStage, duration: Duration) {
    ACTIVE_LEGACY_SOLVER_DIAGNOSTICS.with(|active| {
        if let Some(state) = active.borrow().as_ref() {
            let mut state = state.borrow_mut();
            let field = match stage {
                LegacySolverStage::MatrixOperations => &mut state.snapshot.matrix_operations_micros,
                LegacySolverStage::Norm => &mut state.snapshot.norm_micros,
                LegacySolverStage::TrialVector => &mut state.snapshot.trial_vector_micros,
                LegacySolverStage::Feasibility => &mut state.snapshot.feasibility_micros,
                LegacySolverStage::Globalization => &mut state.snapshot.globalization_micros,
            };
            *field += duration.as_micros();
        }
    });
}

/// Opens one legacy solver scope only inside an enabled calculator capture.
pub(crate) fn begin_legacy_solver_diagnostics() -> Option<LegacySolverDiagnosticsSession> {
    #[cfg(any(test, feature = "equilibrium-instrumentation"))]
    {
        let capture_active = ACTIVE_COUNTERS.with(|active| active.borrow().is_some());
        if !capture_active {
            return None;
        }
        let state = Rc::new(RefCell::new(LegacySolverDiagnosticsState::default()));
        let previous = ACTIVE_LEGACY_SOLVER_DIAGNOSTICS
            .with(|active| active.borrow_mut().replace(Rc::clone(&state)));
        return Some(LegacySolverDiagnosticsSession {
            guard: LegacySolverDiagnosticsGuard { previous },
            state,
        });
    }

    #[cfg(not(any(test, feature = "equilibrium-instrumentation")))]
    {
        None
    }
}

pub(crate) fn finish_legacy_solver_diagnostics(
    session: LegacySolverDiagnosticsSession,
) -> LegacySolverDiagnosticsSnapshot {
    #[cfg(any(test, feature = "equilibrium-instrumentation"))]
    {
        let LegacySolverDiagnosticsSession { guard, state } = session;
        drop(guard);
        return state.borrow().snapshot;
    }

    #[cfg(not(any(test, feature = "equilibrium-instrumentation")))]
    {
        let _ = session;
        LegacySolverDiagnosticsSnapshot::default()
    }
}

impl EquilibriumInstrumentationReport {
    /// Returns the empty report used when instrumentation is disabled.
    pub fn disabled() -> Self {
        Self::default()
    }

    /// Whether this build was compiled with active instrumentation hooks.
    pub fn enabled(&self) -> bool {
        self.enabled
    }

    /// Wall-clock duration of the explicitly captured operation.
    ///
    /// Stage-level timing remains available through `EquilibriumTimingReport`;
    /// this value is useful for correlating counters with an outer benchmark.
    pub fn elapsed(&self) -> Duration {
        self.elapsed
    }

    pub fn counters(&self) -> EquilibriumEvaluationCounters {
        self.counters
    }
}

#[derive(Debug, Clone, Copy)]
enum Counter {
    Context,
    ResidualContext,
    JacobianContext,
    Residual,
    Jacobian,
    CandidateSnapshot,
    GibbsVector,
    EnthalpyVector,
    HeatCapacityVector,
    GibbsComponent,
    EnthalpyComponent,
    HeatCapacityComponent,
    GibbsCoefficientRefresh,
    GibbsIntervalReuse,
    PropertyCoefficientSelection,
    PropertyCalculation,
    PhaseEquilibriumBuild,
}

#[cfg(any(test, feature = "equilibrium-instrumentation"))]
impl EquilibriumEvaluationCounters {
    fn increment(&mut self, counter: Counter) {
        match counter {
            Counter::Context => self.context_evaluations += 1,
            Counter::ResidualContext => self.residual_context_evaluations += 1,
            Counter::JacobianContext => self.jacobian_context_evaluations += 1,
            Counter::Residual => self.residual_evaluations += 1,
            Counter::Jacobian => self.jacobian_evaluations += 1,
            Counter::CandidateSnapshot => self.candidate_snapshot_evaluations += 1,
            Counter::GibbsVector => self.gibbs_vector_evaluations += 1,
            Counter::EnthalpyVector => self.enthalpy_vector_evaluations += 1,
            Counter::HeatCapacityVector => self.heat_capacity_vector_evaluations += 1,
            Counter::GibbsComponent => self.gibbs_component_calls += 1,
            Counter::EnthalpyComponent => self.enthalpy_component_calls += 1,
            Counter::HeatCapacityComponent => self.heat_capacity_component_calls += 1,
            Counter::GibbsCoefficientRefresh => self.gibbs_coefficient_refreshes += 1,
            Counter::GibbsIntervalReuse => self.gibbs_interval_reuses += 1,
            Counter::PropertyCoefficientSelection => self.property_coefficient_selections += 1,
            Counter::PropertyCalculation => self.property_calculations += 1,
            Counter::PhaseEquilibriumBuild => self.phase_equilibrium_builds += 1,
        }
    }
}

#[cfg(any(test, feature = "equilibrium-instrumentation"))]
thread_local! {
    static ACTIVE_COUNTERS: std::cell::RefCell<Option<Arc<Mutex<EquilibriumEvaluationCounters>>>> =
        const { std::cell::RefCell::new(None) };
}

#[cfg(any(test, feature = "equilibrium-instrumentation"))]
struct ActiveCapture {
    previous: Option<Arc<Mutex<EquilibriumEvaluationCounters>>>,
}

#[cfg(any(test, feature = "equilibrium-instrumentation"))]
impl Drop for ActiveCapture {
    fn drop(&mut self) {
        ACTIVE_COUNTERS.with(|active| {
            *active.borrow_mut() = self.previous.take();
        });
    }
}

/// Runs one operation with the shared equilibrium counters enabled.
///
/// This is the single reusable benchmark/test entry point. It aggregates
/// callback calls through a mutex, so an explicitly captured operation can
/// safely use callbacks that run on worker threads. The default library build
/// returns a disabled empty report and has no counter hooks; compile with
/// `--features equilibrium-instrumentation` to collect release measurements.
pub fn capture<T, F>(operation: F) -> (T, EquilibriumInstrumentationReport)
where
    F: FnOnce() -> T,
{
    #[cfg(any(test, feature = "equilibrium-instrumentation"))]
    {
        let started = std::time::Instant::now();
        let counters = Arc::new(Mutex::new(EquilibriumEvaluationCounters::default()));
        let previous =
            ACTIVE_COUNTERS.with(|active| active.borrow_mut().replace(Arc::clone(&counters)));
        let _capture = ActiveCapture { previous };
        let result = operation();
        let snapshot = counters
            .lock()
            .map(|counters| *counters)
            .unwrap_or_default();
        return (
            result,
            EquilibriumInstrumentationReport {
                enabled: true,
                elapsed: started.elapsed(),
                counters: snapshot,
            },
        );
    }

    #[cfg(not(any(test, feature = "equilibrium-instrumentation")))]
    {
        (operation(), EquilibriumInstrumentationReport::default())
    }
}

/// Captures only when the caller explicitly enabled instrumentation.
pub fn capture_if_enabled<T, F>(
    mode: EquilibriumInstrumentationMode,
    operation: F,
) -> (T, EquilibriumInstrumentationReport)
where
    F: FnOnce() -> T,
{
    match mode {
        EquilibriumInstrumentationMode::Disabled => {
            (operation(), EquilibriumInstrumentationReport::disabled())
        }
        EquilibriumInstrumentationMode::Enabled => capture(operation),
    }
}

fn record(counter: Counter) {
    #[cfg(any(test, feature = "equilibrium-instrumentation"))]
    ACTIVE_COUNTERS.with(|active| {
        if let Some(counters) = active.borrow().as_ref() {
            if let Ok(mut counters) = counters.lock() {
                counters.increment(counter);
            }
        }
    });
    #[cfg(not(any(test, feature = "equilibrium-instrumentation")))]
    let _ = counter;
}

pub(crate) fn record_context_evaluation() {
    record(Counter::Context);
}

pub(crate) fn record_residual_context_evaluation() {
    record(Counter::ResidualContext);
}

pub(crate) fn record_jacobian_context_evaluation() {
    record(Counter::JacobianContext);
}

pub(crate) fn record_residual_evaluation() {
    record(Counter::Residual);
}

pub(crate) fn record_jacobian_evaluation() {
    record(Counter::Jacobian);
}

pub(crate) fn record_candidate_snapshot_evaluation() {
    record(Counter::CandidateSnapshot);
}

pub(crate) fn record_vector_evaluation(field: &'static str) {
    match field {
        "standard_gibbs" => record(Counter::GibbsVector),
        "molar_enthalpy" => record(Counter::EnthalpyVector),
        "heat_capacity" => record(Counter::HeatCapacityVector),
        _ => {}
    }
}

pub(crate) fn record_component_call(field: &'static str) {
    match field {
        "gibbs" => record(Counter::GibbsComponent),
        "enthalpy" => record(Counter::EnthalpyComponent),
        "heat_capacity" => record(Counter::HeatCapacityComponent),
        _ => {}
    }
}

pub(crate) fn record_gibbs_coefficient_refresh() {
    record(Counter::GibbsCoefficientRefresh);
}

pub(crate) fn record_gibbs_interval_reuse() {
    record(Counter::GibbsIntervalReuse);
}

pub(crate) fn record_property_coefficient_selection() {
    record(Counter::PropertyCoefficientSelection);
}

pub(crate) fn record_property_calculation() {
    record(Counter::PropertyCalculation);
}

pub(crate) fn record_phase_equilibrium_build() {
    record(Counter::PhaseEquilibriumBuild);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn capture_starts_with_an_empty_counter_snapshot() {
        let (_, report) = capture(|| 7_u8);
        assert!(report.enabled());
        assert_eq!(report.counters(), EquilibriumEvaluationCounters::default());
    }

    #[test]
    fn capture_aggregates_named_evaluations() {
        let (_, report) = capture(|| {
            record_context_evaluation();
            record_vector_evaluation("molar_enthalpy");
            record_component_call("enthalpy");
            record_property_calculation();
        });
        let counters = report.counters();
        assert_eq!(counters.context_evaluations(), 1);
        assert_eq!(counters.enthalpy_vector_evaluations(), 1);
        assert_eq!(counters.enthalpy_component_calls(), 1);
        assert_eq!(counters.property_calculations(), 1);
    }

    #[test]
    fn legacy_stage_capture_preserves_the_internal_breakdown() {
        let (snapshot, report) = capture(|| {
            let session =
                begin_legacy_solver_diagnostics().expect("capture must enable Legacy scope");
            record_legacy_stage(
                LegacySolverStage::MatrixOperations,
                Duration::from_micros(3),
            );
            record_legacy_stage(LegacySolverStage::Norm, Duration::from_micros(5));
            record_legacy_stage(LegacySolverStage::TrialVector, Duration::from_micros(7));
            record_legacy_stage(LegacySolverStage::Feasibility, Duration::from_micros(11));
            record_legacy_stage(LegacySolverStage::Globalization, Duration::from_micros(13));
            finish_legacy_solver_diagnostics(session)
        });

        assert!(report.enabled());
        assert_eq!(snapshot.matrix_operations_micros, 3);
        assert_eq!(snapshot.norm_micros, 5);
        assert_eq!(snapshot.trial_vector_micros, 7);
        assert_eq!(snapshot.feasibility_micros, 11);
        assert_eq!(snapshot.globalization_micros, 13);
    }

    #[test]
    fn legacy_callback_capture_counts_only_adjacent_matching_iterates() {
        let (snapshot, report) = capture(|| {
            let session =
                begin_legacy_solver_diagnostics().expect("capture must enable Legacy scope");
            record_legacy_callback(LegacyCallbackKind::Residual, &[1.0, 2.0]);
            record_legacy_callback(LegacyCallbackKind::Jacobian, &[1.0, 2.0]);
            record_legacy_callback(LegacyCallbackKind::Residual, &[1.0, 2.1]);
            record_legacy_callback(LegacyCallbackKind::Jacobian, &[1.0, 2.0]);
            finish_legacy_solver_diagnostics(session)
        });

        assert!(report.enabled());
        assert_eq!(snapshot.same_iterate_callback_pairs, 1);
    }
}
