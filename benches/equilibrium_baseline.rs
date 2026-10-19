//! Criterion baseline for the production equilibrium facade.
//!
//! Run with:
//!
//! ```text
//! cargo bench --bench equilibrium_baseline --features equilibrium-instrumentation
//! ```
//!
//! Add `equilibrium-benchmark-prettyprint` for a stage-attribution table with
//! wall time and percentages. The pretty report is hierarchical: macro
//! partition, facade/P,H/inner detail, nonlinear coverage, every backend
//! attempt, and accepted callback detail use separate tables and denominators.
//!
//! Criterion owns the statistically measured wall time. The one-time table
//! printed before Criterion starts comes from the same accepted facade
//! results, using the library's opt-in stage timing and evaluation counters.
//! It is intentionally a report, not a benchmark assertion: hardware,
//! compiler, and repository changes can legitimately move the numbers.

use std::time::Duration;

use criterion::{BenchmarkId, Criterion, black_box, criterion_group, criterion_main};

#[cfg(feature = "equilibrium-benchmark-prettyprint")]
use tabled::{Table, Tabled};

use KiThe::Thermodynamics::ChemEquilibrium::prelude::*;

mod support;

use support::{
    PRESSURE_PA, REFERENCE_PRESSURE_PA, REFERENCE_TEMPERATURE_K, RP1311_40, STANJAN_15,
    TP1907_GAS_17, component, tp1907_element_seed,
};

const RANGE_TEMPERATURES: [f64; 3] = [2_400.0, 2_500.0, 2_600.0];
const BASELINE_REPETITIONS: usize = 10;
const CEA_RP1311_40: [&str; 40] = RP1311_40;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Universe {
    SmallH2O,
    Stanjan15,
    Tp1907,
    Rp1311_40,
}

impl Universe {
    fn label(self) -> &'static str {
        match self {
            Self::SmallH2O => "small3",
            Self::Stanjan15 => "stanjan15",
            Self::Tp1907 => "tp1907_18",
            Self::Rp1311_40 => "rp1311_40",
        }
    }
}

#[derive(Clone, Copy, Debug)]
enum ScenarioTask {
    PtMolecular,
    PtElementInventory,
    PhMolecular,
    PhElementInventory,
    PtElementInventoryRange,
    PhElementInventoryRange,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum BackendRoute {
    Production,
    LegacyLm,
    LegacyNr,
    LegacyTr,
}

impl BackendRoute {
    fn label(self) -> &'static str {
        match self {
            Self::Production => "production",
            Self::LegacyLm => "legacy_lm",
            Self::LegacyNr => "legacy_nr",
            Self::LegacyTr => "legacy_tr",
        }
    }
}

#[derive(Clone, Copy, Debug)]
struct Scenario {
    universe: Universe,
    task: ScenarioTask,
    route: BackendRoute,
}

impl Scenario {
    fn label(self) -> String {
        let task = match self.task {
            ScenarioTask::PtMolecular => "pt_molecular_point",
            ScenarioTask::PtElementInventory => "pt_element_inventory_point",
            ScenarioTask::PhMolecular => "ph_molecular_point",
            ScenarioTask::PhElementInventory => "ph_element_inventory_point",
            ScenarioTask::PtElementInventoryRange => "pt_element_inventory_range",
            ScenarioTask::PhElementInventoryRange => "ph_element_inventory_range",
        };
        let base = if matches!(self.universe, Universe::SmallH2O) {
            task.to_string()
        } else {
            format!("{}_{}", self.universe.label(), task)
        };
        if self.route == BackendRoute::Production {
            base
        } else {
            format!("{}_{}", base, self.route.label())
        }
    }
}

/// Numeric values kept deliberately flat so the report remains easy to paste
/// into STORY_TESTS.md and compare with a later optimization run.
#[derive(Clone, Debug, Default)]
struct Measurement {
    facade_us: f64,
    inner_total_us: f64,
    workflow_total_us: f64,
    facade_resolution_us: f64,
    request_preparation_us: f64,
    route_selection_us: f64,
    normalization_us: f64,
    phase_control_us: f64,
    facade_thermochemistry_us: f64,
    ph_formulation_preparation_us: f64,
    recovery_us: f64,
    result_assembly_us: f64,
    scalar_us: f64,
    enthalpy_us: f64,
    ph_trial_total_us: f64,
    ph_trial_inner_execution_us: f64,
    ph_trial_preparation_us: f64,
    ph_trial_template_preparation_us: f64,
    ph_trial_postprocessing_us: f64,
    ph_bracket_orchestration_us: f64,
    recovery_physical_failure_us: f64,
    recovery_normalization_us: f64,
    recovery_normalized_discovery_us: f64,
    recovery_physical_retry_us: f64,
    recovery_reconstruction_us: f64,
    ph_monolithic_preparation_us: f64,
    ph_monolithic_rst_preparation_us: f64,
    ph_monolithic_backend_execution_us: f64,
    ph_monolithic_validation_us: f64,
    ph_monolithic_execution_us: f64,
    ph_monolithic_publication_us: f64,
    repository_lookup_us: f64,
    thermochemistry_us: f64,
    numeric_closure_us: f64,
    symbolic_us: f64,
    equation_us: f64,
    numerical_problem_us: f64,
    active_set_formulation_build_us: f64,
    rst_preparation_us: f64,
    projection_build_us: f64,
    nonlinear_us: f64,
    postprocessing_us: f64,
    validation_us: f64,
    residual_evaluations: f64,
    jacobian_evaluations: f64,
    gibbs_vectors: f64,
    enthalpy_vectors: f64,
    heat_capacity_vectors: f64,
    formulation_builds: f64,
    backend_attempts: f64,
    nonlinear_iterations: f64,
    backend_us: f64,
    backend_callback_us: f64,
    backend_overhead_us: f64,
    solver_iterations: f64,
    solver_residual_calls: f64,
    solver_jacobian_calls: f64,
    solver_linear_system_solves: f64,
    solver_residual_us: f64,
    solver_jacobian_us: f64,
    solver_linear_system_us: f64,
    solver_overhead_us: f64,
    solver_linear_timing_available: bool,
    accepted_backend: &'static str,
    solver_name: &'static str,
    ph_route: &'static str,
    ph_route_decisions: f64,
    #[cfg(feature = "equilibrium-benchmark-prettyprint")]
    backend_attempt_rows: Vec<BackendAttemptMeasurement>,
}

#[cfg(feature = "equilibrium-benchmark-prettyprint")]
#[derive(Clone, Debug, Default)]
struct BackendAttemptMeasurement {
    scope: String,
    attempt_index: usize,
    backend: String,
    outcome: String,
    samples: f64,
    elapsed_us: f64,
    iterations: f64,
    residual_calls: f64,
    jacobian_calls: f64,
    linear_solves: f64,
    residual_us: f64,
    jacobian_us: f64,
    linear_us: f64,
    overhead_us: f64,
    metrics_available: bool,
    wall_time_available: bool,
    evaluation_timing_available: bool,
    linear_timing_available: bool,
}

impl Measurement {
    fn add_assign(&mut self, other: Self) {
        self.facade_us += other.facade_us;
        self.inner_total_us += other.inner_total_us;
        self.workflow_total_us += other.workflow_total_us;
        self.facade_resolution_us += other.facade_resolution_us;
        self.request_preparation_us += other.request_preparation_us;
        self.route_selection_us += other.route_selection_us;
        self.normalization_us += other.normalization_us;
        self.phase_control_us += other.phase_control_us;
        self.facade_thermochemistry_us += other.facade_thermochemistry_us;
        self.ph_formulation_preparation_us += other.ph_formulation_preparation_us;
        self.recovery_us += other.recovery_us;
        self.result_assembly_us += other.result_assembly_us;
        self.scalar_us += other.scalar_us;
        self.enthalpy_us += other.enthalpy_us;
        self.ph_trial_total_us += other.ph_trial_total_us;
        self.ph_trial_inner_execution_us += other.ph_trial_inner_execution_us;
        self.ph_trial_preparation_us += other.ph_trial_preparation_us;
        self.ph_trial_template_preparation_us += other.ph_trial_template_preparation_us;
        self.ph_trial_postprocessing_us += other.ph_trial_postprocessing_us;
        self.ph_bracket_orchestration_us += other.ph_bracket_orchestration_us;
        self.recovery_physical_failure_us += other.recovery_physical_failure_us;
        self.recovery_normalization_us += other.recovery_normalization_us;
        self.recovery_normalized_discovery_us += other.recovery_normalized_discovery_us;
        self.recovery_physical_retry_us += other.recovery_physical_retry_us;
        self.recovery_reconstruction_us += other.recovery_reconstruction_us;
        self.ph_monolithic_preparation_us += other.ph_monolithic_preparation_us;
        self.ph_monolithic_rst_preparation_us += other.ph_monolithic_rst_preparation_us;
        self.ph_monolithic_backend_execution_us += other.ph_monolithic_backend_execution_us;
        self.ph_monolithic_validation_us += other.ph_monolithic_validation_us;
        self.ph_monolithic_execution_us += other.ph_monolithic_execution_us;
        self.ph_monolithic_publication_us += other.ph_monolithic_publication_us;
        self.repository_lookup_us += other.repository_lookup_us;
        self.thermochemistry_us += other.thermochemistry_us;
        self.numeric_closure_us += other.numeric_closure_us;
        self.symbolic_us += other.symbolic_us;
        self.equation_us += other.equation_us;
        self.numerical_problem_us += other.numerical_problem_us;
        self.active_set_formulation_build_us += other.active_set_formulation_build_us;
        self.rst_preparation_us += other.rst_preparation_us;
        self.projection_build_us += other.projection_build_us;
        self.nonlinear_us += other.nonlinear_us;
        self.postprocessing_us += other.postprocessing_us;
        self.validation_us += other.validation_us;
        self.residual_evaluations += other.residual_evaluations;
        self.jacobian_evaluations += other.jacobian_evaluations;
        self.gibbs_vectors += other.gibbs_vectors;
        self.enthalpy_vectors += other.enthalpy_vectors;
        self.heat_capacity_vectors += other.heat_capacity_vectors;
        self.formulation_builds += other.formulation_builds;
        self.backend_attempts += other.backend_attempts;
        self.nonlinear_iterations += other.nonlinear_iterations;
        self.backend_us += other.backend_us;
        self.backend_callback_us += other.backend_callback_us;
        self.backend_overhead_us += other.backend_overhead_us;
        self.solver_iterations += other.solver_iterations;
        self.solver_residual_calls += other.solver_residual_calls;
        self.solver_jacobian_calls += other.solver_jacobian_calls;
        self.solver_linear_system_solves += other.solver_linear_system_solves;
        self.solver_residual_us += other.solver_residual_us;
        self.solver_jacobian_us += other.solver_jacobian_us;
        self.solver_linear_system_us += other.solver_linear_system_us;
        self.solver_overhead_us += other.solver_overhead_us;
        self.solver_linear_timing_available |= other.solver_linear_timing_available;
        if !other.accepted_backend.is_empty() {
            self.accepted_backend = other.accepted_backend;
        }
        if !other.solver_name.is_empty() {
            self.solver_name = other.solver_name;
        }
        if !other.ph_route.is_empty() {
            self.ph_route = other.ph_route;
        }
        self.ph_route_decisions += other.ph_route_decisions;
        #[cfg(feature = "equilibrium-benchmark-prettyprint")]
        for incoming in other.backend_attempt_rows {
            if let Some(existing) = self.backend_attempt_rows.iter_mut().find(|existing| {
                existing.scope == incoming.scope
                    && existing.attempt_index == incoming.attempt_index
                    && existing.backend == incoming.backend
                    && existing.outcome == incoming.outcome
            }) {
                existing.samples += incoming.samples;
                existing.elapsed_us += incoming.elapsed_us;
                existing.iterations += incoming.iterations;
                existing.residual_calls += incoming.residual_calls;
                existing.jacobian_calls += incoming.jacobian_calls;
                existing.linear_solves += incoming.linear_solves;
                existing.residual_us += incoming.residual_us;
                existing.jacobian_us += incoming.jacobian_us;
                existing.linear_us += incoming.linear_us;
                existing.overhead_us += incoming.overhead_us;
                existing.metrics_available |= incoming.metrics_available;
                existing.wall_time_available |= incoming.wall_time_available;
                existing.evaluation_timing_available |= incoming.evaluation_timing_available;
                existing.linear_timing_available |= incoming.linear_timing_available;
            } else {
                self.backend_attempt_rows.push(incoming);
            }
        }
    }

    fn divided_by(mut self, divisor: f64) -> Self {
        self.facade_us /= divisor;
        self.inner_total_us /= divisor;
        self.workflow_total_us /= divisor;
        self.facade_resolution_us /= divisor;
        self.request_preparation_us /= divisor;
        self.route_selection_us /= divisor;
        self.normalization_us /= divisor;
        self.phase_control_us /= divisor;
        self.facade_thermochemistry_us /= divisor;
        self.ph_formulation_preparation_us /= divisor;
        self.recovery_us /= divisor;
        self.result_assembly_us /= divisor;
        self.scalar_us /= divisor;
        self.enthalpy_us /= divisor;
        self.ph_trial_total_us /= divisor;
        self.ph_trial_inner_execution_us /= divisor;
        self.ph_trial_preparation_us /= divisor;
        self.ph_trial_template_preparation_us /= divisor;
        self.ph_trial_postprocessing_us /= divisor;
        self.ph_bracket_orchestration_us /= divisor;
        self.recovery_physical_failure_us /= divisor;
        self.recovery_normalization_us /= divisor;
        self.recovery_normalized_discovery_us /= divisor;
        self.recovery_physical_retry_us /= divisor;
        self.recovery_reconstruction_us /= divisor;
        self.ph_monolithic_preparation_us /= divisor;
        self.ph_monolithic_rst_preparation_us /= divisor;
        self.ph_monolithic_backend_execution_us /= divisor;
        self.ph_monolithic_validation_us /= divisor;
        self.ph_monolithic_execution_us /= divisor;
        self.ph_monolithic_publication_us /= divisor;
        self.repository_lookup_us /= divisor;
        self.thermochemistry_us /= divisor;
        self.numeric_closure_us /= divisor;
        self.symbolic_us /= divisor;
        self.equation_us /= divisor;
        self.numerical_problem_us /= divisor;
        self.active_set_formulation_build_us /= divisor;
        self.rst_preparation_us /= divisor;
        self.projection_build_us /= divisor;
        self.nonlinear_us /= divisor;
        self.postprocessing_us /= divisor;
        self.validation_us /= divisor;
        self.residual_evaluations /= divisor;
        self.jacobian_evaluations /= divisor;
        self.gibbs_vectors /= divisor;
        self.enthalpy_vectors /= divisor;
        self.heat_capacity_vectors /= divisor;
        self.formulation_builds /= divisor;
        self.backend_attempts /= divisor;
        self.nonlinear_iterations /= divisor;
        self.backend_us /= divisor;
        self.backend_callback_us /= divisor;
        self.backend_overhead_us /= divisor;
        self.solver_iterations /= divisor;
        self.solver_residual_calls /= divisor;
        self.solver_jacobian_calls /= divisor;
        self.solver_linear_system_solves /= divisor;
        self.solver_residual_us /= divisor;
        self.solver_jacobian_us /= divisor;
        self.solver_linear_system_us /= divisor;
        self.solver_overhead_us /= divisor;
        self.ph_route_decisions /= divisor;
        #[cfg(feature = "equilibrium-benchmark-prettyprint")]
        for attempt in &mut self.backend_attempt_rows {
            attempt.samples /= divisor;
            attempt.elapsed_us /= divisor;
            attempt.iterations /= divisor;
            attempt.residual_calls /= divisor;
            attempt.jacobian_calls /= divisor;
            attempt.linear_solves /= divisor;
            attempt.residual_us /= divisor;
            attempt.jacobian_us /= divisor;
            attempt.linear_us /= divisor;
            attempt.overhead_us /= divisor;
        }
        self
    }
}

fn percentage(value_us: f64, wall_us: f64) -> f64 {
    if wall_us <= f64::EPSILON {
        0.0
    } else {
        100.0 * value_us.max(0.0) / wall_us
    }
}

#[cfg(feature = "equilibrium-benchmark-prettyprint")]
fn timing_profile(measurement: &Measurement) -> &'static str {
    if measurement.ph_route.is_empty() || measurement.ph_route == "pt" {
        "P,T accepted"
    } else if measurement.recovery_us > 0.0 {
        "P,H recovery"
    } else if measurement.ph_route.starts_with("monolithic") {
        "P,H monolithic"
    } else {
        "P,H nested"
    }
}

#[cfg(feature = "equilibrium-benchmark-prettyprint")]
#[derive(Tabled)]
struct PrettySolverMeasurement {
    #[tabled(rename = "case")]
    case_name: String,
    #[tabled(rename = "backend")]
    backend: String,
    #[tabled(rename = "route")]
    route: String,
    #[tabled(rename = "profile")]
    profile: &'static str,
    #[tabled(rename = "backend wall")]
    backend_wall: String,
    #[tabled(rename = "iterations")]
    iterations: String,
    #[tabled(rename = "res calls")]
    residual_calls: String,
    #[tabled(rename = "res us")]
    residual_time: String,
    #[tabled(rename = "res %")]
    residual_percent: String,
    #[tabled(rename = "Jac calls")]
    jacobian_calls: String,
    #[tabled(rename = "Jac us")]
    jacobian_time: String,
    #[tabled(rename = "Jac %")]
    jacobian_percent: String,
    #[tabled(rename = "linear solves")]
    linear_system_solves: String,
    #[tabled(rename = "linear us")]
    linear_system_time: String,
    #[tabled(rename = "linear %")]
    linear_system_percent: String,
    #[tabled(rename = "overhead us")]
    solver_overhead_time: String,
    #[tabled(rename = "overhead %")]
    solver_overhead_percent: String,
}

#[cfg(feature = "equilibrium-benchmark-prettyprint")]
impl PrettySolverMeasurement {
    fn from_measurement(case_name: String, measurement: &Measurement) -> Self {
        let wall = measurement.backend_us.max(0.0);
        let linear_time = if measurement.solver_linear_timing_available {
            format!("{:.1}", measurement.solver_linear_system_us)
        } else {
            "n/a".to_string()
        };
        let linear_percent = if measurement.solver_linear_timing_available {
            format!(
                "{:.1}%",
                percentage(measurement.solver_linear_system_us, wall)
            )
        } else {
            "n/a".to_string()
        };
        Self {
            case_name,
            backend: measurement.accepted_backend.to_string(),
            route: if measurement.ph_route.is_empty() {
                "pt".to_string()
            } else {
                measurement.ph_route.to_string()
            },
            profile: timing_profile(measurement),
            backend_wall: format!("{:.2} ms", wall / 1_000.0),
            iterations: format!("{:.1}", measurement.solver_iterations),
            residual_calls: format!("{:.1}", measurement.solver_residual_calls),
            residual_time: format!("{:.1}", measurement.solver_residual_us),
            residual_percent: format!("{:.1}%", percentage(measurement.solver_residual_us, wall)),
            jacobian_calls: format!("{:.1}", measurement.solver_jacobian_calls),
            jacobian_time: format!("{:.1}", measurement.solver_jacobian_us),
            jacobian_percent: format!("{:.1}%", percentage(measurement.solver_jacobian_us, wall)),
            linear_system_solves: format!("{:.1}", measurement.solver_linear_system_solves),
            linear_system_time: linear_time,
            linear_system_percent: linear_percent,
            solver_overhead_time: format!("{:.1}", measurement.solver_overhead_us),
            solver_overhead_percent: format!(
                "{:.1}%",
                percentage(measurement.solver_overhead_us, wall)
            ),
        }
    }
}

#[cfg(feature = "equilibrium-benchmark-prettyprint")]
#[derive(Tabled)]
struct PrettyMacroMeasurement {
    #[tabled(rename = "case")]
    case_name: String,
    #[tabled(rename = "backend")]
    backend: String,
    #[tabled(rename = "route")]
    route: String,
    #[tabled(rename = "profile")]
    profile: &'static str,
    #[tabled(rename = "wall")]
    wall: String,
    #[tabled(rename = "facade outside workflow")]
    facade_outside_workflow: String,
    #[tabled(rename = "solve transaction")]
    workflow: String,
    #[tabled(rename = "partition")]
    partition: String,
}

#[cfg(feature = "equilibrium-benchmark-prettyprint")]
#[derive(Tabled)]
struct PrettyLegacyPhEnvelopeRow {
    #[tabled(rename = "case")]
    case_name: String,
    #[tabled(rename = "backend")]
    backend: String,
    #[tabled(rename = "P,H route")]
    route: String,
    #[tabled(rename = "profile")]
    profile: &'static str,
    #[tabled(rename = "facade total")]
    facade_total: String,
    #[tabled(rename = "P,H workflow")]
    ph_workflow: String,
    #[tabled(rename = "trial total")]
    trial_total: String,
    #[tabled(rename = "nonlinear solve")]
    nonlinear_solve: String,
    #[tabled(rename = "Legacy backend wall")]
    legacy_backend_wall: String,
}

#[cfg(feature = "equilibrium-benchmark-prettyprint")]
impl PrettyLegacyPhEnvelopeRow {
    fn from_measurement(case_name: String, measurement: &Measurement) -> Option<Self> {
        if measurement.ph_route.is_empty() || measurement.solver_name == "" {
            return None;
        }
        let facade_us = measurement.facade_us.max(0.0);
        let workflow_us = measurement.workflow_total_us.max(0.0);
        let nested = measurement.ph_route == "nested";
        let nonlinear_us = if nested {
            measurement.nonlinear_us.max(0.0)
        } else {
            measurement
                .ph_monolithic_backend_execution_us
                .max(measurement.nonlinear_us)
                .max(0.0)
        };
        let backend_us = measurement.backend_us.max(0.0);
        Some(Self {
            case_name,
            backend: measurement.solver_name.to_string(),
            route: measurement.ph_route.to_string(),
            profile: timing_profile(measurement),
            facade_total: format!("{:.2} ms", facade_us / 1_000.0),
            ph_workflow: scoped_metric(workflow_us, facade_us),
            trial_total: if nested {
                scoped_metric(measurement.ph_trial_total_us.max(0.0), workflow_us)
            } else {
                "n/a (monolithic)".to_string()
            },
            nonlinear_solve: scoped_metric(
                nonlinear_us,
                if nested {
                    measurement.ph_trial_total_us.max(0.0)
                } else {
                    measurement.ph_monolithic_execution_us.max(workflow_us)
                },
            ),
            // Accepted backend metrics and P,H timing are intentionally kept
            // as separate observations: nested routes can contain several
            // trial/attempt scopes, so forcing one ratio can be misleading.
            legacy_backend_wall: format!("{:.2} ms", backend_us / 1_000.0),
        })
    }
}

#[cfg(feature = "equilibrium-benchmark-prettyprint")]
fn scoped_metric(value_us: f64, scope_us: f64) -> String {
    format!(
        "{:.2} ms / {:.1}%",
        value_us / 1_000.0,
        percentage(value_us, scope_us)
    )
}

#[cfg(feature = "equilibrium-benchmark-prettyprint")]
impl PrettyMacroMeasurement {
    fn from_measurement(case_name: String, measurement: &Measurement) -> Self {
        let wall = measurement.facade_us.max(0.0);
        let outside = (measurement.facade_us - measurement.workflow_total_us).max(0.0);
        let outside_percent = percentage(outside, wall);
        let workflow_percent = percentage(measurement.workflow_total_us, wall);
        Self {
            case_name,
            backend: measurement.solver_name.to_string(),
            route: if measurement.ph_route.is_empty() {
                "pt".to_string()
            } else {
                measurement.ph_route.to_string()
            },
            profile: timing_profile(measurement),
            wall: format!("{:.2} ms", wall / 1_000.0),
            facade_outside_workflow: format!(
                "{:.2} ms / {:.1}%",
                outside / 1_000.0,
                outside_percent
            ),
            workflow: format!(
                "{:.2} ms / {:.1}%",
                measurement.workflow_total_us / 1_000.0,
                workflow_percent
            ),
            partition: format!("{:.1}%", outside_percent + workflow_percent),
        }
    }
}

#[cfg(feature = "equilibrium-benchmark-prettyprint")]
#[derive(Tabled)]
struct PrettyTimingRow {
    #[tabled(rename = "case")]
    case_name: String,
    #[tabled(rename = "scope")]
    scope: String,
    #[tabled(rename = "profile")]
    profile: &'static str,
    #[tabled(rename = "stage")]
    stage: String,
    #[tabled(rename = "time")]
    duration: String,
    #[tabled(rename = "% of scope wall")]
    percent: String,
    #[tabled(rename = "semantics")]
    semantics: String,
}

#[cfg(feature = "equilibrium-benchmark-prettyprint")]
fn timing_row(
    case_name: &str,
    scope: &str,
    profile: &'static str,
    stage: &str,
    duration_us: f64,
    scope_us: f64,
    semantics: &str,
) -> PrettyTimingRow {
    PrettyTimingRow {
        case_name: case_name.to_string(),
        scope: scope.to_string(),
        profile,
        stage: stage.to_string(),
        duration: format!("{:.2} ms", duration_us.max(0.0) / 1_000.0),
        percent: format!("{:.1}%", percentage(duration_us, scope_us)),
        semantics: semantics.to_string(),
    }
}

#[cfg(feature = "equilibrium-benchmark-prettyprint")]
fn facade_timing_rows(case_name: &str, measurement: &Measurement) -> Vec<PrettyTimingRow> {
    let scope = measurement.facade_us;
    [
        ("resolution", measurement.facade_resolution_us),
        ("request preparation", measurement.request_preparation_us),
        ("route selection", measurement.route_selection_us),
        ("normalization", measurement.normalization_us),
        ("phase control", measurement.phase_control_us),
        ("thermochemistry", measurement.facade_thermochemistry_us),
        ("result assembly", measurement.result_assembly_us),
    ]
    .into_iter()
    .map(|(stage, duration)| {
        timing_row(
            case_name,
            "facade",
            timing_profile(measurement),
            stage,
            duration,
            scope,
            "measured facade stage",
        )
    })
    .collect()
}

#[cfg(feature = "equilibrium-benchmark-prettyprint")]
fn inner_timing_rows(case_name: &str, measurement: &Measurement) -> Vec<PrettyTimingRow> {
    if measurement.ph_route.starts_with("monolithic") {
        return Vec::new();
    }
    let scope_name = if measurement.recovery_us > 0.0 {
        "P,H recovery transaction"
    } else {
        "inner fixed-P,T"
    };
    // A recovered P,H point replaces the accepted inner timing with the
    // complete transaction wall time. The accepted normalized solve remains
    // available as a child observation, but it is not the recovery budget.
    let scope = if measurement.recovery_us > 0.0 {
        measurement.recovery_us.max(measurement.inner_total_us)
    } else {
        measurement.inner_total_us
    };
    let mut rows = vec![timing_row(
        case_name,
        scope_name,
        timing_profile(measurement),
        "outside nonlinear envelope",
        (scope - measurement.nonlinear_us).max(0.0),
        scope,
        "exclusive arithmetic complement of nonlinear solve",
    )];
    if measurement.recovery_us > 0.0 {
        let recovery_children = measurement.recovery_normalization_us.max(0.0)
            + measurement.recovery_physical_failure_us.max(0.0)
            + measurement.recovery_normalized_discovery_us.max(0.0)
            + measurement.recovery_physical_retry_us.max(0.0)
            + measurement.recovery_reconstruction_us.max(0.0);
        let recovery_remainder = (measurement.recovery_us - recovery_children).max(0.0);
        rows.extend(
            [
                (
                    "recovery physical failure",
                    measurement.recovery_physical_failure_us,
                ),
                (
                    "recovery normalization",
                    measurement.recovery_normalization_us,
                ),
                (
                    "recovery normalized discovery",
                    measurement.recovery_normalized_discovery_us,
                ),
                (
                    "recovery physical retry",
                    measurement.recovery_physical_retry_us,
                ),
                (
                    "recovery reconstruction",
                    measurement.recovery_reconstruction_us,
                ),
                ("recovery remainder", recovery_remainder),
            ]
            .into_iter()
            .map(|(stage, duration)| {
                timing_row(
                    case_name,
                    scope_name,
                    timing_profile(measurement),
                    stage,
                    duration,
                    measurement.recovery_us,
                    "exclusive recovery child or residual; recovery scope includes physical failure",
                )
            }),
        );
    }
    rows.extend(
        [
            ("thermochemistry", measurement.thermochemistry_us),
            ("numeric closure", measurement.numeric_closure_us),
            ("normalization", measurement.normalization_us),
            ("recovery", measurement.recovery_us),
            ("symbolic construction", measurement.symbolic_us),
            ("equation construction", measurement.equation_us),
            (
                "numerical problem preparation",
                measurement.numerical_problem_us,
            ),
            (
                "active-set formulation build",
                measurement.active_set_formulation_build_us,
            ),
            ("RST preparation", measurement.rst_preparation_us),
            ("projection build", measurement.projection_build_us),
            ("nonlinear solve", measurement.nonlinear_us),
            ("phase-control envelope", measurement.phase_control_us),
            ("validation", measurement.validation_us),
            ("postprocessing", measurement.postprocessing_us),
            ("repository lookup", measurement.repository_lookup_us),
            ("result assembly", measurement.result_assembly_us),
        ]
        .into_iter()
        .map(|(stage, duration)| {
            timing_row(
                case_name,
                scope_name,
                timing_profile(measurement),
                stage,
                duration,
                scope,
                "inclusive; stages are not additive",
            )
        }),
    );
    rows
}

#[cfg(feature = "equilibrium-benchmark-prettyprint")]
fn ph_timing_rows(case_name: &str, measurement: &Measurement) -> Vec<PrettyTimingRow> {
    let scope = measurement.workflow_total_us;
    let mut rows = Vec::new();
    if measurement.recovery_us > 0.0 {
        rows.push(timing_row(
            case_name,
            "P,H workflow",
            timing_profile(measurement),
            "recovery transaction",
            measurement.recovery_us,
            scope,
            "inclusive complete fallback transaction; child timings overlap",
        ));
    }
    rows.extend(
        [
            ("P,H scalar orchestration", measurement.scalar_us),
            ("enthalpy evaluation", measurement.enthalpy_us),
            ("trial total", measurement.ph_trial_total_us),
            (
                "trial inner execution",
                measurement.ph_trial_inner_execution_us,
            ),
            ("trial preparation", measurement.ph_trial_preparation_us),
            (
                "trial template preparation",
                measurement.ph_trial_template_preparation_us,
            ),
            (
                "trial postprocessing",
                measurement.ph_trial_postprocessing_us,
            ),
            (
                "bracket orchestration",
                measurement.ph_bracket_orchestration_us,
            ),
            (
                "monolithic preparation",
                measurement.ph_monolithic_preparation_us,
            ),
            (
                "monolithic RST preparation",
                measurement.ph_monolithic_rst_preparation_us,
            ),
            (
                "monolithic backend execution",
                measurement.ph_monolithic_backend_execution_us,
            ),
            (
                "monolithic validation",
                measurement.ph_monolithic_validation_us,
            ),
            (
                "monolithic execution",
                measurement.ph_monolithic_execution_us,
            ),
            (
                "monolithic publication",
                measurement.ph_monolithic_publication_us,
            ),
        ]
        .into_iter()
        .filter(|(_, duration)| *duration > 0.0)
        .map(|(stage, duration)| {
            timing_row(
                case_name,
                "P,H workflow",
                timing_profile(measurement),
                stage,
                duration,
                scope,
                "inclusive; parent/child timings overlap",
            )
        }),
    );
    rows
}

#[cfg(feature = "equilibrium-benchmark-prettyprint")]
#[derive(Tabled)]
struct PrettyNonlinearRow {
    #[tabled(rename = "case")]
    case_name: String,
    #[tabled(rename = "parent scope")]
    parent_scope: String,
    #[tabled(rename = "profile")]
    profile: &'static str,
    #[tabled(rename = "parent wall")]
    parent_wall: String,
    #[tabled(rename = "nonlinear envelope")]
    nonlinear_envelope: String,
    #[tabled(rename = "attempts")]
    attempts: String,
    #[tabled(rename = "all attempt wall")]
    all_attempt_wall: String,
    #[tabled(rename = "known preparation")]
    known_preparation: String,
    #[tabled(rename = "runner/orchestration remainder")]
    runner_orchestration_remainder: String,
    #[tabled(rename = "backend coverage")]
    backend_coverage: String,
    #[tabled(rename = "no metrics")]
    no_metrics: String,
    #[tabled(rename = "unknown wall")]
    unknown_wall: String,
}

#[cfg(feature = "equilibrium-benchmark-prettyprint")]
fn nonlinear_rows(case_name: &str, measurement: &Measurement) -> PrettyNonlinearRow {
    let attempt_wall = measurement
        .backend_attempt_rows
        .iter()
        .map(|attempt| attempt.elapsed_us)
        .sum::<f64>();
    let (parent_scope, parent_wall_us, envelope_us) =
        if measurement.ph_route.starts_with("monolithic") {
            (
                "monolithic backend execution",
                measurement.ph_monolithic_backend_execution_us.max(0.0),
                measurement.ph_monolithic_backend_execution_us.max(0.0),
            )
        } else if measurement.recovery_us > 0.0 {
            (
                "inner recovery transaction",
                measurement.inner_total_us.max(0.0),
                measurement.recovery_us.max(0.0),
            )
        } else {
            (
                "inner fixed-P,T",
                measurement.inner_total_us.max(0.0),
                measurement.nonlinear_us.max(0.0),
            )
        };
    let known_preparation_us = if parent_scope == "monolithic backend execution" {
        0.0
    } else {
        measurement.active_set_formulation_build_us
            + measurement.rst_preparation_us
            + measurement.projection_build_us
    };
    let runner_orchestration_remainder_us =
        (envelope_us - attempt_wall - known_preparation_us).max(0.0);
    let attempt_count = measurement.backend_attempt_rows.len();
    let unreported_count = measurement
        .backend_attempt_rows
        .iter()
        .filter(|attempt| !attempt.metrics_available)
        .count();
    let unknown_wall_count = measurement
        .backend_attempt_rows
        .iter()
        .filter(|attempt| !attempt.wall_time_available && attempt.outcome != "skipped")
        .count();
    PrettyNonlinearRow {
        case_name: case_name.to_string(),
        parent_scope: parent_scope.to_string(),
        profile: timing_profile(measurement),
        parent_wall: format!("{:.2} ms", parent_wall_us / 1_000.0),
        nonlinear_envelope: format!(
            "{:.2} ms / {:.1}%",
            envelope_us / 1_000.0,
            percentage(envelope_us, parent_wall_us)
        ),
        attempts: attempt_count.to_string(),
        all_attempt_wall: format!(
            "{:.2} ms / {:.1}%",
            attempt_wall / 1_000.0,
            percentage(attempt_wall, parent_wall_us)
        ),
        known_preparation: format!(
            "{:.2} ms / {:.1}%",
            known_preparation_us / 1_000.0,
            percentage(known_preparation_us, parent_wall_us)
        ),
        runner_orchestration_remainder: format!(
            "{:.2} ms / {:.1}%",
            runner_orchestration_remainder_us / 1_000.0,
            percentage(runner_orchestration_remainder_us, parent_wall_us)
        ),
        backend_coverage: format!("{:.1}%", percentage(attempt_wall, envelope_us)),
        no_metrics: if unreported_count == 0 {
            "0".to_string()
        } else {
            format!("{unreported_count}/{attempt_count}")
        },
        unknown_wall: if unknown_wall_count == 0 {
            "0".to_string()
        } else {
            format!("{unknown_wall_count}/{attempt_count}")
        },
    }
}

#[cfg(feature = "equilibrium-benchmark-prettyprint")]
#[derive(Tabled)]
struct PrettyAttemptRow {
    #[tabled(rename = "case")]
    case_name: String,
    #[tabled(rename = "scope")]
    scope: String,
    #[tabled(rename = "profile")]
    profile: &'static str,
    #[tabled(rename = "attempt")]
    attempt: String,
    #[tabled(rename = "backend")]
    backend: String,
    #[tabled(rename = "outcome")]
    outcome: String,
    #[tabled(rename = "metrics")]
    metrics: String,
    #[tabled(rename = "samples")]
    samples: String,
    #[tabled(rename = "wall")]
    wall: String,
    #[tabled(rename = "iterations")]
    iterations: String,
    #[tabled(rename = "residual")]
    residual: String,
    #[tabled(rename = "Jacobian")]
    jacobian: String,
    #[tabled(rename = "linear")]
    linear: String,
    #[tabled(rename = "solver overhead")]
    overhead: String,
}

#[cfg(feature = "equilibrium-benchmark-prettyprint")]
fn attempt_metric(time_us: f64, wall_us: f64) -> String {
    format!(
        "{:.2} ms / {:.1}%",
        time_us / 1_000.0,
        percentage(time_us, wall_us)
    )
}

#[cfg(feature = "equilibrium-benchmark-prettyprint")]
fn optional_attempt_metric(time_us: f64, wall_us: f64, available: bool) -> String {
    if available {
        attempt_metric(time_us, wall_us)
    } else {
        "n/a".to_string()
    }
}

#[cfg(feature = "equilibrium-benchmark-prettyprint")]
fn backend_attempt_rows(case_name: &str, measurement: &Measurement) -> Vec<PrettyAttemptRow> {
    measurement
        .backend_attempt_rows
        .iter()
        .map(|attempt| {
            let wall = attempt.elapsed_us.max(0.0);
            PrettyAttemptRow {
                case_name: case_name.to_string(),
                scope: attempt.scope.clone(),
                profile: timing_profile(measurement),
                attempt: attempt.attempt_index.to_string(),
                backend: attempt.backend.clone(),
                outcome: attempt.outcome.clone(),
                metrics: if attempt.metrics_available {
                    if attempt.evaluation_timing_available {
                        "reported+callbacks".to_string()
                    } else {
                        "reported".to_string()
                    }
                } else if attempt.wall_time_available {
                    "wall-only".to_string()
                } else {
                    "unreported".to_string()
                },
                samples: format!("{:.1}", attempt.samples),
                wall: if attempt.metrics_available || attempt.wall_time_available {
                    format!("{:.2} ms", wall / 1_000.0)
                } else {
                    "n/a".to_string()
                },
                iterations: if attempt.metrics_available {
                    format!("{:.1}", attempt.iterations)
                } else {
                    "n/a".to_string()
                },
                residual: optional_attempt_metric(
                    attempt.residual_us,
                    wall,
                    attempt.evaluation_timing_available,
                ),
                jacobian: optional_attempt_metric(
                    attempt.jacobian_us,
                    wall,
                    attempt.evaluation_timing_available,
                ),
                linear: if attempt.linear_timing_available {
                    attempt_metric(attempt.linear_us, wall)
                } else {
                    "n/a".to_string()
                },
                overhead: optional_attempt_metric(
                    attempt.overhead_us,
                    wall,
                    attempt.evaluation_timing_available,
                ),
            }
        })
        .collect()
}

fn micros(duration: Duration) -> f64 {
    duration.as_secs_f64() * 1.0e6
}

fn measurement_from_reports(
    facade: &EquilibriumInstrumentationReport,
    timing: EquilibriumTimingReport,
    scalar_us: f64,
    enthalpy_us: f64,
    formulation_builds: f64,
) -> Measurement {
    let counters = facade.counters();
    let mut measurement = measurement_from_timing(timing, scalar_us, enthalpy_us);
    measurement.facade_us = micros(facade.elapsed());
    measurement.residual_evaluations = counters.residual_evaluations() as f64;
    measurement.jacobian_evaluations = counters.jacobian_evaluations() as f64;
    measurement.gibbs_vectors = counters.gibbs_vector_evaluations() as f64;
    measurement.enthalpy_vectors = counters.enthalpy_vector_evaluations() as f64;
    measurement.heat_capacity_vectors = counters.heat_capacity_vector_evaluations() as f64;
    measurement.formulation_builds = formulation_builds;
    measurement
}

fn add_facade_timing(measurement: &mut Measurement, timing: &EquilibriumTimingReport) {
    measurement.facade_resolution_us = micros(timing.facade_resolution());
    measurement.request_preparation_us = micros(timing.request_preparation());
    measurement.route_selection_us = micros(timing.route_selection());
    measurement.normalization_us = micros(timing.normalization());
    measurement.phase_control_us = micros(timing.phase_control());
    measurement.facade_thermochemistry_us = micros(timing.thermochemistry_preparation());
    measurement.result_assembly_us = micros(timing.result_assembly());
}

fn measurement_from_timing(
    timing: EquilibriumTimingReport,
    scalar_us: f64,
    enthalpy_us: f64,
) -> Measurement {
    Measurement {
        facade_us: 0.0,
        inner_total_us: micros(timing.total()),
        workflow_total_us: micros(timing.total()),
        route_selection_us: micros(timing.route_selection()),
        normalization_us: micros(timing.normalization()),
        phase_control_us: micros(timing.phase_control()),
        recovery_us: micros(timing.recovery()),
        scalar_us,
        enthalpy_us,
        thermochemistry_us: micros(timing.thermochemistry_preparation()),
        ph_formulation_preparation_us: micros(timing.ph_formulation_preparation()),
        numeric_closure_us: micros(timing.numeric_closure_construction()),
        symbolic_us: micros(timing.symbolic_construction()),
        equation_us: micros(timing.equation_construction()),
        numerical_problem_us: micros(timing.numerical_problem_preparation()),
        active_set_formulation_build_us: micros(timing.active_set_formulation_build()),
        rst_preparation_us: micros(timing.rst_preparation()),
        projection_build_us: micros(timing.projection_build()),
        nonlinear_us: micros(timing.nonlinear_solve()),
        validation_us: micros(timing.validation()),
        postprocessing_us: micros(timing.postprocessing()),
        repository_lookup_us: micros(timing.repository_lookup()),
        ..Measurement::default()
    }
}

fn instrumented_builder(builder: EquilibriumCalculatorBuilder) -> EquilibriumCalculatorBuilder {
    instrumented_builder_with_reference(builder, REFERENCE_PRESSURE_PA)
}

fn instrumented_builder_with_reference(
    builder: EquilibriumCalculatorBuilder,
    reference_pressure_pa: f64,
) -> EquilibriumCalculatorBuilder {
    instrumented_builder_with_route(builder, reference_pressure_pa, BackendRoute::Production)
}

fn instrumented_builder_with_route(
    builder: EquilibriumCalculatorBuilder,
    reference_pressure_pa: f64,
    route: BackendRoute,
) -> EquilibriumCalculatorBuilder {
    let builder = builder
        .prefer_libraries(["NASA_gas", "NASA_cond"])
        .offline_only()
        .pressure_pa(PRESSURE_PA)
        .reference_pressure_pa(reference_pressure_pa);
    let builder = match route {
        BackendRoute::Production => builder.production_cascade(),
        BackendRoute::LegacyLm => builder
            .solver_policy(SolverPolicy::Single(SolverBackend::Legacy(
                LegacyEquilibriumSolver::LM,
            )))
            .expect("single legacy-LM benchmark policy must validate"),
        BackendRoute::LegacyNr => builder
            .solver_policy(SolverPolicy::Single(SolverBackend::Legacy(
                LegacyEquilibriumSolver::NR,
            )))
            .expect("single legacy-NR benchmark policy must validate"),
        BackendRoute::LegacyTr => builder
            .solver_policy(SolverPolicy::Single(SolverBackend::Legacy(
                LegacyEquilibriumSolver::TR,
            )))
            .expect("single legacy-TR benchmark policy must validate"),
    };
    builder
        .preset(EquilibriumCalculatorPreset::Timed)
        .instrumentation(EquilibriumInstrumentationMode::Enabled)
}

fn reference_pressure(universe: Universe) -> f64 {
    match universe {
        Universe::Tp1907 => PRESSURE_PA,
        _ => REFERENCE_PRESSURE_PA,
    }
}

fn inventory(universe: Universe) -> ElementInventory {
    let amounts = match universe {
        Universe::SmallH2O => vec![("H", 4.0), ("O", 2.0)],
        Universe::Stanjan15 => vec![("C", 6.0), ("H", 16.0), ("O", 20.0), ("N", 75.2)],
        // Source-faithful TP-1907 CHON + Ar inventory reconstructed by the
        // frozen fixture (the printed order in its story report is Ar,C,H,N,O).
        Universe::Tp1907 => vec![
            ("Ar", 0.05364816),
            ("C", 1.00182736),
            ("H", 2.0),
            ("N", 8.9462088),
            ("O", 2.40365472),
        ],
        Universe::Rp1311_40 => vec![("Ar", 0.1), ("C", 0.1), ("H", 2.2), ("O", 4.3), ("N", 6.0)],
    };
    ElementInventory::from_amounts(amounts).expect("benchmark inventory must be valid")
}

fn temperature(universe: Universe) -> f64 {
    match universe {
        Universe::SmallH2O | Universe::Stanjan15 => REFERENCE_TEMPERATURE_K,
        Universe::Tp1907 => 700.0,
        Universe::Rp1311_40 => 2_400.0,
    }
}

fn temperature_bounds(universe: Universe) -> TemperatureBounds {
    match universe {
        // The frozen graphite case is intentionally kept inside one native
        // coefficient interval; crossing it is a capability error, not a
        // useful measurement of the P,H nonlinear path.
        Universe::Tp1907 => TemperatureBounds::new(699.0, 700.0),
        _ => TemperatureBounds::new(
            (temperature(universe) - 400.0).max(300.0),
            temperature(universe) + 400.0,
        ),
    }
    .expect("benchmark temperature bounds must be valid")
}

fn molecular_builder(universe: Universe) -> EquilibriumCalculatorBuilder {
    let builder = base_builder(universe);
    match universe {
        Universe::SmallH2O => builder.initial_moles([0.1, 0.05, 1.9]),
        Universe::Stanjan15 => builder.initial_phase_moles([
            (component(Some("gas"), "CH4"), 4.0),
            (component(Some("gas"), "CO2"), 2.0),
            (component(Some("gas"), "O2"), 8.0),
            (component(Some("gas"), "N2"), 37.6),
        ]),
        Universe::Tp1907 => builder.initial_phase_moles([
            (component(Some("gas"), "CH4"), 1.0),
            (component(Some("gas"), "N2"), 3.0),
            (component(Some("gas"), "O2"), 2.0),
        ]),
        Universe::Rp1311_40 => builder.initial_phase_moles([
            (component(None, "Ar"), 0.1),
            (component(None, "CO2"), 0.1),
            (component(None, "H2"), 1.0),
            (component(None, "H2O"), 0.1),
            (component(None, "N2"), 3.0),
            (component(None, "O2"), 2.0),
        ]),
    }
}

fn base_builder(universe: Universe) -> EquilibriumCalculatorBuilder {
    match universe {
        Universe::SmallH2O => EquilibriumCalculator::builder()
            .single_ideal_gas(["H2", "O2", "H2O"])
            .expect("small benchmark gas phase must be valid"),
        Universe::Stanjan15 => EquilibriumCalculator::builder()
            .ideal_gas_phase("gas", STANJAN_15)
            .expect("STANJAN frozen gas phase must be valid"),
        Universe::Tp1907 => EquilibriumCalculator::builder()
            .ideal_gas_phase("gas", TP1907_GAS_17)
            .expect("TP-1907 frozen gas phase must be valid")
            .pure_condensed_phase("graphite", "C(gr)", PhysicalState::Solid)
            .expect("TP-1907 frozen graphite phase must be valid"),
        Universe::Rp1311_40 => EquilibriumCalculator::builder()
            .single_ideal_gas(CEA_RP1311_40)
            .expect("CEA RP-1311 frozen gas phase must be valid"),
    }
}

fn point_builder(universe: Universe, elemental: bool) -> EquilibriumCalculatorBuilder {
    let builder = if elemental {
        base_builder(universe).element_inventory(inventory(universe))
    } else {
        molecular_builder(universe)
    };
    if elemental {
        builder.at_temperature(temperature(universe))
    } else {
        builder.at_temperature(REFERENCE_TEMPERATURE_K)
    }
}

fn ph_builder(
    universe: Universe,
    elemental: bool,
    target_enthalpy: f64,
) -> EquilibriumCalculatorBuilder {
    let builder = if elemental {
        base_builder(universe).element_inventory(inventory(universe))
    } else {
        molecular_builder(universe)
    };
    let builder = if matches!(universe, Universe::Tp1907 | Universe::Rp1311_40) {
        // The frozen high-dimensional P,H cases use the safeguarded nested
        // route: TP-1907 has a condensed-phase interval boundary, while the
        // 40-component case is better characterized by its scalar accepted
        // P,T manifold than by a monolithic enthalpy tolerance probe.
        builder.ph_solve_mode(PhSolveMode::NestedTemperature)
    } else {
        builder
    };
    let builder = if elemental && matches!(universe, Universe::Tp1907) {
        builder.element_numerical_seed(tp1907_element_seed())
    } else {
        builder
    };
    builder.at_total_enthalpy(
        TotalEnthalpyJoules::new(target_enthalpy).expect("target must be finite"),
        temperature(universe),
        temperature_bounds(universe),
    )
}

fn accepted_backend_name(report: &EquilibriumSolveReport) -> &'static str {
    match report.accepted_backend {
        SolverBackend::RustedSciThe(_) => "rst",
        SolverBackend::Legacy(_) => "legacy",
    }
}

fn accepted_solver_name(report: &EquilibriumSolveReport) -> &'static str {
    match report.accepted_backend {
        SolverBackend::Legacy(LegacyEquilibriumSolver::LM) => "legacy_lm",
        SolverBackend::Legacy(LegacyEquilibriumSolver::NR) => "legacy_nr",
        SolverBackend::Legacy(LegacyEquilibriumSolver::TR) => "legacy_tr",
        SolverBackend::RustedSciThe(solver) => solver.name(),
    }
}

fn add_solver_evidence(measurement: &mut Measurement, report: &EquilibriumSolveReport) {
    measurement.backend_attempts = report.started_attempt_count() as f64;
    measurement.nonlinear_iterations = report.nonlinear_iterations() as f64;
    measurement.accepted_backend = accepted_backend_name(report);
    measurement.solver_name = accepted_solver_name(report);
}

fn add_backend_timing(measurement: &mut Measurement, report: &EquilibriumSolveReport) {
    if let Some(attempt) = report
        .attempts
        .iter()
        .find(|attempt| attempt.backend == report.accepted_backend && attempt.is_accepted())
    {
        if let Some(metrics) = &attempt.metrics {
            measurement.backend_us += metrics.elapsed_micros as f64;
            measurement.solver_iterations += metrics.iterations as f64;
            measurement.solver_residual_calls += metrics.residual_evaluations as f64;
            measurement.solver_jacobian_calls += metrics.jacobian_evaluations as f64;
            measurement.solver_linear_system_solves += metrics.linear_solves as f64;
            if let Some(timing) = &metrics.evaluation_timing {
                measurement.solver_residual_us += timing.residual_evaluation_micros as f64;
                measurement.solver_jacobian_us += timing.jacobian_evaluation_micros as f64;
                if let Some(linear_system_us) = timing.linear_system_micros {
                    measurement.solver_linear_system_us += linear_system_us as f64;
                    measurement.solver_linear_timing_available = true;
                }
                measurement.solver_overhead_us += timing.solver_overhead_micros as f64;
                measurement.backend_callback_us +=
                    (timing.residual_evaluation_micros + timing.jacobian_evaluation_micros) as f64;
                measurement.backend_overhead_us += timing.solver_overhead_micros as f64;
            }
        }
    }
}

#[cfg(feature = "equilibrium-benchmark-prettyprint")]
fn add_backend_attempt_rows(
    measurement: &mut Measurement,
    scope: impl Into<String>,
    report: &EquilibriumSolveReport,
) {
    add_backend_attempt_slice(measurement, scope, &report.attempts);
}

#[cfg(feature = "equilibrium-benchmark-prettyprint")]
fn add_phase_control_attempt_rows(
    measurement: &mut Measurement,
    scope: impl Into<String>,
    report: &PhaseControlledSolveReport,
) {
    let scope = scope.into();
    for (iteration, nonlinear_report) in report.nonlinear_reports.iter().enumerate() {
        add_backend_attempt_rows(
            measurement,
            format!("{scope}/active-set-{iteration}"),
            nonlinear_report,
        );
    }
}

#[cfg(feature = "equilibrium-benchmark-prettyprint")]
fn add_multistart_attempt_rows(
    measurement: &mut Measurement,
    scope: impl Into<String>,
    report: &MultiStartSolveReport,
) {
    let scope = scope.into();
    for (seed_index, attempts) in report.backend_attempts().iter().enumerate() {
        add_backend_attempt_slice(measurement, format!("{scope}/seed-{seed_index}"), attempts);
    }
}

#[cfg(feature = "equilibrium-benchmark-prettyprint")]
fn add_backend_attempt_slice(
    measurement: &mut Measurement,
    scope: impl Into<String>,
    attempts: &[SolverAttemptReport],
) {
    let scope = scope.into();
    for (attempt_index, attempt) in attempts.iter().enumerate() {
        let (outcome, _) = match &attempt.outcome {
            SolverAttemptOutcome::Accepted => ("accepted", None),
            SolverAttemptOutcome::Failed { reason, .. } => ("failed", Some(reason)),
            SolverAttemptOutcome::RejectedCandidate { reason } => {
                ("rejected_candidate", Some(reason))
            }
            SolverAttemptOutcome::Skipped { reason } => ("skipped", Some(reason)),
        };
        let (
            elapsed_us,
            iterations,
            residual_calls,
            jacobian_calls,
            linear_solves,
            residual_us,
            jacobian_us,
            linear_us,
            overhead_us,
            linear_timing_available,
        ) = if let Some(metrics) = &attempt.metrics {
            let timing = metrics.evaluation_timing.as_ref();
            (
                metrics.elapsed_micros as f64,
                metrics.iterations as f64,
                metrics.residual_evaluations as f64,
                metrics.jacobian_evaluations as f64,
                metrics.linear_solves as f64,
                timing.map_or(0.0, |timing| timing.residual_evaluation_micros as f64),
                timing.map_or(0.0, |timing| timing.jacobian_evaluation_micros as f64),
                timing
                    .and_then(|timing| timing.linear_system_micros)
                    .unwrap_or(0) as f64,
                timing.map_or(0.0, |timing| timing.solver_overhead_micros as f64),
                timing
                    .and_then(|timing| timing.linear_system_micros)
                    .is_some(),
            )
        } else {
            let elapsed_us = match &attempt.outcome {
                SolverAttemptOutcome::Failed { elapsed_micros, .. } => *elapsed_micros as f64,
                _ => 0.0,
            };
            (elapsed_us, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, false)
        };
        let metrics_available = attempt.metrics.is_some();
        let wall_time_available =
            metrics_available || matches!(&attempt.outcome, SolverAttemptOutcome::Failed { .. });
        let evaluation_timing_available = attempt
            .metrics
            .as_ref()
            .and_then(|metrics| metrics.evaluation_timing.as_ref())
            .is_some();
        measurement
            .backend_attempt_rows
            .push(BackendAttemptMeasurement {
                scope: scope.clone(),
                attempt_index,
                backend: format!("{:?}", attempt.backend),
                outcome: outcome.to_string(),
                samples: 1.0,
                elapsed_us,
                iterations,
                residual_calls,
                jacobian_calls,
                linear_solves,
                residual_us,
                jacobian_us,
                linear_us,
                overhead_us,
                metrics_available,
                wall_time_available,
                evaluation_timing_available,
                linear_timing_available,
            });
    }
}

fn add_ph_route_evidence(measurement: &mut Measurement, report: &PhTemperatureSolveReport) {
    measurement.ph_route = match report.solve_path() {
        PhSolvePath::MonolithicFixedActiveSet => "monolithic",
        PhSolvePath::MonolithicPhaseControl => "monolithic_phase_control",
        PhSolvePath::NestedTemperature => "nested",
    };
    measurement.ph_route_decisions = report.route_decisions().len() as f64;
}

fn add_ph_timing_evidence(measurement: &mut Measurement, report: &PhTemperatureSolveReport) {
    let timing = report.timing();
    measurement.ph_trial_total_us = micros(timing.trial_total());
    measurement.ph_trial_inner_execution_us = micros(timing.trial_inner_execution());
    measurement.ph_trial_preparation_us = micros(timing.trial_preparation());
    measurement.ph_trial_template_preparation_us = micros(timing.trial_template_preparation());
    measurement.ph_trial_postprocessing_us = micros(timing.trial_postprocessing());
    measurement.ph_bracket_orchestration_us = micros(timing.bracket_orchestration());
    measurement.ph_monolithic_preparation_us = micros(timing.monolithic_preparation());
    measurement.ph_monolithic_rst_preparation_us = micros(timing.monolithic_rst_preparation());
    measurement.ph_monolithic_backend_execution_us = micros(timing.monolithic_backend_execution());
    measurement.ph_monolithic_execution_us = micros(timing.monolithic_execution());
    measurement.ph_monolithic_publication_us = micros(timing.monolithic_publication());
    measurement.ph_monolithic_validation_us = report
        .monolithic_evidence()
        .map_or(0.0, |evidence| micros(evidence.validation_duration()));
    for (trial_index, trial) in report.trials().iter().enumerate() {
        #[cfg(not(feature = "equilibrium-benchmark-prettyprint"))]
        let _ = trial_index;
        if let Some(inner) = trial.inner_evidence() {
            add_backend_timing(measurement, inner.solve_report());
            #[cfg(feature = "equilibrium-benchmark-prettyprint")]
            if let Some(recovery) = inner.extensive_normalization_recovery() {
                add_backend_attempt_slice(
                    measurement,
                    format!("trial-{trial_index}/recovery/physical-failure"),
                    &recovery.physical_failure_attempts,
                );
                add_backend_attempt_slice(
                    measurement,
                    format!("trial-{trial_index}/recovery/normalized-discovery"),
                    &recovery.normalized_discovery_attempts,
                );
                add_backend_attempt_slice(
                    measurement,
                    format!("trial-{trial_index}/recovery/physical-retry"),
                    &recovery.physical_retry_attempts,
                );
            } else if let Some(phase_control) = inner.phase_control_report() {
                add_phase_control_attempt_rows(
                    measurement,
                    format!("trial-{trial_index}/phase-control"),
                    phase_control,
                );
            } else if let Some(multi_start) = inner.multi_start_report() {
                add_multistart_attempt_rows(
                    measurement,
                    format!("trial-{trial_index}"),
                    multi_start,
                );
            } else {
                add_backend_attempt_rows(
                    measurement,
                    format!("trial-{trial_index}"),
                    inner.solve_report(),
                );
            }
        }
    }
}

#[cfg(feature = "equilibrium-benchmark-prettyprint")]
fn add_ph_recovery_evidence(
    measurement: &mut Measurement,
    equilibrium: &MultiphaseEquilibriumSolution,
    recovery_scope: Duration,
) {
    if let Some(recovery) = equilibrium.extensive_normalization_recovery() {
        // The solution evidence owns child timings; the enclosing P,H report
        // owns the complete transaction wall time, including the failed
        // physical attempt that authorized recovery.
        measurement.recovery_us = micros(recovery_scope);
        measurement.recovery_physical_failure_us = micros(recovery.physical_failure_duration);
        measurement.recovery_normalization_us = micros(recovery.normalization_duration);
        measurement.recovery_normalized_discovery_us =
            micros(recovery.normalized_discovery_duration);
        measurement.recovery_physical_retry_us = micros(recovery.physical_retry_duration);
        measurement.recovery_reconstruction_us = micros(recovery.reconstruction_duration);
        add_backend_attempt_slice(
            measurement,
            "P,H recovery/physical-failure",
            &recovery.physical_failure_attempts,
        );
        add_backend_attempt_slice(
            measurement,
            "P,H recovery/normalized-discovery",
            &recovery.normalized_discovery_attempts,
        );
        add_backend_attempt_slice(
            measurement,
            "P,H recovery/physical-retry",
            &recovery.physical_retry_attempts,
        );
    }
}

fn solve_scenario(scenario: Scenario, target_enthalpy: f64) -> Measurement {
    match scenario.task {
        ScenarioTask::PtMolecular | ScenarioTask::PtElementInventory => {
            let elemental = matches!(scenario.task, ScenarioTask::PtElementInventory);
            let outcome = instrumented_builder_with_route(
                point_builder(scenario.universe, elemental),
                reference_pressure(scenario.universe),
                scenario.route,
            )
            .solve()
            .expect("frozen P,T benchmark must solve");
            let EquilibriumCalculatorOutcome::PtPoint(point) = outcome else {
                panic!("expected frozen P,T point")
            };
            let mut measurement = measurement_from_reports(
                point.instrumentation_report(),
                *point.solution().timing_report(),
                0.0,
                0.0,
                1.0,
            );
            add_facade_timing(&mut measurement, point.facade_timing());
            add_solver_evidence(&mut measurement, point.solution().solve_report());
            add_backend_timing(&mut measurement, point.solution().solve_report());
            #[cfg(feature = "equilibrium-benchmark-prettyprint")]
            if let Some(recovery) = point.solution().extensive_normalization_recovery() {
                add_backend_attempt_slice(
                    &mut measurement,
                    "recovery/physical-failure",
                    &recovery.physical_failure_attempts,
                );
                add_backend_attempt_slice(
                    &mut measurement,
                    "recovery/normalized-discovery",
                    &recovery.normalized_discovery_attempts,
                );
                add_backend_attempt_slice(
                    &mut measurement,
                    "recovery/physical-retry",
                    &recovery.physical_retry_attempts,
                );
            } else if let Some(phase_control) = point.solution().phase_control_report() {
                add_phase_control_attempt_rows(
                    &mut measurement,
                    "point/phase-control",
                    phase_control,
                );
            } else if let Some(multi_start) = point.solution().multi_start_report() {
                add_multistart_attempt_rows(&mut measurement, "point", multi_start);
            } else {
                add_backend_attempt_rows(
                    &mut measurement,
                    "point",
                    point.solution().solve_report(),
                );
            }
            measurement
        }
        ScenarioTask::PhMolecular | ScenarioTask::PhElementInventory => {
            let elemental = matches!(scenario.task, ScenarioTask::PhElementInventory);
            let outcome = instrumented_builder_with_route(
                ph_builder(scenario.universe, elemental, target_enthalpy),
                reference_pressure(scenario.universe),
                scenario.route,
            )
            .solve()
            .expect("frozen P,H benchmark must solve");
            let EquilibriumCalculatorOutcome::PhPoint(point) = outcome else {
                panic!("expected frozen P,H point")
            };
            let report = point.solution().report();
            let outer = report.timing();
            let mut measurement = measurement_from_reports(
                point.instrumentation_report(),
                report.inner_timing(),
                micros(outer.scalar_orchestration()),
                micros(outer.enthalpy_evaluation()),
                report.fixed_formulation_builds() as f64,
            );
            // For P,H, the attached inner timing is the fixed-P,T solve. Use
            // the complete P,H route for the top-level facade-gap partition.
            measurement.workflow_total_us = micros(outer.total());
            add_facade_timing(&mut measurement, point.facade_timing());
            add_ph_route_evidence(&mut measurement, report);
            add_ph_timing_evidence(&mut measurement, report);
            #[cfg(feature = "equilibrium-benchmark-prettyprint")]
            add_ph_recovery_evidence(
                &mut measurement,
                point.solution().equilibrium(),
                report.inner_timing().recovery(),
            );
            add_solver_evidence(
                &mut measurement,
                point.solution().equilibrium().solve_report(),
            );
            if !matches!(report.solve_path(), PhSolvePath::NestedTemperature) {
                add_backend_timing(
                    &mut measurement,
                    point.solution().equilibrium().solve_report(),
                );
                #[cfg(feature = "equilibrium-benchmark-prettyprint")]
                if let Some(recovery) = point
                    .solution()
                    .equilibrium()
                    .extensive_normalization_recovery()
                {
                    add_backend_attempt_slice(
                        &mut measurement,
                        "recovery/physical-failure",
                        &recovery.physical_failure_attempts,
                    );
                    add_backend_attempt_slice(
                        &mut measurement,
                        "recovery/normalized-discovery",
                        &recovery.normalized_discovery_attempts,
                    );
                    add_backend_attempt_slice(
                        &mut measurement,
                        "recovery/physical-retry",
                        &recovery.physical_retry_attempts,
                    );
                } else if let Some(phase_control) =
                    point.solution().equilibrium().phase_control_report()
                {
                    add_phase_control_attempt_rows(
                        &mut measurement,
                        "monolithic/phase-control",
                        phase_control,
                    );
                } else if let Some(multi_start) =
                    point.solution().equilibrium().multi_start_report()
                {
                    add_multistart_attempt_rows(&mut measurement, "monolithic", multi_start);
                } else {
                    add_backend_attempt_rows(
                        &mut measurement,
                        "monolithic",
                        point.solution().equilibrium().solve_report(),
                    );
                }
            }
            measurement
        }
        ScenarioTask::PtElementInventoryRange => {
            assert!(matches!(scenario.universe, Universe::SmallH2O));
            let outcome = instrumented_builder(
                EquilibriumCalculator::builder()
                    .single_ideal_gas(["H2", "O2", "H2O"])
                    .expect("benchmark gas phase must be valid")
                    .element_inventory(inventory(Universe::SmallH2O))
                    .over_temperature_range(RANGE_TEMPERATURES.to_vec())
                    .expect("benchmark P,T range must be valid"),
            )
            .solve()
            .expect("elemental P,T range benchmark must solve");
            let EquilibriumCalculatorOutcome::PtRange(range) = outcome else {
                panic!("expected elemental P,T range")
            };
            let points = range.solution().points();
            let divisor = points.len() as f64;
            let mut measurement = Measurement {
                facade_us: micros(range.instrumentation_report().elapsed()),
                workflow_total_us: micros(range.solution().report().total()),
                formulation_builds: range.solution().report().formulation_builds() as f64,
                ..Measurement::default()
            };
            add_facade_timing(&mut measurement, range.facade_timing());
            for (point_index, point) in points.iter().enumerate() {
                #[cfg(not(feature = "equilibrium-benchmark-prettyprint"))]
                let _ = point_index;
                measurement.add_assign(measurement_from_timing(*point.report().timing(), 0.0, 0.0));
                #[cfg(feature = "equilibrium-benchmark-prettyprint")]
                if let Some(recovery) = point.solution().extensive_normalization_recovery() {
                    add_backend_attempt_slice(
                        &mut measurement,
                        format!("point-{point_index}/recovery/physical-failure"),
                        &recovery.physical_failure_attempts,
                    );
                    add_backend_attempt_slice(
                        &mut measurement,
                        format!("point-{point_index}/recovery/normalized-discovery"),
                        &recovery.normalized_discovery_attempts,
                    );
                    add_backend_attempt_slice(
                        &mut measurement,
                        format!("point-{point_index}/recovery/physical-retry"),
                        &recovery.physical_retry_attempts,
                    );
                } else if let Some(phase_control) = point.solution().phase_control_report() {
                    add_phase_control_attempt_rows(
                        &mut measurement,
                        format!("point-{point_index}/phase-control"),
                        phase_control,
                    );
                } else if let Some(multi_start) = point.solution().multi_start_report() {
                    add_multistart_attempt_rows(
                        &mut measurement,
                        format!("point-{point_index}"),
                        multi_start,
                    );
                } else {
                    add_backend_attempt_rows(
                        &mut measurement,
                        format!("point-{point_index}"),
                        point.solution().solve_report(),
                    );
                }
            }
            // Point reports provide stage detail; the range report owns the
            // complete workflow wall time used for the facade-gap partition.
            measurement.workflow_total_us = micros(range.solution().report().total());
            let counters = range.instrumentation_report().counters();
            measurement.residual_evaluations = counters.residual_evaluations() as f64;
            measurement.jacobian_evaluations = counters.jacobian_evaluations() as f64;
            measurement.gibbs_vectors = counters.gibbs_vector_evaluations() as f64;
            measurement.enthalpy_vectors = counters.enthalpy_vector_evaluations() as f64;
            measurement.heat_capacity_vectors = counters.heat_capacity_vector_evaluations() as f64;
            for point in range.solution().points() {
                add_solver_evidence(&mut measurement, point.solution().solve_report());
                add_backend_timing(&mut measurement, point.solution().solve_report());
            }
            measurement.divided_by(divisor)
        }
        ScenarioTask::PhElementInventoryRange => {
            assert!(matches!(scenario.universe, Universe::SmallH2O));
            let outcome = instrumented_builder(
                EquilibriumCalculator::builder()
                    .single_ideal_gas(["H2", "O2", "H2O"])
                    .expect("benchmark gas phase must be valid")
                    .element_inventory(inventory(Universe::SmallH2O))
                    .over_enthalpy_range(
                        [
                            target_enthalpy,
                            target_enthalpy + 1_000.0,
                            target_enthalpy + 2_000.0,
                        ],
                        2_300.0,
                        TemperatureBounds::new(2_100.0, 2_900.0)
                            .expect("benchmark temperature bounds must be valid"),
                    )
                    .expect("benchmark P,H range must be valid"),
            )
            .solve()
            .expect("elemental P,H range benchmark must solve");
            let EquilibriumCalculatorOutcome::PhRange(range) = outcome else {
                panic!("expected elemental P,H range")
            };
            let points = range.solution().points();
            let divisor = points.len() as f64;
            let mut measurement = Measurement {
                facade_us: micros(range.instrumentation_report().elapsed()),
                workflow_total_us: micros(range.solution().report().total()),
                formulation_builds: range.solution().report().formulation_builds() as f64,
                ..Measurement::default()
            };
            add_facade_timing(&mut measurement, range.facade_timing());
            for (point_index, point) in points.iter().enumerate() {
                #[cfg(not(feature = "equilibrium-benchmark-prettyprint"))]
                let _ = point_index;
                let report = point.solution().report();
                let outer = report.timing();
                add_ph_route_evidence(&mut measurement, report);
                measurement.ph_trial_total_us += micros(outer.trial_total());
                measurement.ph_trial_inner_execution_us += micros(outer.trial_inner_execution());
                measurement.ph_trial_preparation_us += micros(outer.trial_preparation());
                measurement.ph_bracket_orchestration_us += micros(outer.bracket_orchestration());
                measurement.add_assign(measurement_from_timing(
                    report.inner_timing(),
                    micros(outer.scalar_orchestration()),
                    micros(outer.enthalpy_evaluation()),
                ));
                let mut point_timing = Measurement::default();
                add_ph_timing_evidence(&mut point_timing, report);
                #[cfg(feature = "equilibrium-benchmark-prettyprint")]
                add_ph_recovery_evidence(
                    &mut point_timing,
                    point.solution().equilibrium(),
                    report.inner_timing().recovery(),
                );
                measurement.add_assign(point_timing);

                // Monolithic P,H has no scalar trials, so its backend cascade
                // must be attached directly to the point scope.
                #[cfg(feature = "equilibrium-benchmark-prettyprint")]
                if report.trials().is_empty() {
                    let equilibrium = point.solution().equilibrium();
                    add_backend_timing(&mut measurement, equilibrium.solve_report());
                    if let Some(recovery) = equilibrium.extensive_normalization_recovery() {
                        add_backend_attempt_slice(
                            &mut measurement,
                            format!("point-{point_index}/recovery/physical-failure"),
                            &recovery.physical_failure_attempts,
                        );
                        add_backend_attempt_slice(
                            &mut measurement,
                            format!("point-{point_index}/recovery/normalized-discovery"),
                            &recovery.normalized_discovery_attempts,
                        );
                        add_backend_attempt_slice(
                            &mut measurement,
                            format!("point-{point_index}/recovery/physical-retry"),
                            &recovery.physical_retry_attempts,
                        );
                    } else if let Some(phase_control) = equilibrium.phase_control_report() {
                        add_phase_control_attempt_rows(
                            &mut measurement,
                            format!("point-{point_index}/monolithic/phase-control"),
                            phase_control,
                        );
                    } else if let Some(multi_start) = equilibrium.multi_start_report() {
                        add_multistart_attempt_rows(
                            &mut measurement,
                            format!("point-{point_index}/monolithic"),
                            multi_start,
                        );
                    } else {
                        add_backend_attempt_rows(
                            &mut measurement,
                            format!("point-{point_index}/monolithic"),
                            equilibrium.solve_report(),
                        );
                    }
                }
                measurement.formulation_builds += point.report().formulation_builds() as f64;
            }
            // The range report includes continuation orchestration around all
            // accepted points, which no individual P,H point report can show.
            measurement.workflow_total_us = micros(range.solution().report().total());
            let counters = range.instrumentation_report().counters();
            measurement.residual_evaluations = counters.residual_evaluations() as f64;
            measurement.jacobian_evaluations = counters.jacobian_evaluations() as f64;
            measurement.gibbs_vectors = counters.gibbs_vector_evaluations() as f64;
            measurement.enthalpy_vectors = counters.enthalpy_vector_evaluations() as f64;
            measurement.heat_capacity_vectors = counters.heat_capacity_vector_evaluations() as f64;
            for point in range.solution().points() {
                add_solver_evidence(
                    &mut measurement,
                    point.solution().equilibrium().solve_report(),
                );
            }
            measurement.divided_by(divisor)
        }
    }
}

fn reference_enthalpy_for(universe: Universe, elemental: bool) -> f64 {
    if matches!(universe, Universe::Tp1907) {
        // Frozen TP-1906 source target used by the accepted P,H story. It
        // deliberately recovers about 699.509 K rather than defining H as
        // the equilibrium value at exactly 700 K.
        return -420_039.2761;
    }
    let outcome = instrumented_builder_with_reference(
        point_builder(universe, elemental),
        reference_pressure(universe),
    )
    .solve()
    .expect("reference P,T solve must provide the P,H target");
    let EquilibriumCalculatorOutcome::PtPoint(point) = outcome else {
        panic!("expected reference P,T point")
    };
    let thermochemistry =
        ResolvedThermochemistry::from_resolved_system(point.resolved()).expect("thermochemistry");
    thermochemistry
        .enthalpy_model()
        .evaluate_total(point.solution().component_moles(), temperature(universe))
        .expect("reference enthalpy must evaluate")
}

fn target_enthalpy(universe: Universe, targets: &[(Universe, f64)]) -> f64 {
    targets
        .iter()
        .find(|(candidate, _)| *candidate == universe)
        .map(|(_, target)| *target)
        .expect("target enthalpy must exist for benchmark universe")
}

fn baseline_scenarios() -> Vec<Scenario> {
    let mut scenarios = Vec::with_capacity(23);
    scenarios.extend([
        Scenario {
            universe: Universe::SmallH2O,
            task: ScenarioTask::PtMolecular,
            route: BackendRoute::Production,
        },
        Scenario {
            universe: Universe::SmallH2O,
            task: ScenarioTask::PtElementInventory,
            route: BackendRoute::Production,
        },
        Scenario {
            universe: Universe::SmallH2O,
            task: ScenarioTask::PhMolecular,
            route: BackendRoute::Production,
        },
        Scenario {
            universe: Universe::SmallH2O,
            task: ScenarioTask::PhElementInventory,
            route: BackendRoute::Production,
        },
        // Fixed legacy policies are deliberately retained as diagnostic
        // controls. They expose method-specific regressions that a healthy
        // production cascade can hide by accepting an RST attempt first.
        Scenario {
            universe: Universe::SmallH2O,
            task: ScenarioTask::PtElementInventory,
            route: BackendRoute::LegacyLm,
        },
        Scenario {
            universe: Universe::SmallH2O,
            task: ScenarioTask::PtElementInventory,
            route: BackendRoute::LegacyNr,
        },
        Scenario {
            universe: Universe::SmallH2O,
            task: ScenarioTask::PtElementInventory,
            route: BackendRoute::LegacyTr,
        },
        Scenario {
            universe: Universe::SmallH2O,
            task: ScenarioTask::PhElementInventory,
            route: BackendRoute::LegacyLm,
        },
        Scenario {
            universe: Universe::Stanjan15,
            task: ScenarioTask::PhElementInventory,
            route: BackendRoute::LegacyLm,
        },
        Scenario {
            universe: Universe::Stanjan15,
            task: ScenarioTask::PhElementInventory,
            route: BackendRoute::LegacyTr,
        },
        Scenario {
            universe: Universe::Tp1907,
            task: ScenarioTask::PhElementInventory,
            route: BackendRoute::LegacyTr,
        },
        Scenario {
            universe: Universe::Rp1311_40,
            task: ScenarioTask::PhElementInventory,
            route: BackendRoute::LegacyLm,
        },
        Scenario {
            universe: Universe::Rp1311_40,
            task: ScenarioTask::PhElementInventory,
            route: BackendRoute::LegacyNr,
        },
        Scenario {
            universe: Universe::Rp1311_40,
            task: ScenarioTask::PhElementInventory,
            route: BackendRoute::LegacyTr,
        },
    ]);
    for universe in [Universe::Stanjan15, Universe::Tp1907, Universe::Rp1311_40] {
        scenarios.extend([
            Scenario {
                universe,
                task: ScenarioTask::PtElementInventory,
                route: BackendRoute::Production,
            },
            Scenario {
                universe,
                task: ScenarioTask::PhElementInventory,
                route: BackendRoute::Production,
            },
        ]);
    }
    scenarios.extend([
        Scenario {
            universe: Universe::Stanjan15,
            task: ScenarioTask::PtElementInventory,
            route: BackendRoute::LegacyLm,
        },
        Scenario {
            universe: Universe::Tp1907,
            task: ScenarioTask::PtElementInventory,
            route: BackendRoute::LegacyLm,
        },
        Scenario {
            universe: Universe::Rp1311_40,
            task: ScenarioTask::PtElementInventory,
            route: BackendRoute::LegacyLm,
        },
        Scenario {
            universe: Universe::Rp1311_40,
            task: ScenarioTask::PtElementInventory,
            route: BackendRoute::LegacyNr,
        },
        Scenario {
            universe: Universe::Rp1311_40,
            task: ScenarioTask::PtElementInventory,
            route: BackendRoute::LegacyTr,
        },
    ]);
    // Ranges measure continuation reuse and remain on the small fixture so
    // the expanded component-count baseline does not hide that separate fact.
    scenarios.extend([
        Scenario {
            universe: Universe::SmallH2O,
            task: ScenarioTask::PtElementInventoryRange,
            route: BackendRoute::Production,
        },
        Scenario {
            universe: Universe::SmallH2O,
            task: ScenarioTask::PhElementInventoryRange,
            route: BackendRoute::Production,
        },
    ]);
    scenarios
}

#[cfg(feature = "equilibrium-benchmark-prettyprint")]
fn assert_pretty_coverage(
    scenarios: &[Scenario],
    macro_rows: &[PrettyMacroMeasurement],
    facade_rows: &[PrettyTimingRow],
    ph_rows: &[PrettyTimingRow],
    inner_rows: &[PrettyTimingRow],
    nonlinear_rows: &[PrettyNonlinearRow],
    attempt_rows: &[PrettyAttemptRow],
    solver_rows: &[PrettySolverMeasurement],
    legacy_ph_rows: &[PrettyLegacyPhEnvelopeRow],
) {
    for scenario in scenarios {
        let case = scenario.label();
        assert_eq!(
            macro_rows
                .iter()
                .filter(|row| row.case_name == case)
                .count(),
            1,
            "pretty macro table must cover {case}"
        );
        assert!(
            facade_rows.iter().any(|row| row.case_name == case),
            "pretty facade table must cover {case}"
        );
        assert_eq!(
            nonlinear_rows
                .iter()
                .filter(|row| row.case_name == case)
                .count(),
            1,
            "pretty nonlinear table must cover {case}"
        );
        assert!(
            attempt_rows.iter().any(|row| row.case_name == case),
            "pretty backend-attempt table must cover {case}"
        );
        assert_eq!(
            solver_rows
                .iter()
                .filter(|row| row.case_name == case)
                .count(),
            1,
            "pretty solver table must cover {case}"
        );
        if !matches!(scenario.route, BackendRoute::Production)
            && matches!(
                scenario.task,
                ScenarioTask::PhMolecular | ScenarioTask::PhElementInventory
            )
        {
            assert_eq!(
                legacy_ph_rows
                    .iter()
                    .filter(|row| row.case_name == case)
                    .count(),
                1,
                "Legacy P,H envelope table must cover {case}"
            );
        }

        let has_inner = inner_rows.iter().any(|row| row.case_name == case);
        let has_trial = ph_rows
            .iter()
            .any(|row| row.case_name == case && row.stage == "trial total");
        let has_monolithic_execution = ph_rows
            .iter()
            .any(|row| row.case_name == case && row.stage == "monolithic execution");
        let has_monolithic_backend = ph_rows
            .iter()
            .any(|row| row.case_name == case && row.stage == "monolithic backend execution");
        match scenario.task {
            ScenarioTask::PtMolecular
            | ScenarioTask::PtElementInventory
            | ScenarioTask::PtElementInventoryRange => {
                assert!(has_inner, "fixed-P,T detail must cover {case}")
            }
            ScenarioTask::PhMolecular
            | ScenarioTask::PhElementInventory
            | ScenarioTask::PhElementInventoryRange => {
                assert!(
                    has_trial || has_monolithic_execution || has_monolithic_backend,
                    "P,H detail must expose nested trials or monolithic backend execution for {case}"
                );
                if has_monolithic_execution || has_monolithic_backend {
                    assert!(
                        !has_inner,
                        "monolithic P,H must not fabricate an inner fixed-P,T table for {case}"
                    );
                } else {
                    assert!(
                        has_inner,
                        "nested P,H must expose its fixed-P,T inner detail for {case}"
                    );
                }
            }
        }
    }
}

fn print_baseline_table(targets: &[(Universe, f64)]) {
    println!(
        "\nEquilibrium facade baseline ({} repetitions)",
        BASELINE_REPETITIONS
    );
    #[cfg(not(feature = "equilibrium-benchmark-prettyprint"))]
    println!("case backend route wall_ms outside_ms outside_% workflow_ms workflow_% partition");
    #[cfg(feature = "equilibrium-benchmark-prettyprint")]
    let mut pretty_solver_rows = Vec::new();
    #[cfg(feature = "equilibrium-benchmark-prettyprint")]
    let mut pretty_legacy_ph_rows = Vec::new();
    #[cfg(feature = "equilibrium-benchmark-prettyprint")]
    let mut pretty_macro_rows = Vec::new();
    #[cfg(feature = "equilibrium-benchmark-prettyprint")]
    let mut pretty_facade_rows = Vec::new();
    #[cfg(feature = "equilibrium-benchmark-prettyprint")]
    let mut pretty_ph_rows = Vec::new();
    #[cfg(feature = "equilibrium-benchmark-prettyprint")]
    let mut pretty_inner_rows = Vec::new();
    #[cfg(feature = "equilibrium-benchmark-prettyprint")]
    let mut pretty_nonlinear_rows = Vec::new();
    #[cfg(feature = "equilibrium-benchmark-prettyprint")]
    let mut pretty_attempt_rows = Vec::new();
    let scenarios = baseline_scenarios();
    println!(
        "Dashboard scenarios: {} ({} repetitions per scenario; range rows are per point)",
        scenarios.len(),
        BASELINE_REPETITIONS
    );
    for scenario in scenarios.iter().copied() {
        let mut total = Measurement::default();
        for _ in 0..BASELINE_REPETITIONS {
            total.add_assign(solve_scenario(
                scenario,
                target_enthalpy(scenario.universe, targets),
            ));
        }
        let mean = total.divided_by(BASELINE_REPETITIONS as f64);
        #[cfg(feature = "equilibrium-benchmark-prettyprint")]
        {
            let case_name = scenario.label();
            pretty_macro_rows.push(PrettyMacroMeasurement::from_measurement(
                case_name.clone(),
                &mean,
            ));
            pretty_solver_rows.push(PrettySolverMeasurement::from_measurement(
                case_name.clone(),
                &mean,
            ));
            pretty_facade_rows.extend(facade_timing_rows(&case_name, &mean));
            pretty_inner_rows.extend(inner_timing_rows(&case_name, &mean));
            pretty_nonlinear_rows.push(nonlinear_rows(&case_name, &mean));
            pretty_ph_rows.extend(ph_timing_rows(&case_name, &mean));
            pretty_attempt_rows.extend(backend_attempt_rows(&case_name, &mean));
            if !matches!(scenario.route, BackendRoute::Production)
                && matches!(
                    scenario.task,
                    ScenarioTask::PhMolecular | ScenarioTask::PhElementInventory
                )
            {
                pretty_legacy_ph_rows.push(
                    PrettyLegacyPhEnvelopeRow::from_measurement(case_name.clone(), &mean)
                        .expect("Legacy P,H benchmark must expose an accepted backend"),
                );
            }
        }
        #[cfg(not(feature = "equilibrium-benchmark-prettyprint"))]
        {
            let facade_outside_us = (mean.facade_us - mean.workflow_total_us).max(0.0);
            let outside_percent = percentage(facade_outside_us, mean.facade_us);
            let workflow_percent = percentage(mean.workflow_total_us, mean.facade_us);
            println!(
                "{:<44} {:<24} {:<12} {:>8.2} {:>10.2} {:>8.1}% {:>11.2} {:>9.1}% {:>9.1}%",
                scenario.label(),
                mean.accepted_backend,
                mean.ph_route,
                mean.facade_us / 1_000.0,
                facade_outside_us / 1_000.0,
                outside_percent,
                mean.workflow_total_us / 1_000.0,
                workflow_percent,
                outside_percent + workflow_percent,
            );
        }
    }
    println!(
        "All range rows are per accepted point; wall and workflow times are complete-facade measurements."
    );
    #[cfg(feature = "equilibrium-benchmark-prettyprint")]
    {
        assert_pretty_coverage(
            &scenarios,
            &pretty_macro_rows,
            &pretty_facade_rows,
            &pretty_ph_rows,
            &pretty_inner_rows,
            &pretty_nonlinear_rows,
            &pretty_attempt_rows,
            &pretty_solver_rows,
            &pretty_legacy_ph_rows,
        );
        println!("\nEquilibrium macro time budget");
        println!(
            "Only facade-outside-workflow + workflow are shown here; they are mutually exclusive and partition the facade wall time."
        );
        println!("{}", Table::new(pretty_macro_rows));

        println!("\nFacade detail (percent of facade wall time)");
        println!(
            "Rows are measurements of facade-owned operations; this table is a detail view, not an additive budget."
        );
        println!("{}", Table::new(pretty_facade_rows));

        println!("\nP,H workflow detail (percent of complete P,H workflow wall time)");
        println!(
            "Rows are inclusive parent/child observations. Trial total contains trial inner execution, and monolithic execution contains backend work."
        );
        println!("{}", Table::new(pretty_ph_rows));

        println!("\nInner / recovery workflow detail (percent of local parent wall time)");
        println!(
            "Every canonical timing stage is listed, including zero values. The phase-control envelope contains its nonlinear solve when active; rows are inclusive and intentionally are not summed."
        );
        println!("{}", Table::new(pretty_inner_rows));

        println!("\nNonlinear envelope versus backend attempts");
        println!(
            "The parent envelope is compared with all attempt wall times, including wall-only failures. `known preparation` is the explicit active-set/RST/projection stage inside that parent; the remainder is runner/orchestration work still worth profiling. `no metrics` means internal counters are absent, while `unknown wall` means elapsed time itself is absent."
        );
        println!("{}", Table::new(pretty_nonlinear_rows));

        println!(
            "\nLegacy/ RST backend attempt detail (one row per attempt, including fallback and recovery attempts)"
        );
        println!(
            "Attempt percentages use that attempt's own backend wall time. `wall-only` preserves failed-call elapsed time without inventing iterations or callback splits; `n/a` means only that internal detail was not published."
        );
        println!("{}", Table::new(pretty_attempt_rows));

        println!("\nAccepted solver callback summary (percent of accepted backend wall time)");
        println!(
            "Residual/Jacobian/linear timings are measured inside legacy LM/NR/TR callbacks and linear solves. RST reports n/a for its internal linear split."
        );
        println!("{}", Table::new(pretty_solver_rows));

        println!("\nLegacy P,H end-to-end envelope");
        println!(
            "This chain uses explicit parent denominators: P,H workflow/facade, trial total/workflow, and nonlinear solve/trial (or monolithic execution). Legacy backend wall is shown separately because nested routes can contain several trial/attempt scopes. `n/a (monolithic)` is intentional."
        );
        println!("{}", Table::new(pretty_legacy_ph_rows));
    }
}

fn criterion_benchmark(c: &mut Criterion) {
    let targets = [
        (
            Universe::SmallH2O,
            reference_enthalpy_for(Universe::SmallH2O, false),
        ),
        (
            Universe::Stanjan15,
            reference_enthalpy_for(Universe::Stanjan15, true),
        ),
        (
            Universe::Tp1907,
            reference_enthalpy_for(Universe::Tp1907, true),
        ),
        (
            Universe::Rp1311_40,
            reference_enthalpy_for(Universe::Rp1311_40, true),
        ),
    ];
    print_baseline_table(&targets);

    let mut group = c.benchmark_group("equilibrium_facade");
    // P,H ranges are deliberately part of the baseline and are much more
    // expensive than a point. Ten samples keep this release target practical;
    // the internal table above remains the stable fixed-repeat comparison.
    group.sample_size(10);
    group.warm_up_time(Duration::from_secs(1));
    group.measurement_time(Duration::from_secs(2));
    for scenario in baseline_scenarios() {
        group.bench_function(BenchmarkId::from_parameter(scenario.label()), |bencher| {
            bencher.iter(|| {
                black_box(solve_scenario(
                    scenario,
                    target_enthalpy(scenario.universe, &targets),
                ))
            });
        });
    }
    group.finish();
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
