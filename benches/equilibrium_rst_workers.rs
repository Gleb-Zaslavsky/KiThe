//! Historical runtime experiment for the RST symbolic callback scheduler.
//!
//! This benchmark predates the prepared adapter's Sequential production
//! default. The current adapter owns an explicit Lambdify execution policy;
//! changing `RAYON_NUM_THREADS` alone does not turn the production callback
//! path into Parallel mode. The benchmark is retained for the old
//! request-level worker-scaling evidence, while the current policy crossover
//! is measured by `equilibrium_rst_policy_crossover.rs`.
//!
//! The historical experiment keeps the facade, equations, seeds, and solver
//! policy fixed and compares separate processes configured with
//! `RAYON_NUM_THREADS=1/2/4` (or the default pool). A separate process is
//! required because Rayon initializes its global pool only once.
//!
//! Run the diagnostic report in an optimized profile, one worker setting per
//! process:
//!
//! ```text
//! $env:RAYON_NUM_THREADS="1"
//! cargo bench --bench equilibrium_rst_workers --no-default-features --features equilibrium-instrumentation -- --noplot
//! ```
//!
//! Criterion still measures the same public facade calls. The printed table
//! is the evidence table: facade/workflow/backend denominators are separate,
//! and RST callback time is kept separate from its remaining solver overhead.

use std::env;
use std::time::{Duration, Instant};

use criterion::{BenchmarkId, Criterion, black_box, criterion_group, criterion_main};
use tabled::{Table, Tabled};

use KiThe::Thermodynamics::ChemEquilibrium::prelude::*;

mod support;

use support::{
    PRESSURE_PA, REFERENCE_PRESSURE_PA, REFERENCE_TEMPERATURE_K, RP1311_40, STANJAN_15,
    TP1907_GAS_17, tp1907_element_seed,
};

const REPORT_REPETITIONS: usize = 3;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Case {
    Stanjan15Ph,
    Tp1907Pt,
    Rp1311Pt,
}

impl Case {
    fn label(self) -> &'static str {
        match self {
            Self::Stanjan15Ph => "stanjan15_ph_element_inventory_monolithic",
            Self::Tp1907Pt => "tp1907_18_pt_element_inventory",
            Self::Rp1311Pt => "rp1311_40_pt_element_inventory",
        }
    }

    fn temperature(self) -> f64 {
        match self {
            Self::Stanjan15Ph => REFERENCE_TEMPERATURE_K,
            Self::Tp1907Pt => 700.0,
            Self::Rp1311Pt => 2_400.0,
        }
    }

    fn reference_pressure(self) -> f64 {
        match self {
            Self::Tp1907Pt => PRESSURE_PA,
            _ => REFERENCE_PRESSURE_PA,
        }
    }

    fn inventory(self) -> ElementInventory {
        let values = match self {
            Self::Stanjan15Ph => vec![("C", 6.0), ("H", 16.0), ("O", 20.0), ("N", 75.2)],
            Self::Tp1907Pt => vec![
                ("Ar", 0.05364816),
                ("C", 1.00182736),
                ("H", 2.0),
                ("N", 8.9462088),
                ("O", 2.40365472),
            ],
            Self::Rp1311Pt => vec![("Ar", 0.1), ("C", 0.1), ("H", 2.2), ("O", 4.3), ("N", 6.0)],
        };
        ElementInventory::from_amounts(values).expect("benchmark inventory must be valid")
    }

    fn base_builder(self) -> EquilibriumCalculatorBuilder {
        match self {
            Self::Stanjan15Ph => EquilibriumCalculator::builder()
                .ideal_gas_phase("gas", STANJAN_15)
                .expect("STANJAN phase"),
            Self::Tp1907Pt => EquilibriumCalculator::builder()
                .ideal_gas_phase("gas", TP1907_GAS_17)
                .expect("TP-1907 gas phase")
                .pure_condensed_phase("graphite", "C(gr)", PhysicalState::Solid)
                .expect("graphite phase"),
            Self::Rp1311Pt => EquilibriumCalculator::builder()
                .single_ideal_gas(RP1311_40)
                .expect("RP-1311 phase"),
        }
    }

    fn builder(
        self,
        solver: RustedSciTheSolver,
        target_enthalpy: f64,
    ) -> EquilibriumCalculatorBuilder {
        let builder = self
            .base_builder()
            .element_inventory(self.inventory())
            .prefer_libraries(["NASA_gas", "NASA_cond"])
            .offline_only()
            .pressure_pa(PRESSURE_PA)
            .reference_pressure_pa(self.reference_pressure())
            .solver_policy(SolverPolicy::Single(SolverBackend::RustedSciThe(solver)))
            .expect("single RST policy must validate")
            .preset(EquilibriumCalculatorPreset::Timed)
            .instrumentation(EquilibriumInstrumentationMode::Enabled);

        match self {
            Self::Stanjan15Ph => builder
                .ph_solve_mode(PhSolveMode::Monolithic)
                .at_total_enthalpy(
                    TotalEnthalpyJoules::new(target_enthalpy).expect("finite target"),
                    self.temperature(),
                    TemperatureBounds::new(2_000.0, 2_800.0).expect("valid bounds"),
                ),
            Self::Tp1907Pt => builder
                .element_numerical_seed(tp1907_element_seed())
                .at_temperature(self.temperature()),
            Self::Rp1311Pt => builder.at_temperature(self.temperature()),
        }
    }
}

#[derive(Clone, Copy, Debug)]
enum RstSolver {
    Lm,
    BacktrackingLm,
    MinpackLm,
    NielsenLm,
    DampedNewton,
}

impl RstSolver {
    fn label(self) -> &'static str {
        match self {
            Self::Lm => "rst_lm",
            Self::BacktrackingLm => "rst_backtracking_lm",
            Self::MinpackLm => "rst_minpack_lm",
            Self::NielsenLm => "rst_nielsen_lm",
            Self::DampedNewton => "rst_damped_newton",
        }
    }

    fn backend(self) -> RustedSciTheSolver {
        match self {
            Self::Lm => RustedSciTheSolver::LevenbergMarquardt,
            Self::BacktrackingLm => RustedSciTheSolver::BacktrackingLevenbergMarquardt,
            Self::MinpackLm => RustedSciTheSolver::MinpackLevenbergMarquardt,
            Self::NielsenLm => RustedSciTheSolver::NielsenLevenbergMarquardt,
            Self::DampedNewton => RustedSciTheSolver::DampedNewton,
        }
    }
}

fn reference_enthalpy() -> f64 {
    let point = Case::Stanjan15Ph
        .base_builder()
        .element_inventory(Case::Stanjan15Ph.inventory())
        .at_temperature(Case::Stanjan15Ph.temperature())
        .prefer_libraries(["NASA_gas", "NASA_cond"])
        .offline_only()
        .pressure_pa(PRESSURE_PA)
        .reference_pressure_pa(REFERENCE_PRESSURE_PA)
        .production_cascade()
        .preset(EquilibriumCalculatorPreset::Silent)
        .solve()
        .expect("reference P,T solve");
    let EquilibriumCalculatorOutcome::PtPoint(point) = point else {
        panic!("reference must be a P,T point")
    };
    let thermo = ResolvedThermochemistry::from_resolved_system(point.resolved())
        .expect("reference thermochemistry");
    thermo
        .enthalpy_model()
        .evaluate_total(
            point.solution().component_moles(),
            Case::Stanjan15Ph.temperature(),
        )
        .expect("reference enthalpy")
}

#[derive(Clone, Debug)]
struct AttemptMetricsRow {
    index: usize,
    backend: String,
    outcome: String,
    backend_us: f64,
    residual_calls: usize,
    residual_us: f64,
    jacobian_calls: usize,
    jacobian_us: f64,
    overhead_us: f64,
    state_residual_calls: Option<usize>,
    trial_residual_calls: Option<usize>,
    state_jacobian_calls: Option<usize>,
    trial_jacobian_calls: Option<usize>,
    factorization_us: Option<f64>,
    system_solve_us: Option<f64>,
    accepted_steps: Option<usize>,
    rejected_steps: Option<usize>,
    trace_len: Option<usize>,
}

#[derive(Clone, Debug)]
struct StageMetricsRow {
    scope: String,
    stage: String,
    time_us: f64,
    scope_us: f64,
}

#[derive(Clone, Debug)]
struct Evidence {
    case: &'static str,
    solver: &'static str,
    workers: usize,
    requested_workers: String,
    status: String,
    facade_us: f64,
    workflow_us: f64,
    attempts: Vec<AttemptMetricsRow>,
    stages: Vec<StageMetricsRow>,
}

fn outcome_label(outcome: &SolverAttemptOutcome) -> &'static str {
    match outcome {
        SolverAttemptOutcome::Accepted => "accepted",
        SolverAttemptOutcome::Failed { .. } => "failed",
        SolverAttemptOutcome::RejectedCandidate { .. } => "rejected",
        SolverAttemptOutcome::Skipped { .. } => "skipped",
    }
}

fn attempt_metrics(attempts: &[SolverAttemptReport]) -> Vec<AttemptMetricsRow> {
    attempts
        .iter()
        .enumerate()
        .filter_map(|(index, attempt)| {
            let metrics = attempt.metrics.as_ref()?;
            let timing = metrics.evaluation_timing.as_ref()?;
            Some(AttemptMetricsRow {
                index,
                backend: format!("{:?}", attempt.backend),
                outcome: outcome_label(&attempt.outcome).to_owned(),
                backend_us: metrics.elapsed_micros as f64,
                residual_calls: metrics.residual_evaluations,
                residual_us: timing.residual_evaluation_micros as f64,
                jacobian_calls: metrics.jacobian_evaluations,
                jacobian_us: timing.jacobian_evaluation_micros as f64,
                overhead_us: timing.solver_overhead_micros as f64,
                state_residual_calls: timing.state_residual_evaluations,
                trial_residual_calls: timing.trial_residual_evaluations,
                state_jacobian_calls: timing.state_jacobian_evaluations,
                trial_jacobian_calls: timing.trial_jacobian_evaluations,
                factorization_us: timing.linear_factorization_micros.map(|value| value as f64),
                system_solve_us: timing.linear_system_solve_micros.map(|value| value as f64),
                accepted_steps: timing.accepted_steps,
                rejected_steps: timing.rejected_steps,
                trace_len: timing.iteration_telemetry.as_ref().map(Vec::len),
            })
        })
        .collect()
}

fn stage_metrics(
    stages: &mut Vec<StageMetricsRow>,
    scope: &str,
    name: &str,
    duration: Duration,
    scope_us: f64,
) {
    stages.push(StageMetricsRow {
        scope: scope.to_owned(),
        stage: name.to_owned(),
        time_us: duration.as_secs_f64() * 1.0e6,
        scope_us,
    });
}

fn add_equilibrium_stages(
    stages: &mut Vec<StageMetricsRow>,
    scope: &str,
    timing: EquilibriumTimingReport,
) {
    let scope_us = timing.total().as_secs_f64() * 1.0e6;
    stage_metrics(
        stages,
        scope,
        "repository lookup",
        timing.repository_lookup(),
        scope_us,
    );
    stage_metrics(
        stages,
        scope,
        "thermochemistry preparation",
        timing.thermochemistry_preparation(),
        scope_us,
    );
    stage_metrics(
        stages,
        scope,
        "numeric closure construction",
        timing.numeric_closure_construction(),
        scope_us,
    );
    stage_metrics(
        stages,
        scope,
        "symbolic construction",
        timing.symbolic_construction(),
        scope_us,
    );
    stage_metrics(
        stages,
        scope,
        "equation construction",
        timing.equation_construction(),
        scope_us,
    );
    stage_metrics(
        stages,
        scope,
        "numerical problem preparation",
        timing.numerical_problem_preparation(),
        scope_us,
    );
    stage_metrics(
        stages,
        scope,
        "active-set formulation build",
        timing.active_set_formulation_build(),
        scope_us,
    );
    stage_metrics(
        stages,
        scope,
        "RST preparation",
        timing.rst_preparation(),
        scope_us,
    );
    stage_metrics(
        stages,
        scope,
        "projection build",
        timing.projection_build(),
        scope_us,
    );
    stage_metrics(
        stages,
        scope,
        "nonlinear solve",
        timing.nonlinear_solve(),
        scope_us,
    );
    stage_metrics(
        stages,
        scope,
        "phase-control envelope",
        timing.phase_control(),
        scope_us,
    );
    stage_metrics(stages, scope, "validation", timing.validation(), scope_us);
    stage_metrics(
        stages,
        scope,
        "postprocessing",
        timing.postprocessing(),
        scope_us,
    );
}

fn add_facade_stages(stages: &mut Vec<StageMetricsRow>, timing: EquilibriumTimingReport) {
    let scope_us = timing.total().as_secs_f64() * 1.0e6;
    let scope = "facade";
    stage_metrics(
        stages,
        scope,
        "resolution",
        timing.facade_resolution(),
        scope_us,
    );
    stage_metrics(
        stages,
        scope,
        "request preparation",
        timing.request_preparation(),
        scope_us,
    );
    stage_metrics(
        stages,
        scope,
        "route selection",
        timing.route_selection(),
        scope_us,
    );
    stage_metrics(
        stages,
        scope,
        "normalization",
        timing.normalization(),
        scope_us,
    );
    stage_metrics(
        stages,
        scope,
        "phase control",
        timing.phase_control(),
        scope_us,
    );
    stage_metrics(
        stages,
        scope,
        "thermochemistry preparation",
        timing.thermochemistry_preparation(),
        scope_us,
    );
    stage_metrics(stages, scope, "recovery", timing.recovery(), scope_us);
    stage_metrics(
        stages,
        scope,
        "result assembly",
        timing.result_assembly(),
        scope_us,
    );
}

fn add_ph_stages(stages: &mut Vec<StageMetricsRow>, timing: PhTemperatureTimingReport) {
    let scope_us = timing.total().as_secs_f64() * 1.0e6;
    let scope = "P,H transaction";
    stage_metrics(
        stages,
        scope,
        "scalar orchestration",
        timing.scalar_orchestration(),
        scope_us,
    );
    stage_metrics(
        stages,
        scope,
        "enthalpy evaluation",
        timing.enthalpy_evaluation(),
        scope_us,
    );
    stage_metrics(stages, scope, "trial total", timing.trial_total(), scope_us);
    stage_metrics(
        stages,
        scope,
        "trial inner execution",
        timing.trial_inner_execution(),
        scope_us,
    );
    stage_metrics(
        stages,
        scope,
        "trial preparation",
        timing.trial_preparation(),
        scope_us,
    );
    stage_metrics(
        stages,
        scope,
        "trial template preparation",
        timing.trial_template_preparation(),
        scope_us,
    );
    stage_metrics(
        stages,
        scope,
        "trial postprocessing",
        timing.trial_postprocessing(),
        scope_us,
    );
    stage_metrics(
        stages,
        scope,
        "bracket orchestration",
        timing.bracket_orchestration(),
        scope_us,
    );
    stage_metrics(
        stages,
        scope,
        "monolithic preparation",
        timing.monolithic_preparation(),
        scope_us,
    );
    stage_metrics(
        stages,
        scope,
        "request validation",
        timing.request_validation(),
        scope_us,
    );
    stage_metrics(
        stages,
        scope,
        "initial enthalpy",
        timing.initial_enthalpy(),
        scope_us,
    );
    stage_metrics(
        stages,
        scope,
        "prepared problem retarget",
        timing.prepared_problem_retarget(),
        scope_us,
    );
    stage_metrics(
        stages,
        scope,
        "Gibbs snapshot",
        timing.gibbs_snapshot(),
        scope_us,
    );
    stage_metrics(
        stages,
        scope,
        "monolithic RST preparation",
        timing.monolithic_rst_preparation(),
        scope_us,
    );
    stage_metrics(
        stages,
        scope,
        "monolithic backend execution",
        timing.monolithic_backend_execution(),
        scope_us,
    );
    stage_metrics(
        stages,
        scope,
        "monolithic execution",
        timing.monolithic_execution(),
        scope_us,
    );
    stage_metrics(
        stages,
        scope,
        "monolithic publication",
        timing.monolithic_publication(),
        scope_us,
    );
}

fn add_rst_preparation_stages(
    stages: &mut Vec<StageMetricsRow>,
    scope: &str,
    timing: RstPreparationTimingReport,
) {
    let scope_us = timing.total().as_secs_f64() * 1.0e6;
    stage_metrics(
        stages,
        scope,
        "equation generation",
        timing.equation_generation(),
        scope_us,
    );
    stage_metrics(
        stages,
        scope,
        "temperature substitution",
        timing.temperature_substitution(),
        scope_us,
    );
    stage_metrics(stages, scope, "row scaling", timing.row_scaling(), scope_us);
    stage_metrics(
        stages,
        scope,
        "enthalpy equation",
        timing.enthalpy_equation(),
        scope_us,
    );
    stage_metrics(
        stages,
        scope,
        "variable materialization",
        timing.variable_materialization(),
        scope_us,
    );
    stage_metrics(
        stages,
        scope,
        "parameter evaluation",
        timing.parameter_evaluation(),
        scope_us,
    );
    stage_metrics(
        stages,
        scope,
        "log-mole bounds",
        timing.log_mole_bounds(),
        scope_us,
    );
    stage_metrics(
        stages,
        scope,
        if timing.backend_stage_breakdown_available() {
            "RST preparation aggregate"
        } else {
            "RST backend preparation (aggregate; internal split unavailable)"
        },
        timing.problem_materialization(),
        scope_us,
    );
    if timing.backend_stage_breakdown_available() {
        for (label, duration) in [
            ("RST input validation", timing.rst_input_validation()),
            ("RST expression parsing", timing.rst_expression_parsing()),
            (
                "RST expression graph materialization",
                timing.rst_expression_graph_materialization(),
            ),
            (
                "RST symbolic differentiation",
                Some(timing.symbolic_differentiation()),
            ),
            (
                "RST residual lambdification",
                Some(timing.residual_lambdification()),
            ),
            (
                "RST Jacobian lambdification",
                Some(timing.jacobian_lambdification()),
            ),
            ("RST parameter binding", timing.rst_parameter_binding()),
            (
                "RST prepared problem assembly",
                timing.rst_problem_assembly(),
            ),
            ("RST AOT materialization", timing.rst_aot_materialization()),
            ("RST final validation", timing.rst_final_validation()),
            ("RST unattributed preparation", timing.rst_unattributed()),
        ] {
            if let Some(duration) = duration {
                stage_metrics(stages, scope, label, duration, scope_us);
            }
        }
    }
}

fn solve_case(
    case: Case,
    solver: RstSolver,
    target_enthalpy: f64,
) -> Result<EquilibriumCalculatorOutcome, EquilibriumCalculatorError> {
    case.builder(solver.backend(), target_enthalpy).solve()
}

fn collect_case(case: Case, solver: RstSolver, target_enthalpy: f64) -> Evidence {
    let requested_workers = env::var("RAYON_NUM_THREADS").unwrap_or_else(|_| "default".to_owned());
    let workers = rayon::current_num_threads();
    let started = Instant::now();
    let outcome = solve_case(case, solver, target_enthalpy);
    let facade_us = started.elapsed().as_secs_f64() * 1.0e6;
    match outcome {
        Ok(EquilibriumCalculatorOutcome::PtPoint(point)) => Evidence {
            case: case.label(),
            solver: solver.label(),
            workers,
            requested_workers,
            status: "accepted".to_owned(),
            facade_us,
            workflow_us: point.solution().timing_report().total().as_secs_f64() * 1.0e6,
            attempts: attempt_metrics(&point.solution().solve_report().attempts),
            stages: {
                let mut stages = Vec::new();
                add_facade_stages(&mut stages, *point.facade_timing());
                add_equilibrium_stages(
                    &mut stages,
                    "P,T workflow",
                    *point.solution().timing_report(),
                );
                add_rst_preparation_stages(
                    &mut stages,
                    "P,T RST preparation",
                    point.solution().timing_report().rst_preparation_timing(),
                );
                stages
            },
        },
        Ok(EquilibriumCalculatorOutcome::PhPoint(point)) => {
            let report = point.solution().report();
            Evidence {
                case: case.label(),
                solver: solver.label(),
                workers,
                requested_workers,
                status: format!("accepted:{:?}", report.solve_path()),
                facade_us,
                workflow_us: report.timing().total().as_secs_f64() * 1.0e6,
                attempts: report
                    .monolithic_evidence()
                    .map(|evidence| attempt_metrics(&evidence.solve_report().attempts))
                    .unwrap_or_default(),
                stages: {
                    let mut stages = Vec::new();
                    add_facade_stages(&mut stages, *point.facade_timing());
                    add_ph_stages(&mut stages, report.timing());
                    add_equilibrium_stages(
                        &mut stages,
                        "P,H inner preparation",
                        report.inner_timing(),
                    );
                    add_rst_preparation_stages(
                        &mut stages,
                        "P,H RST preparation",
                        report.timing().rst_preparation_timing(),
                    );
                    stages
                },
            }
        }
        Ok(other) => Evidence {
            case: case.label(),
            solver: solver.label(),
            workers,
            requested_workers,
            status: format!("unexpected:{other:?}"),
            facade_us,
            workflow_us: 0.0,
            attempts: Vec::new(),
            stages: Vec::new(),
        },
        Err(error) => {
            let attempts = match &error {
                EquilibriumCalculatorError::Preparation(cause)
                | EquilibriumCalculatorError::Solve(cause) => {
                    attempt_metrics(&cause.backend_attempt_reports())
                }
                _ => Vec::new(),
            };
            Evidence {
                case: case.label(),
                solver: solver.label(),
                workers,
                requested_workers,
                status: "failed".to_owned(),
                facade_us,
                workflow_us: 0.0,
                attempts,
                stages: Vec::new(),
            }
        }
    }
}

fn ms(value_us: f64) -> String {
    format!("{:.3}", value_us / 1_000.0)
}

fn share(value_us: f64, denominator_us: f64) -> String {
    if denominator_us <= f64::EPSILON {
        "n/a".to_owned()
    } else {
        format!("{:.1}%", 100.0 * value_us.max(0.0) / denominator_us)
    }
}

#[derive(Tabled)]
struct Row {
    #[tabled(rename = "case")]
    case: String,
    #[tabled(rename = "solver")]
    solver: String,
    #[tabled(rename = "requested")]
    requested_workers: String,
    #[tabled(rename = "actual workers")]
    workers: usize,
    #[tabled(rename = "status")]
    status: String,
    #[tabled(rename = "facade ms")]
    facade: String,
    #[tabled(rename = "workflow ms")]
    workflow: String,
    #[tabled(rename = "attempts")]
    attempts: usize,
}

#[derive(Tabled)]
struct AttemptRow {
    #[tabled(rename = "case")]
    case: String,
    #[tabled(rename = "solver")]
    solver: String,
    #[tabled(rename = "#")]
    index: usize,
    #[tabled(rename = "backend")]
    backend: String,
    #[tabled(rename = "outcome")]
    outcome: String,
    #[tabled(rename = "backend ms")]
    backend_ms: String,
    #[tabled(rename = "residual ms / %")]
    residual: String,
    #[tabled(rename = "Jacobian ms / %")]
    jacobian: String,
    #[tabled(rename = "overhead ms / %")]
    overhead: String,
    #[tabled(rename = "res calls")]
    residual_calls: usize,
    #[tabled(rename = "Jac calls")]
    jacobian_calls: usize,
}

#[derive(Tabled)]
struct NativeAttemptRow {
    #[tabled(rename = "case")]
    case: String,
    #[tabled(rename = "solver")]
    solver: String,
    #[tabled(rename = "#")]
    index: usize,
    #[tabled(rename = "backend")]
    backend: String,
    #[tabled(rename = "state R/J")]
    state_callbacks: String,
    #[tabled(rename = "trial R/J")]
    trial_callbacks: String,
    #[tabled(rename = "factor ms")]
    factorization_ms: String,
    #[tabled(rename = "system ms")]
    system_solve_ms: String,
    #[tabled(rename = "accepted/rejected")]
    steps: String,
    #[tabled(rename = "trace rows")]
    trace_len: String,
}

#[derive(Tabled)]
struct StageRow {
    #[tabled(rename = "case")]
    case: String,
    #[tabled(rename = "solver")]
    solver: String,
    #[tabled(rename = "requested")]
    requested_workers: String,
    #[tabled(rename = "workers")]
    workers: usize,
    #[tabled(rename = "scope")]
    scope: String,
    #[tabled(rename = "stage")]
    stage: String,
    #[tabled(rename = "stage ms")]
    stage_ms: String,
    #[tabled(rename = "scope ms")]
    scope_ms: String,
    #[tabled(rename = "% of scope")]
    scope_share: String,
}

fn print_report(rows: &[Evidence]) {
    let table_rows = rows
        .iter()
        .map(|evidence| Row {
            case: evidence.case.to_owned(),
            solver: evidence.solver.to_owned(),
            requested_workers: evidence.requested_workers.clone(),
            workers: evidence.workers,
            status: evidence.status.clone(),
            facade: ms(evidence.facade_us),
            workflow: ms(evidence.workflow_us),
            attempts: evidence.attempts.len(),
        })
        .collect::<Vec<_>>();
    let attempt_rows = rows
        .iter()
        .flat_map(|evidence| {
            evidence.attempts.iter().map(|metrics| AttemptRow {
                case: evidence.case.to_owned(),
                solver: evidence.solver.to_owned(),
                index: metrics.index,
                backend: metrics.backend.clone(),
                outcome: metrics.outcome.clone(),
                backend_ms: ms(metrics.backend_us),
                residual: format!(
                    "{} / {}",
                    ms(metrics.residual_us),
                    share(metrics.residual_us, metrics.backend_us)
                ),
                jacobian: format!(
                    "{} / {}",
                    ms(metrics.jacobian_us),
                    share(metrics.jacobian_us, metrics.backend_us)
                ),
                overhead: format!(
                    "{} / {}",
                    ms(metrics.overhead_us),
                    share(metrics.overhead_us, metrics.backend_us)
                ),
                residual_calls: metrics.residual_calls,
                jacobian_calls: metrics.jacobian_calls,
            })
        })
        .collect::<Vec<_>>();
    let stage_rows = rows
        .iter()
        .flat_map(|evidence| {
            evidence.stages.iter().map(|stage| StageRow {
                case: evidence.case.to_owned(),
                solver: evidence.solver.to_owned(),
                requested_workers: evidence.requested_workers.clone(),
                workers: evidence.workers,
                scope: stage.scope.clone(),
                stage: stage.stage.clone(),
                stage_ms: ms(stage.time_us),
                scope_ms: ms(stage.scope_us),
                scope_share: share(stage.time_us, stage.scope_us),
            })
        })
        .collect::<Vec<_>>();
    let native_attempt_rows = rows
        .iter()
        .flat_map(|evidence| {
            evidence.attempts.iter().filter_map(|metrics| {
                if metrics.state_residual_calls.is_none() && metrics.trial_residual_calls.is_none()
                {
                    return None;
                }
                Some(NativeAttemptRow {
                    case: evidence.case.to_owned(),
                    solver: evidence.solver.to_owned(),
                    index: metrics.index,
                    backend: metrics.backend.clone(),
                    state_callbacks: format!(
                        "{}/{}",
                        metrics
                            .state_residual_calls
                            .map_or_else(|| "-".to_owned(), |value| value.to_string()),
                        metrics
                            .state_jacobian_calls
                            .map_or_else(|| "-".to_owned(), |value| value.to_string())
                    ),
                    trial_callbacks: format!(
                        "{}/{}",
                        metrics
                            .trial_residual_calls
                            .map_or_else(|| "-".to_owned(), |value| value.to_string()),
                        metrics
                            .trial_jacobian_calls
                            .map_or_else(|| "-".to_owned(), |value| value.to_string())
                    ),
                    factorization_ms: metrics.factorization_us.map_or_else(|| "-".to_owned(), ms),
                    system_solve_ms: metrics.system_solve_us.map_or_else(|| "-".to_owned(), ms),
                    steps: format!(
                        "{}/{}",
                        metrics
                            .accepted_steps
                            .map_or_else(|| "-".to_owned(), |value| value.to_string()),
                        metrics
                            .rejected_steps
                            .map_or_else(|| "-".to_owned(), |value| value.to_string())
                    ),
                    trace_len: metrics
                        .trace_len
                        .map_or_else(|| "-".to_owned(), |value| value.to_string()),
                })
            })
        })
        .collect::<Vec<_>>();
    println!(
        "\nPERF-RST-02 RST worker/scheduling experiment\n{}",
        Table::new(table_rows)
    );
    println!(
        "\nPERF-RST-02 worker-attributed stage timing\n{}",
        Table::new(stage_rows)
    );
    println!(
        "\nPERF-RST-02 every RST attempt\n{}",
        Table::new(attempt_rows)
    );
    println!(
        "\nPERF-RST-02 RST native solver telemetry\n{}",
        Table::new(native_attempt_rows)
    );
    println!(
        "PERF-RST-02 report. Stage percentages use the local scope total shown in the row. All stage timings are inclusive: nonlinear solve contains runner work. RST 0.4.19 preparation children are printed in the stage table, while native solver counters are kept in the separate telemetry table. RST callback percentages use backend wall as denominator; run each worker setting in a fresh process. Use PERF-RST-06 for the current Sequential/Parallel policy crossover."
    );
}

fn criterion_benchmark(c: &mut Criterion) {
    let target_enthalpy = reference_enthalpy();
    let cases = [Case::Stanjan15Ph, Case::Tp1907Pt, Case::Rp1311Pt];
    let solvers = [
        RstSolver::Lm,
        RstSolver::BacktrackingLm,
        RstSolver::MinpackLm,
        RstSolver::NielsenLm,
        RstSolver::DampedNewton,
    ];
    let filter = env::var("RST_WORKERS_FILTER").ok();
    let report_repetitions = env::var("RST_WORKERS_REPORT_REPETITIONS")
        .ok()
        .and_then(|value| value.parse::<usize>().ok())
        .unwrap_or(REPORT_REPETITIONS);
    let selected = |case: Case, solver: RstSolver| {
        filter.as_ref().is_none_or(|filter| {
            let label = format!("{}_{}", case.label(), solver.label());
            label.contains(filter)
        })
    };
    let mut evidence = Vec::new();
    for case in cases {
        for solver in solvers {
            if !selected(case, solver) {
                continue;
            }
            for _ in 0..report_repetitions {
                evidence.push(collect_case(case, solver, target_enthalpy));
            }
        }
    }
    print_report(&evidence);

    let mut group = c.benchmark_group("equilibrium_rst_workers");
    group.sample_size(10);
    group.warm_up_time(Duration::from_secs(1));
    group.measurement_time(Duration::from_secs(2));
    for case in cases {
        for solver in solvers {
            if !selected(case, solver) {
                continue;
            }
            let id = format!("{}_{}", case.label(), solver.label());
            group.bench_function(BenchmarkId::from_parameter(id), |bencher| {
                bencher.iter(|| black_box(solve_case(case, solver, target_enthalpy)));
            });
        }
    }
    group.finish();
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
