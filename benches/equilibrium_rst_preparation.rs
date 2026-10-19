//! PERF-RST-01: detailed preparation/execution baseline for the RST backend.
//!
//! This benchmark deliberately uses one fixed RST backend and never enables a
//! legacy fallback. Criterion measures the solve itself; the tables printed
//! before Criterion starts are diagnostic evidence from the same public
//! `EquilibriumCalculator` facade.
//!
//! Run in an optimized profile with:
//!
//! ```text
//! cargo bench --bench equilibrium_rst_preparation --features equilibrium-instrumentation -- --noplot
//! ```
//!
//! Set `RST_PREPARATION_REPORT_ONLY=1` to collect and print the complete
//! diagnostic matrix without entering the slower Criterion phase.
//! Set `RST_PREPARATION_COMPACT=1` together with it to print one median row
//! per route and solver instead of the full repeated stage/attempt tables.
//!
//! The macro table contains only the facade/workflow partition. The following
//! tables use their own scope denominators and list concrete operations. Stage
//! values are inclusive when a stage contains nested work and must not be
//! summed as if they were mutually exclusive.

use std::time::{Duration, Instant};

use criterion::{BenchmarkId, Criterion, black_box, criterion_group, criterion_main};
use tabled::{Table, Tabled};

use KiThe::Thermodynamics::ChemEquilibrium::prelude::*;

const REPORT_REPETITIONS: usize = 3;
const RANGE_TEMPERATURES: [f64; 3] = [2_400.0, 2_500.0, 2_600.0];

mod support;

#[cfg(feature = "equilibrium-benchmark-rst-frontends")]
use RustedSciThe::numerical::Nonlinear_systems::prelude::SymbolicLambdifyFrontend;
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
            Self::Small3 | Self::Stanjan15 => REFERENCE_TEMPERATURE_K,
            Self::Tp1907_18 => 700.0,
            Self::Rp1311_40 => 2_400.0,
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
    PtElementInventoryRange,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum RstSolver {
    Lm,
    BacktrackingLm,
    MinpackLm,
    NielsenLm,
    // RST's Newton/NR path; the exposed method is line-search damped.
    DampedNewton,
}

impl RstSolver {
    fn label(self) -> &'static str {
        match self {
            Self::Lm => "rst_levenberg_marquardt",
            Self::BacktrackingLm => "rst_backtracking_levenberg_marquardt",
            Self::MinpackLm => "rst_minpack_levenberg_marquardt",
            Self::NielsenLm => "rst_nielsen_levenberg_marquardt",
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

impl Task {
    fn label(self) -> &'static str {
        match self {
            Self::PtMolecular => "pt_molecular",
            Self::PtElementInventory => "pt_element_inventory",
            Self::PhMolecular => "ph_molecular",
            Self::PhElementInventory => "ph_element_inventory",
            Self::PtElementInventoryRange => "pt_element_inventory_range",
        }
    }

    fn is_elemental(self) -> bool {
        matches!(
            self,
            Self::PtElementInventory | Self::PhElementInventory | Self::PtElementInventoryRange
        )
    }
}

#[derive(Clone, Copy, Debug)]
struct Case {
    universe: Universe,
    task: Task,
    rst_solver: RstSolver,
}

impl Case {
    fn label(self) -> String {
        format!(
            "{}_{}_{}",
            self.universe.label(),
            self.task.label(),
            self.rst_solver.label()
        )
    }
}

#[derive(Clone, Debug)]
struct StageSample {
    case: String,
    scope: String,
    stage: String,
    time_us: f64,
    scope_us: f64,
    semantics: String,
}

#[derive(Clone, Debug)]
struct AttemptSample {
    case: String,
    scope: String,
    index: usize,
    backend: String,
    outcome: String,
    wall_us: f64,
    iterations: usize,
    residual_calls: usize,
    jacobian_calls: usize,
    linear_solves: usize,
    residual_us: Option<f64>,
    jacobian_us: Option<f64>,
    linear_us: Option<f64>,
    overhead_us: Option<f64>,
}

#[derive(Clone, Debug)]
struct RstPreparationSample {
    case: String,
    scope: String,
    report: RstPreparationTimingReport,
}

#[derive(Clone, Debug, Default)]
struct Evidence {
    case: String,
    route: String,
    policy: String,
    status: String,
    wall_us: f64,
    workflow_us: f64,
    stages: Vec<StageSample>,
    attempts: Vec<AttemptSample>,
    rst_preparations: Vec<RstPreparationSample>,
    residual_evaluations: usize,
    jacobian_evaluations: usize,
    formulation_builds: usize,
    formulation_reuses: usize,
    symbolic_parameter_updates: usize,
    phase_transitions: usize,
    accepted_backend: String,
    iterations: usize,
}

fn micros(duration: Duration) -> f64 {
    duration.as_secs_f64() * 1.0e6
}

fn stage(
    rows: &mut Vec<StageSample>,
    case: &str,
    scope: &str,
    name: &str,
    duration: Duration,
    scope_us: f64,
    semantics: &str,
) {
    rows.push(StageSample {
        case: case.to_owned(),
        scope: scope.to_owned(),
        stage: name.to_owned(),
        time_us: micros(duration),
        scope_us,
        semantics: semantics.to_owned(),
    });
}

fn equilibrium_stages(
    rows: &mut Vec<StageSample>,
    case: &str,
    scope: &str,
    timing: EquilibriumTimingReport,
    scope_us: f64,
) {
    let inclusive = "inclusive timing; nested stages are not additive";
    stage(
        rows,
        case,
        scope,
        "repository lookup",
        timing.repository_lookup(),
        scope_us,
        inclusive,
    );
    stage(
        rows,
        case,
        scope,
        "thermochemistry preparation",
        timing.thermochemistry_preparation(),
        scope_us,
        inclusive,
    );
    stage(
        rows,
        case,
        scope,
        "numeric closure construction",
        timing.numeric_closure_construction(),
        scope_us,
        inclusive,
    );
    stage(
        rows,
        case,
        scope,
        "symbolic construction",
        timing.symbolic_construction(),
        scope_us,
        inclusive,
    );
    stage(
        rows,
        case,
        scope,
        "equation construction",
        timing.equation_construction(),
        scope_us,
        inclusive,
    );
    stage(
        rows,
        case,
        scope,
        "numerical problem preparation",
        timing.numerical_problem_preparation(),
        scope_us,
        inclusive,
    );
    stage(
        rows,
        case,
        scope,
        "active-set formulation build",
        timing.active_set_formulation_build(),
        scope_us,
        inclusive,
    );
    stage(
        rows,
        case,
        scope,
        "RST preparation",
        timing.rst_preparation(),
        scope_us,
        inclusive,
    );
    stage(
        rows,
        case,
        scope,
        "projection build",
        timing.projection_build(),
        scope_us,
        inclusive,
    );
    stage(
        rows,
        case,
        scope,
        "nonlinear solve",
        timing.nonlinear_solve(),
        scope_us,
        inclusive,
    );
    stage(
        rows,
        case,
        scope,
        "phase-control envelope",
        timing.phase_control(),
        scope_us,
        inclusive,
    );
    stage(
        rows,
        case,
        scope,
        "validation",
        timing.validation(),
        scope_us,
        inclusive,
    );
    stage(
        rows,
        case,
        scope,
        "postprocessing",
        timing.postprocessing(),
        scope_us,
        inclusive,
    );
}

fn rst_preparation_stages(
    rows: &mut Vec<StageSample>,
    case: &str,
    scope: &str,
    timing: RstPreparationTimingReport,
) {
    if timing.total() == Duration::ZERO {
        return;
    }
    let scope_us = micros(timing.total());
    let inclusive = "measured child stage; RST aggregate is the backend-call scope";
    stage(
        rows,
        case,
        scope,
        "KiThe symbolic input validation",
        timing.symbolic_input_validation(),
        scope_us,
        inclusive,
    );
    stage(
        rows,
        case,
        scope,
        "KiThe symbolic phase topology",
        timing.symbolic_phase_topology(),
        scope_us,
        inclusive,
    );
    stage(
        rows,
        case,
        scope,
        "KiThe symbolic variable expansion",
        timing.symbolic_variable_expansion(),
        scope_us,
        inclusive,
    );
    stage(
        rows,
        case,
        scope,
        "KiThe symbolic phase totals",
        timing.symbolic_phase_totals(),
        scope_us,
        inclusive,
    );
    stage(
        rows,
        case,
        scope,
        "KiThe symbolic activity offsets",
        timing.symbolic_activity_offsets(),
        scope_us,
        inclusive,
    );
    stage(
        rows,
        case,
        scope,
        "KiThe symbolic reaction rows",
        timing.symbolic_reaction_assembly(),
        scope_us,
        inclusive,
    );
    stage(
        rows,
        case,
        scope,
        "KiThe symbolic element-balance rows",
        timing.symbolic_element_balance_assembly(),
        scope_us,
        inclusive,
    );
    stage(
        rows,
        case,
        scope,
        "KiThe symbolic residual simplification",
        timing.symbolic_residual_simplification(),
        scope_us,
        inclusive,
    );
    stage(
        rows,
        case,
        scope,
        "equation generation",
        timing.equation_generation(),
        scope_us,
        inclusive,
    );
    stage(
        rows,
        case,
        scope,
        "temperature substitution",
        timing.temperature_substitution(),
        scope_us,
        inclusive,
    );
    stage(
        rows,
        case,
        scope,
        "row scaling",
        timing.row_scaling(),
        scope_us,
        inclusive,
    );
    stage(
        rows,
        case,
        scope,
        "enthalpy equation",
        timing.enthalpy_equation(),
        scope_us,
        inclusive,
    );
    stage(
        rows,
        case,
        scope,
        "variable materialization",
        timing.variable_materialization(),
        scope_us,
        inclusive,
    );
    stage(
        rows,
        case,
        scope,
        "parameter evaluation",
        timing.parameter_evaluation(),
        scope_us,
        inclusive,
    );
    stage(
        rows,
        case,
        scope,
        "log-mole bounds",
        timing.log_mole_bounds(),
        scope_us,
        inclusive,
    );
    stage(
        rows,
        case,
        scope,
        "RST backend preparation (aggregate)",
        timing.problem_materialization(),
        scope_us,
        "RST 0.4.19 aggregate; child stages are reported below",
    );
    if timing.backend_stage_breakdown_available() {
        stage(
            rows,
            case,
            scope,
            "RST symbolic differentiation",
            timing.symbolic_differentiation(),
            scope_us,
            "RST-reported child stage",
        );
        stage(
            rows,
            case,
            scope,
            "RST residual lambdification",
            timing.residual_lambdification(),
            scope_us,
            "RST-reported child stage",
        );
        stage(
            rows,
            case,
            scope,
            "RST Jacobian lambdification",
            timing.jacobian_lambdification(),
            scope_us,
            "RST-reported child stage",
        );
        for (label, duration) in [
            ("RST input validation", timing.rst_input_validation()),
            ("RST expression parsing", timing.rst_expression_parsing()),
            (
                "RST expression graph materialization",
                timing.rst_expression_graph_materialization(),
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
                stage(
                    rows,
                    case,
                    scope,
                    label,
                    duration,
                    scope_us,
                    "RST-reported child stage",
                );
            }
        }
    }
    if let Some(build_duration) = timing.backend_build_duration() {
        stage(
            rows,
            case,
            scope,
            "RST artifact build",
            build_duration,
            scope_us,
            "optional RST artifact lifecycle stage",
        );
    }
}

fn add_rst_preparation_sample(
    evidence: &mut Evidence,
    case: &str,
    scope: &str,
    timing: RstPreparationTimingReport,
) {
    if timing.total() > Duration::ZERO {
        evidence.rst_preparations.push(RstPreparationSample {
            case: case.to_owned(),
            scope: scope.to_owned(),
            report: timing,
        });
    }
}

fn facade_stages(
    rows: &mut Vec<StageSample>,
    case: &str,
    timing: EquilibriumTimingReport,
    scope_us: f64,
) {
    let semantics = "facade-owned stage; request/workflow timing is separate";
    stage(
        rows,
        case,
        "facade",
        "resolution",
        timing.facade_resolution(),
        scope_us,
        semantics,
    );
    stage(
        rows,
        case,
        "facade",
        "request preparation",
        timing.request_preparation(),
        scope_us,
        semantics,
    );
    stage(
        rows,
        case,
        "facade",
        "route selection",
        timing.route_selection(),
        scope_us,
        semantics,
    );
    stage(
        rows,
        case,
        "facade",
        "normalization",
        timing.normalization(),
        scope_us,
        semantics,
    );
    stage(
        rows,
        case,
        "facade",
        "phase control",
        timing.phase_control(),
        scope_us,
        semantics,
    );
    stage(
        rows,
        case,
        "facade",
        "thermochemistry preparation",
        timing.thermochemistry_preparation(),
        scope_us,
        semantics,
    );
    stage(
        rows,
        case,
        "facade",
        "recovery",
        timing.recovery(),
        scope_us,
        semantics,
    );
    stage(
        rows,
        case,
        "facade",
        "result assembly",
        timing.result_assembly(),
        scope_us,
        semantics,
    );
}

fn ph_stages(
    rows: &mut Vec<StageSample>,
    case: &str,
    timing: PhTemperatureTimingReport,
    scope_us: f64,
) {
    let semantics = "P,H route stage; nested inner timing is reported separately";
    stage(
        rows,
        case,
        "P,H transaction",
        "scalar orchestration",
        timing.scalar_orchestration(),
        scope_us,
        semantics,
    );
    stage(
        rows,
        case,
        "P,H transaction",
        "enthalpy evaluation",
        timing.enthalpy_evaluation(),
        scope_us,
        semantics,
    );
    stage(
        rows,
        case,
        "P,H transaction",
        "trial total",
        timing.trial_total(),
        scope_us,
        semantics,
    );
    stage(
        rows,
        case,
        "P,H transaction",
        "trial inner execution",
        timing.trial_inner_execution(),
        scope_us,
        semantics,
    );
    stage(
        rows,
        case,
        "P,H transaction",
        "trial preparation",
        timing.trial_preparation(),
        scope_us,
        semantics,
    );
    stage(
        rows,
        case,
        "P,H transaction",
        "trial template preparation",
        timing.trial_template_preparation(),
        scope_us,
        semantics,
    );
    stage(
        rows,
        case,
        "P,H transaction",
        "trial postprocessing",
        timing.trial_postprocessing(),
        scope_us,
        semantics,
    );
    stage(
        rows,
        case,
        "P,H transaction",
        "bracket orchestration",
        timing.bracket_orchestration(),
        scope_us,
        semantics,
    );
    stage(
        rows,
        case,
        "P,H transaction",
        "monolithic preparation",
        timing.monolithic_preparation(),
        scope_us,
        semantics,
    );
    stage(
        rows,
        case,
        "P,H transaction",
        "monolithic RST preparation",
        timing.monolithic_rst_preparation(),
        scope_us,
        semantics,
    );
    stage(
        rows,
        case,
        "P,H transaction",
        "monolithic backend execution",
        timing.monolithic_backend_execution(),
        scope_us,
        semantics,
    );
    stage(
        rows,
        case,
        "P,H transaction",
        "monolithic execution",
        timing.monolithic_execution(),
        scope_us,
        semantics,
    );
    stage(
        rows,
        case,
        "P,H transaction",
        "monolithic publication",
        timing.monolithic_publication(),
        scope_us,
        semantics,
    );
}

fn outcome_label(outcome: &SolverAttemptOutcome) -> String {
    match outcome {
        SolverAttemptOutcome::Accepted => "accepted".to_owned(),
        SolverAttemptOutcome::Failed { .. } => "failed".to_owned(),
        SolverAttemptOutcome::RejectedCandidate { .. } => "rejected_candidate".to_owned(),
        SolverAttemptOutcome::Skipped { .. } => "skipped".to_owned(),
    }
}

fn attempt_rows(
    rows: &mut Vec<AttemptSample>,
    case: &str,
    scope: &str,
    report: &EquilibriumSolveReport,
) {
    attempt_rows_slice(rows, case, scope, &report.attempts);
}

fn attempt_rows_slice(
    rows: &mut Vec<AttemptSample>,
    case: &str,
    scope: &str,
    attempts: &[SolverAttemptReport],
) {
    for (index, attempt) in attempts.iter().enumerate() {
        let metrics = attempt.metrics.as_ref();
        let timing = metrics.and_then(|metrics| metrics.evaluation_timing.as_ref());
        rows.push(AttemptSample {
            case: case.to_owned(),
            scope: scope.to_owned(),
            index,
            backend: format!("{:?}", attempt.backend),
            outcome: outcome_label(&attempt.outcome),
            wall_us: metrics.map_or(0.0, |metrics| metrics.elapsed_micros as f64),
            iterations: metrics.map_or(0, |metrics| metrics.iterations),
            residual_calls: metrics.map_or(0, |metrics| metrics.residual_evaluations),
            jacobian_calls: metrics.map_or(0, |metrics| metrics.jacobian_evaluations),
            linear_solves: metrics.map_or(0, |metrics| metrics.linear_solves),
            residual_us: timing.map(|timing| timing.residual_evaluation_micros as f64),
            jacobian_us: timing.map(|timing| timing.jacobian_evaluation_micros as f64),
            linear_us: timing
                .and_then(|timing| timing.linear_system_micros.map(|value| value as f64)),
            overhead_us: timing.map(|timing| timing.solver_overhead_micros as f64),
        });
    }
}

fn add_solution_evidence(
    evidence: &mut Evidence,
    solution: &MultiphaseEquilibriumSolution,
    scope: &str,
) {
    let report = solution.solve_report();
    summarize_solution(evidence, solution);
    if let Some(phase_control) = solution.phase_control_report() {
        for (iteration, nonlinear_report) in phase_control.nonlinear_reports.iter().enumerate() {
            attempt_rows(
                &mut evidence.attempts,
                &evidence.case,
                &format!("{scope}/active-set-{iteration}"),
                nonlinear_report,
            );
        }
    } else {
        attempt_rows(&mut evidence.attempts, &evidence.case, scope, report);
    }
}

fn summarize_solution(evidence: &mut Evidence, solution: &MultiphaseEquilibriumSolution) {
    let report = solution.solve_report();
    evidence.accepted_backend = format!("{:?}", report.accepted_backend);
    evidence.iterations += solution.nonlinear_iterations();
    evidence.phase_transitions += solution.phase_control_transitions();
}

fn inventory(universe: Universe) -> ElementInventory {
    let amounts = match universe {
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
    ElementInventory::from_amounts(amounts).expect("benchmark inventory must be valid")
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
        Universe::Tp1907_18 => base_builder(universe).initial_phase_moles([
            (component(Some("gas"), "Ar"), 0.05364816),
            (component(Some("gas"), "CH4"), 1.0),
            (component(Some("gas"), "N2"), 3.0),
            (component(Some("gas"), "O2"), 2.0),
            (component(Some("graphite"), "C(gr)"), 0.1),
        ]),
        Universe::Rp1311_40 => base_builder(universe).initial_phase_moles([
            (component(None, "Ar"), 0.1),
            (component(None, "CO2"), 0.1),
            (component(None, "H2"), 1.0),
            (component(None, "H2O"), 0.1),
            (component(None, "N2"), 3.0),
            (component(None, "O2"), 2.0),
        ]),
    }
}

fn rst_builder(
    builder: EquilibriumCalculatorBuilder,
    reference_pressure_pa: f64,
    rst_solver: RstSolver,
) -> EquilibriumCalculatorBuilder {
    builder
        .prefer_libraries(["NASA_gas", "NASA_cond"])
        .offline_only()
        .pressure_pa(PRESSURE_PA)
        .reference_pressure_pa(reference_pressure_pa)
        .solver_policy(SolverPolicy::Single(SolverBackend::RustedSciThe(
            rst_solver.backend(),
        )))
        .expect("single RST policy must validate")
        .preset(EquilibriumCalculatorPreset::Timed)
        .instrumentation(EquilibriumInstrumentationMode::Enabled)
}

fn builder_for(case: Case, target_enthalpy: f64) -> EquilibriumCalculatorBuilder {
    let builder = if case.task.is_elemental() {
        base_builder(case.universe).element_inventory(inventory(case.universe))
    } else {
        molecular_builder(case.universe)
    };
    let builder = if matches!(case.task, Task::PhElementInventory)
        && matches!(case.universe, Universe::Tp1907_18)
    {
        builder.element_numerical_seed(tp1907_element_seed())
    } else {
        builder
    };
    let builder = match case.task {
        Task::PtMolecular | Task::PtElementInventory => {
            builder.at_temperature(case.universe.temperature())
        }
        Task::PhMolecular | Task::PhElementInventory => {
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
        }
        Task::PtElementInventoryRange => builder
            .over_temperature_range(RANGE_TEMPERATURES.to_vec())
            .expect("range must be valid"),
    };
    rst_builder(builder, case.universe.reference_pressure(), case.rst_solver)
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
        panic!("expected P,T point")
    };
    let thermo =
        ResolvedThermochemistry::from_resolved_system(point.resolved()).expect("thermochemistry");
    thermo
        .enthalpy_model()
        .evaluate_total(point.solution().component_moles(), universe.temperature())
        .expect("reference enthalpy")
}

fn target_for(case: Case) -> f64 {
    reference_enthalpy(case.universe, case.task.is_elemental())
}

fn solve_case(
    case: Case,
    target_enthalpy: f64,
) -> Result<EquilibriumCalculatorOutcome, EquilibriumCalculatorError> {
    builder_for(case, target_enthalpy).solve()
}

#[cfg(feature = "equilibrium-benchmark-rst-frontends")]
#[derive(Debug, Clone)]
struct FrontendSignature {
    temperature_k: f64,
    mole_fractions: Vec<f64>,
    phase_count: usize,
    accepted_backend: String,
    rst_preparation: Option<RstPreparationTimingReport>,
}

#[cfg(feature = "equilibrium-benchmark-rst-frontends")]
fn frontend_signature(outcome: &EquilibriumCalculatorOutcome) -> FrontendSignature {
    let (temperature_k, solution, rst_preparation) = match outcome {
        EquilibriumCalculatorOutcome::PtPoint(point) => (
            point.solution().conditions().temperature(),
            point.solution(),
            point.solution().timing_report().rst_preparation_timing(),
        ),
        EquilibriumCalculatorOutcome::PhPoint(point) => (
            point.solution().temperature(),
            point.solution().equilibrium(),
            point.solution().report().timing().rst_preparation_timing(),
        ),
        other => panic!("frontend A/B only accepts point outcomes, got {other:?}"),
    };
    let moles = solution.component_moles();
    let total = moles.iter().sum::<f64>();
    FrontendSignature {
        temperature_k,
        mole_fractions: moles.iter().map(|amount| amount / total).collect(),
        phase_count: solution.phases().len(),
        accepted_backend: format!("{:?}", solution.solve_report().accepted_backend),
        rst_preparation: (!rst_preparation.total().is_zero()).then_some(rst_preparation),
    }
}

#[cfg(feature = "equilibrium-benchmark-rst-frontends")]
#[derive(Tabled)]
struct FrontendPreparationStageRow {
    #[tabled(rename = "case")]
    case: String,
    #[tabled(rename = "stage")]
    stage: String,
    #[tabled(rename = "ExprLegacy ms (% prep)")]
    expr_legacy: String,
    #[tabled(rename = "AtomView ms (% prep)")]
    atom_view: String,
    #[tabled(rename = "AtomView vs ExprLegacy")]
    delta: String,
    #[tabled(rename = "scope")]
    scope: String,
}

#[cfg(feature = "equilibrium-benchmark-rst-frontends")]
#[derive(Tabled)]
struct NativeSolverControlRow {
    #[tabled(rename = "case")]
    case: String,
    #[tabled(rename = "status")]
    status: String,
    #[tabled(rename = "accepted backend")]
    backend: String,
    #[tabled(rename = "facade ms")]
    facade_ms: String,
    #[tabled(rename = "RST prep")]
    rst_preparation: String,
    #[tabled(rename = "T K")]
    temperature_k: String,
    #[tabled(rename = "delta T vs ExprLegacy K")]
    temperature_delta_k: String,
    #[tabled(rename = "phase count")]
    phase_count: usize,
    #[tabled(rename = "max fraction delta vs ExprLegacy")]
    fraction_delta: String,
    #[tabled(rename = "detail")]
    detail: String,
}

#[cfg(feature = "equilibrium-benchmark-rst-frontends")]
fn frontend_stage_rows(
    case: &str,
    expr: RstPreparationTimingReport,
    atom: RstPreparationTimingReport,
) -> Vec<FrontendPreparationStageRow> {
    let detailed =
        expr.backend_stage_breakdown_available() && atom.backend_stage_breakdown_available();
    let mut stages = vec![
        (
            "KiThe symbolic input validation",
            Some(expr.symbolic_input_validation()),
            Some(atom.symbolic_input_validation()),
            "KiThe preparation child",
        ),
        (
            "KiThe phase topology",
            Some(expr.symbolic_phase_topology()),
            Some(atom.symbolic_phase_topology()),
            "KiThe preparation child",
        ),
        (
            "KiThe variable expansion",
            Some(expr.symbolic_variable_expansion()),
            Some(atom.symbolic_variable_expansion()),
            "KiThe preparation child",
        ),
        (
            "KiThe phase totals",
            Some(expr.symbolic_phase_totals()),
            Some(atom.symbolic_phase_totals()),
            "KiThe preparation child",
        ),
        (
            "KiThe activity offsets",
            Some(expr.symbolic_activity_offsets()),
            Some(atom.symbolic_activity_offsets()),
            "KiThe preparation child",
        ),
        (
            "KiThe reaction assembly",
            Some(expr.symbolic_reaction_assembly()),
            Some(atom.symbolic_reaction_assembly()),
            "KiThe preparation child",
        ),
        (
            "KiThe element-balance assembly",
            Some(expr.symbolic_element_balance_assembly()),
            Some(atom.symbolic_element_balance_assembly()),
            "KiThe preparation child",
        ),
        (
            "KiThe residual simplification",
            Some(expr.symbolic_residual_simplification()),
            Some(atom.symbolic_residual_simplification()),
            "KiThe preparation child",
        ),
        (
            "KiThe equation generation",
            Some(expr.equation_generation()),
            Some(atom.equation_generation()),
            "KiThe preparation child",
        ),
        (
            "KiThe temperature substitution",
            Some(expr.temperature_substitution()),
            Some(atom.temperature_substitution()),
            "KiThe preparation child",
        ),
        (
            "KiThe row scaling",
            Some(expr.row_scaling()),
            Some(atom.row_scaling()),
            "KiThe preparation child",
        ),
        (
            "KiThe enthalpy equation",
            Some(expr.enthalpy_equation()),
            Some(atom.enthalpy_equation()),
            "KiThe preparation child; zero for P,T",
        ),
        (
            "KiThe variable materialization",
            Some(expr.variable_materialization()),
            Some(atom.variable_materialization()),
            "KiThe preparation child",
        ),
        (
            "KiThe parameter evaluation",
            Some(expr.parameter_evaluation()),
            Some(atom.parameter_evaluation()),
            "KiThe preparation child",
        ),
        (
            "KiThe log-mole bounds",
            Some(expr.log_mole_bounds()),
            Some(atom.log_mole_bounds()),
            "KiThe preparation child",
        ),
        (
            "RST problem materialization aggregate",
            Some(expr.problem_materialization()),
            Some(atom.problem_materialization()),
            "aggregate; do not add its RST child rows",
        ),
        (
            "RST symbolic differentiation",
            detailed.then_some(expr.symbolic_differentiation()),
            detailed.then_some(atom.symbolic_differentiation()),
            "RST-reported child of aggregate",
        ),
        (
            "RST residual lambdification",
            detailed.then_some(expr.residual_lambdification()),
            detailed.then_some(atom.residual_lambdification()),
            "RST-reported child of aggregate",
        ),
        (
            "RST Jacobian lambdification",
            detailed.then_some(expr.jacobian_lambdification()),
            detailed.then_some(atom.jacobian_lambdification()),
            "RST-reported child of aggregate",
        ),
    ];
    for (name, expr_duration, atom_duration) in [
        (
            "RST input validation",
            expr.rst_input_validation(),
            atom.rst_input_validation(),
        ),
        (
            "RST expression parsing",
            expr.rst_expression_parsing(),
            atom.rst_expression_parsing(),
        ),
        (
            "RST expression graph materialization",
            expr.rst_expression_graph_materialization(),
            atom.rst_expression_graph_materialization(),
        ),
        (
            "RST parameter binding",
            expr.rst_parameter_binding(),
            atom.rst_parameter_binding(),
        ),
        (
            "RST prepared problem assembly",
            expr.rst_problem_assembly(),
            atom.rst_problem_assembly(),
        ),
        (
            "RST AOT materialization",
            expr.rst_aot_materialization(),
            atom.rst_aot_materialization(),
        ),
        (
            "RST final validation",
            expr.rst_final_validation(),
            atom.rst_final_validation(),
        ),
        (
            "RST unattributed preparation",
            expr.rst_unattributed(),
            atom.rst_unattributed(),
        ),
    ] {
        stages.push((
            name,
            expr_duration,
            atom_duration,
            "RST-reported child of aggregate",
        ));
    }

    let format_value = |value: Option<Duration>, total: Duration| {
        value.map_or_else(
            || "n/a".to_owned(),
            |duration| {
                format!(
                    "{:.3} ({})",
                    duration.as_secs_f64() * 1.0e3,
                    percent(micros(duration), micros(total))
                )
            },
        )
    };
    stages
        .into_iter()
        .map(|(stage, expr_duration, atom_duration, scope)| {
            let delta = match (expr_duration, atom_duration) {
                (Some(expr), Some(atom)) if !expr.is_zero() => format!(
                    "{:+.1}%",
                    (atom.as_secs_f64() / expr.as_secs_f64() - 1.0) * 100.0
                ),
                (Some(expr), Some(atom)) if expr.is_zero() && atom.is_zero() => "0.0%".to_owned(),
                (Some(_), Some(_)) => "n/a (Expr=0)".to_owned(),
                _ => "n/a".to_owned(),
            };
            FrontendPreparationStageRow {
                case: case.to_owned(),
                stage: stage.to_owned(),
                expr_legacy: format_value(expr_duration, expr.total()),
                atom_view: format_value(atom_duration, atom.total()),
                delta,
                scope: scope.to_owned(),
            }
        })
        .collect()
}

#[cfg(feature = "equilibrium-benchmark-rst-frontends")]
fn frontend_cases() -> [Case; 6] {
    [
        Case {
            universe: Universe::Small3,
            task: Task::PtMolecular,
            rst_solver: RstSolver::MinpackLm,
        },
        Case {
            universe: Universe::Stanjan15,
            task: Task::PtElementInventory,
            rst_solver: RstSolver::MinpackLm,
        },
        Case {
            universe: Universe::Tp1907_18,
            task: Task::PtElementInventory,
            rst_solver: RstSolver::MinpackLm,
        },
        Case {
            universe: Universe::Rp1311_40,
            task: Task::PtElementInventory,
            rst_solver: RstSolver::MinpackLm,
        },
        Case {
            universe: Universe::Stanjan15,
            task: Task::PhMolecular,
            rst_solver: RstSolver::MinpackLm,
        },
        Case {
            universe: Universe::Rp1311_40,
            task: Task::PhElementInventory,
            rst_solver: RstSolver::MinpackLm,
        },
    ]
}

#[cfg(feature = "equilibrium-benchmark-rst-frontends")]
fn run_frontend(
    case: Case,
    frontend: SymbolicLambdifyFrontend,
) -> Result<
    (Duration, FrontendSignature),
    KiThe::Thermodynamics::ChemEquilibrium::equilibrium_calculator::EquilibriumCalculatorError,
> {
    let target = target_for(case);
    let started = Instant::now();
    let frontend = match frontend {
        SymbolicLambdifyFrontend::ExprLegacy => {
            KiThe::Thermodynamics::ChemEquilibrium::equilibrium_rst_backend::RstSymbolicFrontend::ExprLegacy
        }
        SymbolicLambdifyFrontend::AtomViewNative => {
            KiThe::Thermodynamics::ChemEquilibrium::equilibrium_rst_backend::RstSymbolicFrontend::AtomView
        }
    };
    let outcome = builder_for(case, target)
        .rst_symbolic_frontend(frontend)
        .solve()?;
    let elapsed = started.elapsed();
    Ok((elapsed, frontend_signature(&outcome)))
}

#[cfg(feature = "equilibrium-benchmark-rst-frontends")]
fn run_native_legacy_lm(
    case: Case,
) -> Result<
    (Duration, FrontendSignature),
    KiThe::Thermodynamics::ChemEquilibrium::equilibrium_calculator::EquilibriumCalculatorError,
> {
    let target = target_for(case);
    let started = Instant::now();
    let builder = builder_for(case, target)
        .solver_policy(SolverPolicy::Single(SolverBackend::Legacy(
            LegacyEquilibriumSolver::LM,
        )))
        .expect("single KiThe-native Legacy LM policy must validate");
    let outcome = builder.solve()?;
    let elapsed = started.elapsed();
    Ok((elapsed, frontend_signature(&outcome)))
}

#[cfg(feature = "equilibrium-benchmark-rst-frontends")]
fn summarize_frontend_error(
    error: &KiThe::Thermodynamics::ChemEquilibrium::equilibrium_calculator::EquilibriumCalculatorError,
) -> String {
    use KiThe::Thermodynamics::ChemEquilibrium::equilibrium_calculator::EquilibriumCalculatorError;

    let mut summary = error.to_string();
    let source = match error {
        EquilibriumCalculatorError::Solve(source)
        | EquilibriumCalculatorError::Preparation(source) => Some(source),
        _ => None,
    };
    if let Some(source) = source {
        for attempt in source.backend_attempt_reports() {
            let status = match &attempt.outcome {
                KiThe::Thermodynamics::ChemEquilibrium::equilibrium_solver_policy::SolverAttemptOutcome::Accepted => "accepted".to_owned(),
                KiThe::Thermodynamics::ChemEquilibrium::equilibrium_solver_policy::SolverAttemptOutcome::Failed { kind, reason, .. } => format!("failed:{kind:?}:{reason}"),
                KiThe::Thermodynamics::ChemEquilibrium::equilibrium_solver_policy::SolverAttemptOutcome::RejectedCandidate { reason } => format!("rejected:{reason}"),
                KiThe::Thermodynamics::ChemEquilibrium::equilibrium_solver_policy::SolverAttemptOutcome::Skipped { reason } => format!("skipped:{reason}"),
            };
            let metrics = attempt.metrics.as_ref().map_or_else(
                || "metrics=n/a".to_owned(),
                |metrics| format!(
                    "termination={:?},iterations={},residual_calls={},jacobian_calls={},elapsed_ms={}",
                    metrics.termination,
                    metrics.iterations,
                    metrics.residual_evaluations,
                    metrics.jacobian_evaluations,
                    metrics.elapsed_millis,
                ),
            );
            summary.push_str(&format!(
                "\n  backend={:?},outcome={status},{metrics}",
                attempt.backend,
            ));
            if let Some(iterations) = attempt
                .metrics
                .as_ref()
                .and_then(|metrics| metrics.evaluation_timing.as_ref())
                .and_then(|timing| timing.iteration_telemetry.as_ref())
            {
                summary.push_str("\n  iteration | residual calls | jacobian calls | linear solves | accepted | rejected | residual us | jacobian us | linear us");
                for iteration in iterations {
                    summary.push_str(&format!(
                        "\n  {} | {} | {} | {} | {} | {} | {} | {} | {}",
                        iteration.iteration,
                        iteration.residual_evaluations,
                        iteration.jacobian_evaluations,
                        iteration.linear_solves,
                        iteration.accepted_steps,
                        iteration.rejected_steps,
                        iteration.residual_micros,
                        iteration.jacobian_micros,
                        iteration.linear_solve_micros,
                    ));
                }
            }
        }
    }
    summary
}

#[cfg(feature = "equilibrium-benchmark-rst-frontends")]
fn report_frontend_comparison(c: &mut Criterion) {
    let frontends = [
        ("expr_legacy", SymbolicLambdifyFrontend::ExprLegacy),
        ("atom_view", SymbolicLambdifyFrontend::AtomViewNative),
    ];
    let mut group = c.benchmark_group("equilibrium_rst_frontend_ab");
    group.sample_size(10);
    group.warm_up_time(Duration::from_secs(1));
    group.measurement_time(Duration::from_secs(2));

    println!(
        "\nPERF-RST-FRONTEND matched facade routes; fixed RST solver=Minpack LM, plus KiThe Legacy LM controls"
    );
    println!(
        "Three passes: accepted full-solve A/B; matched preparation-stage attribution; same-input KiThe Legacy LM control. Criterion timings are emitted only for accepted, output-matched frontend pairs."
    );
    println!(
        "case | frontend | accepted_backend | status | facade_ms | rst_preparation_ms | temperature_K | delta_T_K | phase_count | max_fraction_delta_vs_expr"
    );
    let mut preparation_stage_rows = Vec::new();
    let mut preparation_stage_unavailable = Vec::new();
    let mut native_lm_rows = Vec::new();
    let selected_case = std::env::var("RST_FRONTEND_CASE").ok();
    let selected_frontend = std::env::var("RST_FRONTEND_ONLY").ok();
    let report_only = std::env::var_os("RST_FRONTEND_REPORT_ONLY").is_some();
    for case in frontend_cases().into_iter().filter(|case| {
        selected_case.as_deref().is_none_or(|selected| {
            format!("{}_{}", case.universe.label(), case.task.label()) == selected
        })
    }) {
        let selected_frontends = frontends.iter().filter(|(label, _)| {
            selected_frontend
                .as_deref()
                .is_none_or(|selected| *label == selected)
        });
        let outcomes = selected_frontends
            .map(|(label, frontend)| (*label, *frontend, run_frontend(case, *frontend)))
            .collect::<Vec<_>>();
        let expr_preparation = outcomes
            .iter()
            .find(|(label, _, result)| *label == "expr_legacy" && result.is_ok())
            .and_then(|(_, _, result)| result.as_ref().ok().map(|(_, signature)| signature));
        let atom_preparation = outcomes
            .iter()
            .find(|(label, _, result)| *label == "atom_view" && result.is_ok())
            .and_then(|(_, _, result)| result.as_ref().ok().map(|(_, signature)| signature));
        if let (Some(expr), Some(atom)) = (expr_preparation, atom_preparation) {
            match (expr.rst_preparation, atom.rst_preparation) {
                (Some(expr_timing), Some(atom_timing))
                    if expr_timing.backend_stage_breakdown_available()
                        && atom_timing.backend_stage_breakdown_available() =>
                {
                    preparation_stage_rows.extend(frontend_stage_rows(
                        &format!("{}_{}", case.universe.label(), case.task.label()),
                        expr_timing,
                        atom_timing,
                    ));
                }
                _ => preparation_stage_unavailable.push(format!(
                    "{}_{}: detailed RST preparation timing unavailable for at least one frontend",
                    case.universe.label(),
                    case.task.label()
                )),
            }
        }
        if case.universe == Universe::Stanjan15 && case.task == Task::PhMolecular {
            for (label, _, outcome) in &outcomes {
                assert!(
                    outcome.is_ok(),
                    "STANJAN-15 monolithic P,H must be accepted with {label}: {}",
                    outcome
                        .as_ref()
                        .err()
                        .map(summarize_frontend_error)
                        .unwrap_or_default(),
                );
            }
        }
        let reference = outcomes
            .iter()
            .find(|(label, _, result)| *label == "expr_legacy" && result.is_ok())
            .and_then(|(_, _, result)| {
                result.as_ref().ok().map(|(_, signature)| signature.clone())
            });
        let atom = outcomes
            .iter()
            .find(|(label, _, _)| *label == "atom_view")
            .and_then(|(_, _, result)| {
                result.as_ref().ok().map(|(_, signature)| signature.clone())
            });
        let parity = match (&reference, &atom) {
            (Some(expr), Some(atom)) => {
                assert_eq!(
                    expr.accepted_backend,
                    atom.accepted_backend,
                    "{} symbolic frontends were accepted by different nonlinear backends",
                    case.label()
                );
                assert_eq!(expr.mole_fractions.len(), atom.mole_fractions.len());
                let fraction_delta = expr
                    .mole_fractions
                    .iter()
                    .zip(&atom.mole_fractions)
                    .map(|(left, right)| (left - right).abs())
                    .fold(0.0, f64::max);
                let temperature_delta = (expr.temperature_k - atom.temperature_k).abs();
                let topology_matches = expr.phase_count == atom.phase_count;
                assert!(
                    temperature_delta <= 1.0e-3,
                    "{} temperature mismatch: {temperature_delta:e} K",
                    case.label()
                );
                assert!(
                    topology_matches,
                    "{} phase topology differs between symbolic frontends",
                    case.label()
                );
                assert!(
                    fraction_delta <= 1.0e-8,
                    "{} published mole-fraction mismatch: {fraction_delta:e}",
                    case.label()
                );
                Some((temperature_delta, fraction_delta))
            }
            _ => None,
        };

        for (label, frontend, outcome) in &outcomes {
            let case_name = format!("{}_{}", case.universe.label(), case.task.label());
            match outcome {
                Ok((elapsed, signature)) => {
                    let prep = signature.rst_preparation.map_or_else(
                        || "n/a".to_owned(),
                        |timing| format!("{:.3}", timing.total().as_secs_f64() * 1.0e3),
                    );
                    let (temperature_delta, fraction_delta) = parity.map_or_else(
                        || ("n/a".to_owned(), "n/a".to_owned()),
                        |(temperature, fraction)| {
                            (format!("{temperature:.3e}"), format!("{fraction:.6e}"))
                        },
                    );
                    println!(
                        "{case_name} | {label} | {} | accepted | {:.3} | {prep} | {:.9} | {temperature_delta} | {} | {fraction_delta}",
                        signature.accepted_backend,
                        elapsed.as_secs_f64() * 1.0e3,
                        signature.temperature_k,
                        signature.phase_count,
                    );
                    if !report_only && parity.is_some() {
                        group.bench_function(BenchmarkId::new(*label, &case_name), |bencher| {
                            bencher.iter(|| {
                                black_box(run_frontend(case, *frontend).unwrap_or_else(|error| {
                                    panic!("{case_name} {label} failed during sampling: {error}")
                                }));
                            });
                        });
                    }
                }
                Err(error) => {
                    println!(
                        "{case_name} | {label} | n/a | FAILED | n/a | n/a | n/a | n/a | n/a | n/a\n  {}",
                        summarize_frontend_error(&error)
                    );
                }
            }
        }

        match run_native_legacy_lm(case) {
            Ok((elapsed, signature)) => {
                let (status, fraction_delta) = if let Some(expr) = &reference {
                    assert_eq!(
                        expr.mole_fractions.len(),
                        signature.mole_fractions.len(),
                        "{} native LM component dimension differs from RST ExprLegacy",
                        case.label()
                    );
                    let fraction_delta = expr
                        .mole_fractions
                        .iter()
                        .zip(&signature.mole_fractions)
                        .map(|(left, right)| (left - right).abs())
                        .fold(0.0_f64, f64::max);
                    let temperature_delta = (expr.temperature_k - signature.temperature_k).abs();
                    let parity = fraction_delta <= 1.0e-8
                        && temperature_delta <= 1.0e-3
                        && expr.phase_count == signature.phase_count;
                    (
                        if parity {
                            "accepted/parity"
                        } else {
                            "accepted/mismatch"
                        },
                        format!("{fraction_delta:.6e}"),
                    )
                } else {
                    ("accepted/no-RST-reference", "n/a".to_owned())
                };
                native_lm_rows.push(NativeSolverControlRow {
                    case: format!("{}_{}", case.universe.label(), case.task.label()),
                    status: status.to_owned(),
                    backend: signature.accepted_backend.clone(),
                    facade_ms: format!("{:.3}", elapsed.as_secs_f64() * 1.0e3),
                    rst_preparation: "n/a (KiThe native)".to_owned(),
                    temperature_k: format!("{:.9}", signature.temperature_k),
                    temperature_delta_k: reference.as_ref().map_or_else(
                        || "n/a".to_owned(),
                        |expr| {
                            format!(
                                "{:.6e}",
                                (expr.temperature_k - signature.temperature_k).abs()
                            )
                        },
                    ),
                    phase_count: signature.phase_count,
                    fraction_delta,
                    detail: "same facade input; compared against ExprLegacy RST result".to_owned(),
                });
                if status == "accepted/parity" && !report_only {
                    let case_name = format!("{}_{}", case.universe.label(), case.task.label());
                    group.bench_function(
                            BenchmarkId::new("ki_the_legacy_lm", &case_name),
                            |bencher| {
                                bencher.iter(|| {
                                    black_box(run_native_legacy_lm(case).unwrap_or_else(|error| {
                                        panic!("{case_name} KiThe Legacy LM failed during sampling: {error}")
                                    }));
                                });
                            },
                        );
                }
            }
            Err(error) => native_lm_rows.push(NativeSolverControlRow {
                case: format!("{}_{}", case.universe.label(), case.task.label()),
                status: "failed".to_owned(),
                backend: "n/a".to_owned(),
                facade_ms: "n/a".to_owned(),
                rst_preparation: "n/a (KiThe native)".to_owned(),
                temperature_k: "n/a".to_owned(),
                temperature_delta_k: "n/a".to_owned(),
                phase_count: 0,
                fraction_delta: "n/a".to_owned(),
                detail: summarize_frontend_error(&error),
            }),
        }
    }
    if !preparation_stage_rows.is_empty() {
        println!(
            "\nPERF-RST-FRONTEND matched preparation stages (one diagnostic pass per frontend; child rows are nested and not additive)\n{}",
            Table::new(preparation_stage_rows)
        );
    }
    if !preparation_stage_unavailable.is_empty() {
        println!(
            "\nPERF-RST-FRONTEND preparation-stage detail unavailable\n{}",
            preparation_stage_unavailable.join("\n")
        );
    }
    if !native_lm_rows.is_empty() {
        println!(
            "\nPERF-RST-FRONTEND KiThe-native Legacy LM controls (same public facade inputs; no RST preparation)\n{}",
            Table::new(native_lm_rows)
        );
    }
    if !report_only {
        group.finish();
    }
}

fn collect_case(case: Case, target_enthalpy: f64) -> Evidence {
    let started = Instant::now();
    let case_name = case.label();
    let mut evidence = Evidence {
        case: case_name.clone(),
        route: match case.task {
            Task::PhMolecular | Task::PhElementInventory => "P,H".to_owned(),
            _ => "P,T".to_owned(),
        },
        policy: case.rst_solver.label().to_owned(),
        status: "accepted".to_owned(),
        wall_us: 0.0,
        accepted_backend: "not accepted".to_owned(),
        ..Evidence::default()
    };

    let outcome = match solve_case(case, target_enthalpy) {
        Ok(outcome) => outcome,
        Err(error) => {
            evidence.wall_us = micros(started.elapsed());
            evidence.status = match &error {
                EquilibriumCalculatorError::Preparation(cause)
                | EquilibriumCalculatorError::Solve(cause) => {
                    format!("failed:{:?}", cause.kind())
                }
                _ => "failed:facade".to_owned(),
            };
            let attempts = match &error {
                EquilibriumCalculatorError::Preparation(cause)
                | EquilibriumCalculatorError::Solve(cause) => cause.backend_attempt_reports(),
                _ => Vec::new(),
            };
            attempt_rows_slice(
                &mut evidence.attempts,
                &case_name,
                "failed RST solve",
                &attempts,
            );
            if let EquilibriumCalculatorError::Preparation(cause)
            | EquilibriumCalculatorError::Solve(cause) = &error
            {
                evidence.iterations = cause.nonlinear_iterations();
            }
            return evidence;
        }
    };

    match outcome {
        EquilibriumCalculatorOutcome::PtPoint(point) => {
            evidence.wall_us = micros(point.instrumentation_report().elapsed());
            let facade = *point.facade_timing();
            let inner = *point.solution().timing_report();
            evidence.workflow_us = micros(inner.total());
            facade_stages(&mut evidence.stages, &case_name, facade, evidence.wall_us);
            equilibrium_stages(
                &mut evidence.stages,
                &case_name,
                "P,T solve",
                inner,
                evidence.workflow_us,
            );
            rst_preparation_stages(
                &mut evidence.stages,
                &case_name,
                "P,T RST preparation",
                inner.rst_preparation_timing(),
            );
            add_rst_preparation_sample(
                &mut evidence,
                &case_name,
                "P,T",
                inner.rst_preparation_timing(),
            );
            add_solution_evidence(&mut evidence, point.solution(), "P,T backend");
            evidence.residual_evaluations = point
                .instrumentation_report()
                .counters()
                .residual_evaluations();
            evidence.jacobian_evaluations = point
                .instrumentation_report()
                .counters()
                .jacobian_evaluations();
        }
        EquilibriumCalculatorOutcome::PhPoint(point) => {
            evidence.wall_us = micros(point.instrumentation_report().elapsed());
            let facade = *point.facade_timing();
            let report = point.solution().report();
            let outer = report.timing();
            let inner = report.inner_timing();
            evidence.workflow_us = micros(outer.total());
            facade_stages(&mut evidence.stages, &case_name, facade, evidence.wall_us);
            ph_stages(
                &mut evidence.stages,
                &case_name,
                outer,
                evidence.workflow_us,
            );
            equilibrium_stages(
                &mut evidence.stages,
                &case_name,
                "P,H inner solve",
                inner,
                micros(inner.total()),
            );
            rst_preparation_stages(
                &mut evidence.stages,
                &case_name,
                "P,H RST preparation",
                outer.rst_preparation_timing(),
            );
            add_rst_preparation_sample(
                &mut evidence,
                &case_name,
                "P,H outer",
                outer.rst_preparation_timing(),
            );
            if report.monolithic_evidence().is_some() {
                add_solution_evidence(
                    &mut evidence,
                    point.solution().equilibrium(),
                    "P,H monolithic backend",
                );
            } else {
                for (index, trial) in report.trials().iter().enumerate() {
                    let trial_timing = trial.timing();
                    let trial_scope = micros(trial_timing.total());
                    let template_timing = trial_timing.template_preparation();
                    rst_preparation_stages(
                        &mut evidence.stages,
                        &case_name,
                        &format!("trial-{index} template RST preparation"),
                        template_timing.rst_preparation_timing(),
                    );
                    add_rst_preparation_sample(
                        &mut evidence,
                        &case_name,
                        &format!("trial-{index} template"),
                        template_timing.rst_preparation_timing(),
                    );
                    stage(
                        &mut evidence.stages,
                        &case_name,
                        &format!("trial-{index}"),
                        "trial total",
                        trial_timing.total(),
                        trial_scope,
                        "one scalar P,H trial",
                    );
                    stage(
                        &mut evidence.stages,
                        &case_name,
                        &format!("trial-{index}"),
                        "inner execution",
                        trial_timing.inner_execution(),
                        trial_scope,
                        "trial-local wall stage",
                    );
                    stage(
                        &mut evidence.stages,
                        &case_name,
                        &format!("trial-{index}"),
                        "enthalpy evaluation",
                        trial_timing.enthalpy_evaluation(),
                        trial_scope,
                        "trial-local thermochemistry",
                    );
                    stage(
                        &mut evidence.stages,
                        &case_name,
                        &format!("trial-{index}"),
                        "postprocessing",
                        trial_timing.postprocessing(),
                        trial_scope,
                        "trial-local evidence publication",
                    );
                    equilibrium_stages(
                        &mut evidence.stages,
                        &case_name,
                        &format!("trial-{index} inner"),
                        trial.inner_timing(),
                        micros(trial.inner_timing().total()),
                    );
                    rst_preparation_stages(
                        &mut evidence.stages,
                        &case_name,
                        &format!("trial-{index} RST preparation"),
                        trial.inner_timing().rst_preparation_timing(),
                    );
                    add_rst_preparation_sample(
                        &mut evidence,
                        &case_name,
                        &format!("trial-{index}"),
                        trial.inner_timing().rst_preparation_timing(),
                    );
                    if let Some(inner_evidence) = trial.inner_evidence() {
                        attempt_rows(
                            &mut evidence.attempts,
                            &case_name,
                            &format!("trial-{index}"),
                            inner_evidence.solve_report(),
                        );
                    }
                }
            }
            if report.monolithic_evidence().is_none() {
                summarize_solution(&mut evidence, point.solution().equilibrium());
            }
            evidence.residual_evaluations = point
                .instrumentation_report()
                .counters()
                .residual_evaluations();
            evidence.jacobian_evaluations = point
                .instrumentation_report()
                .counters()
                .jacobian_evaluations();
        }
        EquilibriumCalculatorOutcome::PtRange(range) => {
            evidence.wall_us = micros(range.instrumentation_report().elapsed());
            let report = range.solution().report();
            evidence.workflow_us = micros(report.total());
            facade_stages(
                &mut evidence.stages,
                &case_name,
                *range.facade_timing(),
                evidence.wall_us,
            );
            let initial = report.initial_formulation_timing();
            equilibrium_stages(
                &mut evidence.stages,
                &case_name,
                "range initial formulation",
                *initial,
                micros(initial.total()),
            );
            rst_preparation_stages(
                &mut evidence.stages,
                &case_name,
                "range initial RST preparation",
                initial.rst_preparation_timing(),
            );
            add_rst_preparation_sample(
                &mut evidence,
                &case_name,
                "range initial",
                initial.rst_preparation_timing(),
            );
            for (index, point) in range.solution().points().iter().enumerate() {
                let timing = point.report().timing();
                let scope = micros(timing.total());
                equilibrium_stages(
                    &mut evidence.stages,
                    &case_name,
                    &format!("range point-{index}"),
                    *timing,
                    scope,
                );
                rst_preparation_stages(
                    &mut evidence.stages,
                    &case_name,
                    &format!("range point-{index} RST preparation"),
                    timing.rst_preparation_timing(),
                );
                add_rst_preparation_sample(
                    &mut evidence,
                    &case_name,
                    &format!("range point-{index}"),
                    timing.rst_preparation_timing(),
                );
                add_solution_evidence(
                    &mut evidence,
                    point.solution(),
                    &format!("range point-{index}"),
                );
            }
            evidence.formulation_builds += report.formulation_builds();
            evidence.formulation_reuses = report.formulation_reuses();
            evidence.symbolic_parameter_updates = report.symbolic_parameter_updates();
            evidence.residual_evaluations = range
                .instrumentation_report()
                .counters()
                .residual_evaluations();
            evidence.jacobian_evaluations = range
                .instrumentation_report()
                .counters()
                .jacobian_evaluations();
        }
        EquilibriumCalculatorOutcome::PhRange(_) => {
            panic!("RST preparation bench uses no P,H range case")
        }
    }
    evidence
}

fn percent(value_us: f64, scope_us: f64) -> String {
    if scope_us <= f64::EPSILON {
        "0.0%".to_owned()
    } else {
        format!("{:.1}%", 100.0 * value_us.max(0.0) / scope_us)
    }
}

#[derive(Tabled)]
struct MacroRow {
    #[tabled(rename = "case")]
    case: String,
    #[tabled(rename = "route")]
    route: String,
    #[tabled(rename = "wall")]
    wall: String,
    #[tabled(rename = "outside workflow")]
    outside: String,
    #[tabled(rename = "workflow")]
    workflow: String,
    #[tabled(rename = "backend")]
    backend: String,
    #[tabled(rename = "policy")]
    policy: String,
    #[tabled(rename = "status")]
    status: String,
    #[tabled(rename = "attempts")]
    attempts: usize,
}

#[derive(Tabled)]
struct StageRow {
    #[tabled(rename = "case")]
    case: String,
    #[tabled(rename = "scope")]
    scope: String,
    #[tabled(rename = "stage")]
    stage: String,
    #[tabled(rename = "time")]
    time: String,
    #[tabled(rename = "% of scope")]
    percent: String,
    #[tabled(rename = "semantics")]
    semantics: String,
}

#[derive(Tabled)]
struct AttemptRow {
    #[tabled(rename = "case")]
    case: String,
    #[tabled(rename = "scope")]
    scope: String,
    #[tabled(rename = "#")]
    index: usize,
    #[tabled(rename = "backend")]
    backend: String,
    #[tabled(rename = "outcome")]
    outcome: String,
    #[tabled(rename = "wall")]
    wall: String,
    #[tabled(rename = "iter")]
    iterations: usize,
    #[tabled(rename = "res calls")]
    residual_calls: usize,
    #[tabled(rename = "res")]
    residual: String,
    #[tabled(rename = "Jac calls")]
    jacobian_calls: usize,
    #[tabled(rename = "Jac")]
    jacobian: String,
    #[tabled(rename = "linear calls")]
    linear_solves: usize,
    #[tabled(rename = "linear")]
    linear: String,
    #[tabled(rename = "overhead")]
    overhead: String,
}

#[derive(Tabled)]
struct ExecutionRow {
    #[tabled(rename = "case")]
    case: String,
    #[tabled(rename = "route")]
    route: String,
    #[tabled(rename = "RST backend")]
    backend: String,
    #[tabled(rename = "iterations")]
    iterations: usize,
    #[tabled(rename = "residual evals")]
    residual: usize,
    #[tabled(rename = "Jacobian evals")]
    jacobian: usize,
    #[tabled(rename = "formulation builds")]
    builds: usize,
    #[tabled(rename = "formulation reuses")]
    reuses: usize,
    #[tabled(rename = "symbolic T updates")]
    symbolic_updates: usize,
    #[tabled(rename = "phase transitions")]
    transitions: usize,
}

#[derive(Tabled)]
struct RstPreparationRow {
    #[tabled(rename = "case")]
    case: String,
    #[tabled(rename = "scope")]
    scope: String,
    #[tabled(rename = "effective backend")]
    backend: String,
    #[tabled(rename = "artifact action")]
    artifact_action: String,
    #[tabled(rename = "RST prep aggregate")]
    preparation: String,
    #[tabled(rename = "AOT build")]
    build: String,
    #[tabled(rename = "residual jobs")]
    residual_jobs: String,
    #[tabled(rename = "Jacobian jobs")]
    jacobian_jobs: String,
    #[tabled(rename = "internal split")]
    internal_split: String,
}

#[derive(Tabled)]
struct CompactRow {
    #[tabled(rename = "case")]
    case: String,
    #[tabled(rename = "route")]
    route: String,
    #[tabled(rename = "method")]
    method: String,
    #[tabled(rename = "wall median")]
    wall: String,
    #[tabled(rename = "workflow median")]
    workflow: String,
    #[tabled(rename = "result")]
    result: String,
    #[tabled(rename = "attempts")]
    attempts: String,
}

fn median(values: &mut [f64]) -> f64 {
    values.sort_by(f64::total_cmp);
    values[values.len() / 2]
}

fn compact_method(policy: &str) -> &'static str {
    match policy {
        "rst_levenberg_marquardt" => "LM",
        "rst_backtracking_levenberg_marquardt" => "BT-LM",
        "rst_minpack_levenberg_marquardt" => "Minpack",
        "rst_nielsen_levenberg_marquardt" => "Nielsen",
        "rst_damped_newton" => "Damped NR",
        _ => "RST?",
    }
}

fn compact_report(evidence: &[Evidence]) {
    let mut rows = Vec::new();
    for chunk in evidence.chunks(REPORT_REPETITIONS) {
        if chunk.is_empty() {
            continue;
        }
        let first = &chunk[0];
        let mut wall: Vec<_> = chunk.iter().map(|item| item.wall_us).collect();
        let mut workflow: Vec<_> = chunk.iter().map(|item| item.workflow_us).collect();
        let status = if chunk.iter().all(|item| item.status == first.status) {
            if first.status == "accepted" {
                "PASS".to_owned()
            } else {
                "FAIL".to_owned()
            }
        } else {
            "MIXED".to_owned()
        };
        let attempts: Vec<_> = chunk
            .iter()
            .map(|item| item.attempts.len().to_string())
            .collect();
        let workflow = if first.status == "accepted" {
            format!("{:.2} ms", median(&mut workflow) / 1_000.0)
        } else {
            "n/a (rejected)".to_owned()
        };
        rows.push(CompactRow {
            case: first
                .case
                .split("_rst_")
                .next()
                .unwrap_or(&first.case)
                .to_owned(),
            route: first.route.clone(),
            method: compact_method(&first.policy).to_owned(),
            wall: format!("{:.2} ms", median(&mut wall) / 1_000.0),
            workflow,
            result: status,
            attempts: attempts.join("/") + " (3 reps)",
        });
    }
    println!(
        "\nPERF-RST-01 compact current matrix (median of three repetitions)\n{}",
        Table::new(rows)
    );
    println!(
        "FAIL rows are isolated-method rejection time, not successful-solve performance. Use the full report for stage and attempt attribution."
    );
}

fn print_report(evidence: &[Evidence]) {
    if std::env::var_os("RST_PREPARATION_COMPACT").is_some() {
        compact_report(evidence);
        return;
    }
    let macro_rows: Vec<_> = evidence
        .iter()
        .map(|item| MacroRow {
            case: item.case.clone(),
            route: item.route.clone(),
            wall: format!("{:.2} ms", item.wall_us / 1_000.0),
            outside: format!(
                "{:.2} ms / {}",
                (item.wall_us - item.workflow_us).max(0.0) / 1_000.0,
                percent(item.wall_us - item.workflow_us, item.wall_us)
            ),
            workflow: format!(
                "{:.2} ms / {}",
                item.workflow_us / 1_000.0,
                percent(item.workflow_us, item.wall_us)
            ),
            backend: item.accepted_backend.clone(),
            policy: item.policy.clone(),
            status: item.status.clone(),
            attempts: item.attempts.len(),
        })
        .collect();
    let stage_rows: Vec<_> = evidence
        .iter()
        .flat_map(|item| {
            item.stages.iter().map(|row| StageRow {
                case: row.case.clone(),
                scope: row.scope.clone(),
                stage: row.stage.clone(),
                time: format!("{:.2} ms", row.time_us / 1_000.0),
                percent: percent(row.time_us, row.scope_us),
                semantics: row.semantics.clone(),
            })
        })
        .collect();
    let attempt_rows: Vec<_> = evidence
        .iter()
        .flat_map(|item| {
            item.attempts.iter().map(|row| AttemptRow {
                case: row.case.clone(),
                scope: row.scope.clone(),
                index: row.index,
                backend: row.backend.clone(),
                outcome: row.outcome.clone(),
                wall: format!("{:.2} ms", row.wall_us / 1_000.0),
                iterations: row.iterations,
                residual_calls: row.residual_calls,
                residual: row.residual_us.map_or_else(
                    || "n/a".to_owned(),
                    |value| {
                        format!(
                            "{:.2} ms / {}",
                            value / 1_000.0,
                            percent(value, row.wall_us)
                        )
                    },
                ),
                jacobian_calls: row.jacobian_calls,
                jacobian: row.jacobian_us.map_or_else(
                    || "n/a".to_owned(),
                    |value| {
                        format!(
                            "{:.2} ms / {}",
                            value / 1_000.0,
                            percent(value, row.wall_us)
                        )
                    },
                ),
                linear_solves: row.linear_solves,
                linear: row.linear_us.map_or_else(
                    || "n/a".to_owned(),
                    |value| {
                        format!(
                            "{:.2} ms / {}",
                            value / 1_000.0,
                            percent(value, row.wall_us)
                        )
                    },
                ),
                overhead: row.overhead_us.map_or_else(
                    || "n/a".to_owned(),
                    |value| {
                        format!(
                            "{:.2} ms / {}",
                            value / 1_000.0,
                            percent(value, row.wall_us)
                        )
                    },
                ),
            })
        })
        .collect();
    let execution_rows: Vec<_> = evidence
        .iter()
        .map(|item| ExecutionRow {
            case: item.case.clone(),
            route: item.route.clone(),
            backend: item.accepted_backend.clone(),
            iterations: item.iterations,
            residual: item.residual_evaluations,
            jacobian: item.jacobian_evaluations,
            builds: item.formulation_builds,
            reuses: item.formulation_reuses,
            symbolic_updates: item.symbolic_parameter_updates,
            transitions: item.phase_transitions,
        })
        .collect();
    let rst_preparation_rows: Vec<_> = evidence
        .iter()
        .flat_map(|item| item.rst_preparations.iter())
        .map(|sample| {
            let report = sample.report;
            RstPreparationRow {
                case: sample.case.clone(),
                scope: sample.scope.clone(),
                backend: report.backend_kind().as_str().to_owned(),
                artifact_action: report.artifact_action().as_str().to_owned(),
                preparation: format!(
                    "{:.2} ms",
                    micros(report.problem_materialization()) / 1_000.0
                ),
                build: report.backend_build_duration().map_or_else(
                    || "n/a".to_owned(),
                    |duration| format!("{:.2} ms", micros(duration) / 1_000.0),
                ),
                residual_jobs: report
                    .generated_residual_jobs()
                    .map_or_else(|| "n/a".to_owned(), |value| value.to_string()),
                jacobian_jobs: report
                    .generated_jacobian_jobs()
                    .map_or_else(|| "n/a".to_owned(), |value| value.to_string()),
                internal_split: if report.backend_stage_breakdown_available() {
                    "available".to_owned()
                } else {
                    "unavailable in RST 0.4.19".to_owned()
                },
            }
        })
        .collect();

    let symbolic_simplification =
        match std::env::var("KITHE_RST_SYMBOLIC_SIMPLIFICATION").as_deref() {
            Ok("raw") | Ok("off") => "Raw (production)",
            Ok("simplified") | Ok("on") => "Simplified (baseline override)",
            _ if cfg!(feature = "equilibrium-benchmark-simplified-symbolic") => {
                "Simplified (benchmark baseline)"
            }
            _ => "Raw (production)",
        };
    println!(
        "\nPERF-RST-01 production_lambdify_policy=Sequential rst_symbolic_simplification={symbolic_simplification}\nThe adapter selects Sequential explicitly; Parallel is measured separately by PERF-RST-06."
    );
    println!(
        "\nPERF-RST-01 macro partition (workflow is macro-only)\n{}",
        Table::new(macro_rows)
    );
    println!(
        "\nPERF-RST-01 concrete stage ownership\n{}",
        Table::new(stage_rows)
    );
    println!(
        "\nPERF-RST-01 RST preparation evidence (metadata and detailed-stage availability)\n{}",
        Table::new(rst_preparation_rows)
    );
    println!(
        "\nPERF-RST-01 every RST backend attempt\n{}",
        Table::new(attempt_rows)
    );
    println!(
        "\nPERF-RST-01 execution/reuse evidence\n{}",
        Table::new(execution_rows)
    );
    println!(
        "\nStage percentages use the named scope denominator. Inclusive rows intentionally do not sum to 100%; the attempt ledger accounts for every RST call, including rejected attempts. RST 0.4.19 preparation children are printed below the aggregate row; use those exclusive child stages for bottleneck attribution."
    );
}

fn cases() -> Vec<Case> {
    let mut cases = Vec::new();
    let all_solvers = [
        RstSolver::Lm,
        RstSolver::BacktrackingLm,
        RstSolver::MinpackLm,
        RstSolver::NielsenLm,
        RstSolver::DampedNewton,
    ];
    for rst_solver in all_solvers {
        for (universe, task) in [
            (Universe::Small3, Task::PtMolecular),
            (Universe::Stanjan15, Task::PtMolecular),
            (Universe::Tp1907_18, Task::PtMolecular),
            (Universe::Rp1311_40, Task::PtMolecular),
            (Universe::Small3, Task::PtElementInventory),
            (Universe::Small3, Task::PhElementInventory),
            (Universe::Stanjan15, Task::PtElementInventory),
        ] {
            cases.push(Case {
                universe,
                task,
                rst_solver,
            });
        }
    }
    for rst_solver in [RstSolver::Lm, RstSolver::MinpackLm] {
        for (universe, task) in [
            (Universe::Small3, Task::PhMolecular),
            (Universe::Stanjan15, Task::PhMolecular),
            (Universe::Stanjan15, Task::PhElementInventory),
            (Universe::Tp1907_18, Task::PtElementInventory),
            (Universe::Tp1907_18, Task::PhElementInventory),
            (Universe::Rp1311_40, Task::PtElementInventory),
            (Universe::Rp1311_40, Task::PhElementInventory),
            (Universe::Small3, Task::PtElementInventoryRange),
        ] {
            cases.push(Case {
                universe,
                task,
                rst_solver,
            });
        }
    }
    cases
}

fn criterion_benchmark(c: &mut Criterion) {
    #[cfg(feature = "equilibrium-benchmark-rst-frontends")]
    {
        if std::env::var_os("RST_FRONTEND_AB_ONLY").is_some() {
            report_frontend_comparison(c);
            return;
        }
    }
    let cases = cases();
    let targets: Vec<_> = cases
        .iter()
        .map(|case| (*case, target_for(*case)))
        .collect();
    let mut evidence = Vec::new();
    for (case, target) in &targets {
        for _ in 0..REPORT_REPETITIONS {
            evidence.push(collect_case(*case, *target));
        }
    }
    print_report(&evidence);
    if std::env::var_os("RST_PREPARATION_REPORT_ONLY").is_some() {
        return;
    }

    let mut group = c.benchmark_group("equilibrium_rst_preparation");
    group.sample_size(10);
    group.warm_up_time(Duration::from_secs(1));
    group.measurement_time(Duration::from_secs(2));
    for (case, target) in targets {
        group.bench_function(BenchmarkId::from_parameter(case.label()), |bencher| {
            bencher.iter(|| black_box(solve_case(case, target)));
        });
    }
    group.finish();
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
