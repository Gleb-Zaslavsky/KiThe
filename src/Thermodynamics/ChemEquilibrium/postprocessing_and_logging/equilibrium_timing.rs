//! Optional stage timing for the canonical equilibrium workflow.
//!
//! Timing is deliberately a value object rather than logging side effects.
//! A caller can enable it for characterization or diagnostics, inspect the
//! immutable report attached to a solution, and leave the default production
//! path free from per-stage clock reads.

use std::time::{Duration, Instant};

/// Controls whether the canonical workflow records stage timings.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum EquilibriumTimingMode {
    /// Do not read the clock or populate timing fields.
    #[default]
    Disabled,
    /// Record the requested stage durations in the immutable result report.
    Enabled,
}

/// Effective symbolic backend reported by the RST preparation lifecycle.
///
/// `Mixed` can occur when one outer report aggregates several preparation
/// attempts. It is deliberately separate from a measured duration: backend
/// identity is metadata, not an additional timing stage.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum RstPreparationBackendKind {
    #[default]
    Unknown,
    Lambdify,
    Aot,
    Mixed,
}

impl RstPreparationBackendKind {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Unknown => "unknown",
            Self::Lambdify => "lambdify",
            Self::Aot => "aot",
            Self::Mixed => "mixed",
        }
    }
}

/// Artifact action reported by RST during symbolic preparation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum RstPreparationArtifactAction {
    #[default]
    Unknown,
    NotApplicable,
    Reused,
    Built,
    FallbackToLambdify,
    Mixed,
}

impl RstPreparationArtifactAction {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Unknown => "unknown",
            Self::NotApplicable => "not-applicable",
            Self::Reused => "reused",
            Self::Built => "built",
            Self::FallbackToLambdify => "fallback-to-lambdify",
            Self::Mixed => "mixed",
        }
    }
}

/// Fine-grained preparation evidence for one RST symbolic payload.
///
/// These stages are nested inside the public equilibrium timing report. They
/// are intentionally diagnostic: symbolic differentiation and lambdification
/// are owned by the RST adapter, while equation construction and parameter
/// evaluation remain KiThe-side preparation work.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct RstPreparationTimingReport {
    enabled: bool,
    /// Whether the RST dependency exposed a decomposition below aggregate
    /// problem materialization for this report.
    backend_stage_breakdown_available: bool,
    total: Duration,
    symbolic_input_validation: Duration,
    symbolic_phase_topology: Duration,
    symbolic_variable_expansion: Duration,
    symbolic_phase_totals: Duration,
    symbolic_activity_offsets: Duration,
    symbolic_reaction_assembly: Duration,
    symbolic_element_balance_assembly: Duration,
    symbolic_residual_simplification: Duration,
    equation_generation: Duration,
    temperature_substitution: Duration,
    row_scaling: Duration,
    enthalpy_equation: Duration,
    variable_materialization: Duration,
    parameter_evaluation: Duration,
    log_mole_bounds: Duration,
    symbolic_differentiation: Duration,
    residual_lambdification: Duration,
    jacobian_lambdification: Duration,
    problem_materialization: Duration,
    /// RST stages that may be inapplicable for the selected constructor.
    /// `None` remains distinct from a measured zero duration.
    rst_input_validation: Option<Duration>,
    rst_expression_parsing: Option<Duration>,
    rst_expression_graph_materialization: Option<Duration>,
    rst_parameter_binding: Option<Duration>,
    rst_problem_assembly: Option<Duration>,
    rst_aot_materialization: Option<Duration>,
    rst_final_validation: Option<Duration>,
    rst_unattributed: Option<Duration>,
    /// Effective backend metadata from RST's public preparation report.
    backend_kind: RstPreparationBackendKind,
    /// Artifact lifecycle metadata from RST's public preparation report.
    artifact_action: RstPreparationArtifactAction,
    /// Optional AOT build time. `None` means this lifecycle did not expose or
    /// perform a build, not that a build took zero time.
    backend_build_duration: Option<Duration>,
    /// Optional generated callback job counts, when RST prepared an AOT plan.
    generated_residual_jobs: Option<usize>,
    generated_jacobian_jobs: Option<usize>,
}

impl RstPreparationTimingReport {
    pub fn enabled(&self) -> bool {
        self.enabled
    }

    /// Returns whether symbolic differentiation and lambdification have
    /// independently measured durations.
    ///
    /// `false` means that the dependency returned only its aggregate
    /// preparation duration. When `true`, child stages below the aggregate
    /// were explicitly reported by RST and are safe to use for attribution.
    pub fn backend_stage_breakdown_available(&self) -> bool {
        self.backend_stage_breakdown_available
    }

    pub fn total(&self) -> Duration {
        self.total
    }

    pub fn symbolic_input_validation(&self) -> Duration {
        self.symbolic_input_validation
    }

    pub fn symbolic_phase_topology(&self) -> Duration {
        self.symbolic_phase_topology
    }

    pub fn symbolic_variable_expansion(&self) -> Duration {
        self.symbolic_variable_expansion
    }

    pub fn symbolic_phase_totals(&self) -> Duration {
        self.symbolic_phase_totals
    }

    pub fn symbolic_activity_offsets(&self) -> Duration {
        self.symbolic_activity_offsets
    }

    pub fn symbolic_reaction_assembly(&self) -> Duration {
        self.symbolic_reaction_assembly
    }

    pub fn symbolic_element_balance_assembly(&self) -> Duration {
        self.symbolic_element_balance_assembly
    }

    pub fn symbolic_residual_simplification(&self) -> Duration {
        self.symbolic_residual_simplification
    }

    pub fn equation_generation(&self) -> Duration {
        self.equation_generation
    }

    pub fn temperature_substitution(&self) -> Duration {
        self.temperature_substitution
    }

    pub fn row_scaling(&self) -> Duration {
        self.row_scaling
    }

    pub fn enthalpy_equation(&self) -> Duration {
        self.enthalpy_equation
    }

    pub fn variable_materialization(&self) -> Duration {
        self.variable_materialization
    }

    pub fn parameter_evaluation(&self) -> Duration {
        self.parameter_evaluation
    }

    pub fn log_mole_bounds(&self) -> Duration {
        self.log_mole_bounds
    }

    pub fn symbolic_differentiation(&self) -> Duration {
        self.symbolic_differentiation
    }

    pub fn residual_lambdification(&self) -> Duration {
        self.residual_lambdification
    }

    pub fn jacobian_lambdification(&self) -> Duration {
        self.jacobian_lambdification
    }

    pub fn problem_materialization(&self) -> Duration {
        self.problem_materialization
    }

    pub fn rst_input_validation(&self) -> Option<Duration> {
        self.rst_input_validation
    }

    pub fn rst_expression_parsing(&self) -> Option<Duration> {
        self.rst_expression_parsing
    }

    pub fn rst_expression_graph_materialization(&self) -> Option<Duration> {
        self.rst_expression_graph_materialization
    }

    pub fn rst_parameter_binding(&self) -> Option<Duration> {
        self.rst_parameter_binding
    }

    pub fn rst_problem_assembly(&self) -> Option<Duration> {
        self.rst_problem_assembly
    }

    pub fn rst_aot_materialization(&self) -> Option<Duration> {
        self.rst_aot_materialization
    }

    pub fn rst_final_validation(&self) -> Option<Duration> {
        self.rst_final_validation
    }

    pub fn rst_unattributed(&self) -> Option<Duration> {
        self.rst_unattributed
    }

    pub fn backend_kind(&self) -> RstPreparationBackendKind {
        self.backend_kind
    }

    pub fn artifact_action(&self) -> RstPreparationArtifactAction {
        self.artifact_action
    }

    pub fn backend_build_duration(&self) -> Option<Duration> {
        self.backend_build_duration
    }

    pub fn generated_residual_jobs(&self) -> Option<usize> {
        self.generated_residual_jobs
    }

    pub fn generated_jacobian_jobs(&self) -> Option<usize> {
        self.generated_jacobian_jobs
    }

    pub(crate) fn accumulate(&mut self, nested: Self) {
        if !nested.enabled {
            return;
        }
        self.enabled = true;
        self.backend_stage_breakdown_available |= nested.backend_stage_breakdown_available;
        self.total += nested.total;
        self.symbolic_input_validation += nested.symbolic_input_validation;
        self.symbolic_phase_topology += nested.symbolic_phase_topology;
        self.symbolic_variable_expansion += nested.symbolic_variable_expansion;
        self.symbolic_phase_totals += nested.symbolic_phase_totals;
        self.symbolic_activity_offsets += nested.symbolic_activity_offsets;
        self.symbolic_reaction_assembly += nested.symbolic_reaction_assembly;
        self.symbolic_element_balance_assembly += nested.symbolic_element_balance_assembly;
        self.symbolic_residual_simplification += nested.symbolic_residual_simplification;
        self.equation_generation += nested.equation_generation;
        self.temperature_substitution += nested.temperature_substitution;
        self.row_scaling += nested.row_scaling;
        self.enthalpy_equation += nested.enthalpy_equation;
        self.variable_materialization += nested.variable_materialization;
        self.parameter_evaluation += nested.parameter_evaluation;
        self.log_mole_bounds += nested.log_mole_bounds;
        self.symbolic_differentiation += nested.symbolic_differentiation;
        self.residual_lambdification += nested.residual_lambdification;
        self.jacobian_lambdification += nested.jacobian_lambdification;
        self.problem_materialization += nested.problem_materialization;
        self.rst_input_validation =
            add_optional_duration(self.rst_input_validation, nested.rst_input_validation);
        self.rst_expression_parsing =
            add_optional_duration(self.rst_expression_parsing, nested.rst_expression_parsing);
        self.rst_expression_graph_materialization = add_optional_duration(
            self.rst_expression_graph_materialization,
            nested.rst_expression_graph_materialization,
        );
        self.rst_parameter_binding =
            add_optional_duration(self.rst_parameter_binding, nested.rst_parameter_binding);
        self.rst_problem_assembly =
            add_optional_duration(self.rst_problem_assembly, nested.rst_problem_assembly);
        self.rst_aot_materialization =
            add_optional_duration(self.rst_aot_materialization, nested.rst_aot_materialization);
        self.rst_final_validation =
            add_optional_duration(self.rst_final_validation, nested.rst_final_validation);
        self.rst_unattributed =
            add_optional_duration(self.rst_unattributed, nested.rst_unattributed);
        self.backend_kind = merge_backend_kind(self.backend_kind, nested.backend_kind);
        self.artifact_action = merge_artifact_action(self.artifact_action, nested.artifact_action);
        self.backend_build_duration =
            add_optional_duration(self.backend_build_duration, nested.backend_build_duration);
        self.generated_residual_jobs =
            add_optional_count(self.generated_residual_jobs, nested.generated_residual_jobs);
        self.generated_jacobian_jobs =
            add_optional_count(self.generated_jacobian_jobs, nested.generated_jacobian_jobs);
    }

    fn record_backend_metadata(
        &mut self,
        backend_kind: RstPreparationBackendKind,
        artifact_action: RstPreparationArtifactAction,
        build_duration: Option<Duration>,
        generated_residual_jobs: Option<usize>,
        generated_jacobian_jobs: Option<usize>,
    ) {
        if !self.enabled {
            return;
        }
        self.backend_kind = merge_backend_kind(self.backend_kind, backend_kind);
        self.artifact_action = merge_artifact_action(self.artifact_action, artifact_action);
        self.backend_build_duration =
            add_optional_duration(self.backend_build_duration, build_duration);
        self.generated_residual_jobs =
            add_optional_count(self.generated_residual_jobs, generated_residual_jobs);
        self.generated_jacobian_jobs =
            add_optional_count(self.generated_jacobian_jobs, generated_jacobian_jobs);
    }
}

fn merge_backend_kind(
    current: RstPreparationBackendKind,
    next: RstPreparationBackendKind,
) -> RstPreparationBackendKind {
    match (current, next) {
        (RstPreparationBackendKind::Unknown, value)
        | (value, RstPreparationBackendKind::Unknown) => value,
        (left, right) if left == right => left,
        _ => RstPreparationBackendKind::Mixed,
    }
}

fn merge_artifact_action(
    current: RstPreparationArtifactAction,
    next: RstPreparationArtifactAction,
) -> RstPreparationArtifactAction {
    match (current, next) {
        (RstPreparationArtifactAction::Unknown, value)
        | (value, RstPreparationArtifactAction::Unknown) => value,
        (left, right) if left == right => left,
        _ => RstPreparationArtifactAction::Mixed,
    }
}

fn add_optional_duration(current: Option<Duration>, next: Option<Duration>) -> Option<Duration> {
    match (current, next) {
        (Some(left), Some(right)) => Some(left + right),
        (Some(value), None) | (None, Some(value)) => Some(value),
        (None, None) => None,
    }
}

fn add_optional_count(current: Option<usize>, next: Option<usize>) -> Option<usize> {
    match (current, next) {
        (Some(left), Some(right)) => Some(left.saturating_add(right)),
        (Some(value), None) | (None, Some(value)) => Some(value),
        (None, None) => None,
    }
}

#[allow(dead_code)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum RstPreparationTimingStage {
    SymbolicInputValidation,
    SymbolicPhaseTopology,
    SymbolicVariableExpansion,
    SymbolicPhaseTotals,
    SymbolicActivityOffsets,
    SymbolicReactionAssembly,
    SymbolicElementBalanceAssembly,
    SymbolicResidualSimplification,
    EquationGeneration,
    TemperatureSubstitution,
    RowScaling,
    EnthalpyEquation,
    VariableMaterialization,
    ParameterEvaluation,
    LogMoleBounds,
    SymbolicDifferentiation,
    ResidualLambdification,
    JacobianLambdification,
    ProblemMaterialization,
    RstInputValidation,
    RstExpressionParsing,
    RstExpressionGraphMaterialization,
    RstParameterBinding,
    RstProblemAssembly,
    RstAotMaterialization,
    RstFinalValidation,
    RstUnattributed,
}

#[derive(Debug)]
pub(crate) struct RstPreparationTimingCollector {
    report: RstPreparationTimingReport,
}

impl RstPreparationTimingCollector {
    pub(crate) fn new(mode: EquilibriumTimingMode) -> Self {
        Self {
            report: RstPreparationTimingReport {
                enabled: matches!(mode, EquilibriumTimingMode::Enabled),
                ..RstPreparationTimingReport::default()
            },
        }
    }

    pub(crate) fn measure<T, F>(&mut self, stage: RstPreparationTimingStage, operation: F) -> T
    where
        F: FnOnce() -> T,
    {
        if !self.report.enabled {
            return operation();
        }
        let started = Instant::now();
        let result = operation();
        self.record(stage, started.elapsed());
        result
    }

    pub(crate) fn record(&mut self, stage: RstPreparationTimingStage, duration: Duration) {
        let target = match stage {
            RstPreparationTimingStage::SymbolicInputValidation => {
                &mut self.report.symbolic_input_validation
            }
            RstPreparationTimingStage::SymbolicPhaseTopology => {
                &mut self.report.symbolic_phase_topology
            }
            RstPreparationTimingStage::SymbolicVariableExpansion => {
                &mut self.report.symbolic_variable_expansion
            }
            RstPreparationTimingStage::SymbolicPhaseTotals => {
                &mut self.report.symbolic_phase_totals
            }
            RstPreparationTimingStage::SymbolicActivityOffsets => {
                &mut self.report.symbolic_activity_offsets
            }
            RstPreparationTimingStage::SymbolicReactionAssembly => {
                &mut self.report.symbolic_reaction_assembly
            }
            RstPreparationTimingStage::SymbolicElementBalanceAssembly => {
                &mut self.report.symbolic_element_balance_assembly
            }
            RstPreparationTimingStage::SymbolicResidualSimplification => {
                &mut self.report.symbolic_residual_simplification
            }
            RstPreparationTimingStage::EquationGeneration => &mut self.report.equation_generation,
            RstPreparationTimingStage::TemperatureSubstitution => {
                &mut self.report.temperature_substitution
            }
            RstPreparationTimingStage::RowScaling => &mut self.report.row_scaling,
            RstPreparationTimingStage::EnthalpyEquation => &mut self.report.enthalpy_equation,
            RstPreparationTimingStage::VariableMaterialization => {
                &mut self.report.variable_materialization
            }
            RstPreparationTimingStage::ParameterEvaluation => &mut self.report.parameter_evaluation,
            RstPreparationTimingStage::LogMoleBounds => &mut self.report.log_mole_bounds,
            RstPreparationTimingStage::SymbolicDifferentiation => {
                &mut self.report.symbolic_differentiation
            }
            RstPreparationTimingStage::ResidualLambdification => {
                &mut self.report.residual_lambdification
            }
            RstPreparationTimingStage::JacobianLambdification => {
                &mut self.report.jacobian_lambdification
            }
            RstPreparationTimingStage::ProblemMaterialization => {
                &mut self.report.problem_materialization
            }
            RstPreparationTimingStage::RstInputValidation => {
                return self.record_optional_backend_stage(
                    RstPreparationTimingStage::RstInputValidation,
                    Some(duration),
                );
            }
            RstPreparationTimingStage::RstExpressionParsing => {
                return self.record_optional_backend_stage(
                    RstPreparationTimingStage::RstExpressionParsing,
                    Some(duration),
                );
            }
            RstPreparationTimingStage::RstExpressionGraphMaterialization => {
                return self.record_optional_backend_stage(
                    RstPreparationTimingStage::RstExpressionGraphMaterialization,
                    Some(duration),
                );
            }
            RstPreparationTimingStage::RstParameterBinding => {
                return self.record_optional_backend_stage(
                    RstPreparationTimingStage::RstParameterBinding,
                    Some(duration),
                );
            }
            RstPreparationTimingStage::RstProblemAssembly => {
                return self.record_optional_backend_stage(
                    RstPreparationTimingStage::RstProblemAssembly,
                    Some(duration),
                );
            }
            RstPreparationTimingStage::RstAotMaterialization => {
                return self.record_optional_backend_stage(
                    RstPreparationTimingStage::RstAotMaterialization,
                    Some(duration),
                );
            }
            RstPreparationTimingStage::RstFinalValidation => {
                return self.record_optional_backend_stage(
                    RstPreparationTimingStage::RstFinalValidation,
                    Some(duration),
                );
            }
            RstPreparationTimingStage::RstUnattributed => {
                return self.record_optional_backend_stage(
                    RstPreparationTimingStage::RstUnattributed,
                    Some(duration),
                );
            }
        };
        *target += duration;
    }

    pub(crate) fn mark_backend_stage_breakdown_available(&mut self) {
        if self.report.enabled {
            self.report.backend_stage_breakdown_available = true;
        }
    }

    pub(crate) fn record_optional_backend_stage(
        &mut self,
        stage: RstPreparationTimingStage,
        duration: Option<Duration>,
    ) {
        if !self.report.enabled {
            return;
        }
        self.report.backend_stage_breakdown_available = true;
        let target = match stage {
            RstPreparationTimingStage::RstInputValidation => &mut self.report.rst_input_validation,
            RstPreparationTimingStage::RstExpressionParsing => {
                &mut self.report.rst_expression_parsing
            }
            RstPreparationTimingStage::RstExpressionGraphMaterialization => {
                &mut self.report.rst_expression_graph_materialization
            }
            RstPreparationTimingStage::RstParameterBinding => {
                &mut self.report.rst_parameter_binding
            }
            RstPreparationTimingStage::RstProblemAssembly => &mut self.report.rst_problem_assembly,
            RstPreparationTimingStage::RstAotMaterialization => {
                &mut self.report.rst_aot_materialization
            }
            RstPreparationTimingStage::RstFinalValidation => &mut self.report.rst_final_validation,
            RstPreparationTimingStage::RstUnattributed => &mut self.report.rst_unattributed,
            _ => return,
        };
        *target = add_optional_duration(*target, duration);
    }

    pub(crate) fn accumulate(&mut self, nested: RstPreparationTimingReport) {
        self.report.accumulate(nested);
    }

    pub(crate) fn record_backend_metadata(
        &mut self,
        backend_kind: RstPreparationBackendKind,
        artifact_action: RstPreparationArtifactAction,
        build_duration: Option<Duration>,
        generated_residual_jobs: Option<usize>,
        generated_jacobian_jobs: Option<usize>,
    ) {
        self.report.record_backend_metadata(
            backend_kind,
            artifact_action,
            build_duration,
            generated_residual_jobs,
            generated_jacobian_jobs,
        );
    }

    pub(crate) fn set_total(&mut self, duration: Duration) {
        if self.report.enabled {
            self.report.total = duration;
        }
    }

    pub(crate) fn finish(self) -> RstPreparationTimingReport {
        self.report
    }
}

/// Named stages used by the canonical equilibrium pipeline.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum EquilibriumTimingStage {
    FacadeResolution,
    RequestPreparation,
    RouteSelection,
    Normalization,
    PhFormulationPreparation,
    Recovery,
    ResultAssembly,
    RepositoryLookup,
    ThermochemistryPreparation,
    NumericClosureConstruction,
    SymbolicConstruction,
    EquationConstruction,
    NumericalProblemPreparation,
    ActiveSetFormulationBuild,
    RstPreparation,
    ProjectionBuild,
    NonlinearSolve,
    PhaseControl,
    Validation,
    Postprocessing,
}

/// Immutable timing evidence attached to an accepted equilibrium result.
///
/// Durations are kept as `Duration` values so callers can choose their own
/// display unit. The `total` value is wall-clock time for the measured public
/// operation; stage values are inclusive where a stage contains nested solver
/// work and therefore should not be summed as an independent wall-clock total.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct EquilibriumTimingReport {
    enabled: bool,
    total: Duration,
    facade_resolution: Duration,
    request_preparation: Duration,
    route_selection: Duration,
    normalization: Duration,
    ph_formulation_preparation: Duration,
    recovery: Duration,
    result_assembly: Duration,
    repository_lookup: Duration,
    thermochemistry_preparation: Duration,
    numeric_closure_construction: Duration,
    symbolic_construction: Duration,
    equation_construction: Duration,
    numerical_problem_preparation: Duration,
    active_set_formulation_build: Duration,
    rst_preparation: Duration,
    projection_build: Duration,
    nonlinear_solve: Duration,
    phase_control: Duration,
    validation: Duration,
    postprocessing: Duration,
    rst_preparation_timing: RstPreparationTimingReport,
}

impl EquilibriumTimingReport {
    /// Whether this report contains measured values rather than zero defaults.
    pub fn enabled(&self) -> bool {
        self.enabled
    }

    /// Wall-clock time for the measured public operation.
    pub fn total(&self) -> Duration {
        self.total
    }

    pub fn facade_resolution(&self) -> Duration {
        self.facade_resolution
    }

    pub fn request_preparation(&self) -> Duration {
        self.request_preparation
    }

    pub fn route_selection(&self) -> Duration {
        self.route_selection
    }

    pub fn normalization(&self) -> Duration {
        self.normalization
    }

    pub fn ph_formulation_preparation(&self) -> Duration {
        self.ph_formulation_preparation
    }

    pub fn recovery(&self) -> Duration {
        self.recovery
    }

    pub fn result_assembly(&self) -> Duration {
        self.result_assembly
    }

    pub fn repository_lookup(&self) -> Duration {
        self.repository_lookup
    }

    pub fn thermochemistry_preparation(&self) -> Duration {
        self.thermochemistry_preparation
    }

    pub fn numeric_closure_construction(&self) -> Duration {
        self.numeric_closure_construction
    }

    pub fn symbolic_construction(&self) -> Duration {
        self.symbolic_construction
    }

    pub fn equation_construction(&self) -> Duration {
        self.equation_construction
    }

    pub fn numerical_problem_preparation(&self) -> Duration {
        self.numerical_problem_preparation
    }

    pub fn active_set_formulation_build(&self) -> Duration {
        self.active_set_formulation_build
    }

    pub fn rst_preparation(&self) -> Duration {
        self.rst_preparation
    }

    pub fn projection_build(&self) -> Duration {
        self.projection_build
    }

    pub fn nonlinear_solve(&self) -> Duration {
        self.nonlinear_solve
    }

    pub fn phase_control(&self) -> Duration {
        self.phase_control
    }

    pub fn validation(&self) -> Duration {
        self.validation
    }

    pub fn postprocessing(&self) -> Duration {
        self.postprocessing
    }

    /// Fine-grained stages inside RST symbolic/problem preparation.
    pub fn rst_preparation_timing(&self) -> RstPreparationTimingReport {
        self.rst_preparation_timing
    }

    /// Adds one nested timing report to an outer aggregate.
    ///
    /// Stage values are inclusive, so callers must interpret the aggregate as
    /// accumulated evidence rather than summing it with the outer wall clock.
    pub(crate) fn accumulate(&mut self, nested: Self) {
        if !nested.enabled {
            return;
        }
        if !self.enabled {
            *self = nested;
            return;
        }
        self.total += nested.total;
        self.facade_resolution += nested.facade_resolution;
        self.request_preparation += nested.request_preparation;
        self.route_selection += nested.route_selection;
        self.normalization += nested.normalization;
        self.ph_formulation_preparation += nested.ph_formulation_preparation;
        self.recovery += nested.recovery;
        self.result_assembly += nested.result_assembly;
        self.repository_lookup += nested.repository_lookup;
        self.thermochemistry_preparation += nested.thermochemistry_preparation;
        self.numeric_closure_construction += nested.numeric_closure_construction;
        self.symbolic_construction += nested.symbolic_construction;
        self.equation_construction += nested.equation_construction;
        self.numerical_problem_preparation += nested.numerical_problem_preparation;
        self.active_set_formulation_build += nested.active_set_formulation_build;
        self.rst_preparation += nested.rst_preparation;
        self.projection_build += nested.projection_build;
        self.nonlinear_solve += nested.nonlinear_solve;
        self.phase_control += nested.phase_control;
        self.validation += nested.validation;
        self.postprocessing += nested.postprocessing;
        self.rst_preparation_timing
            .accumulate(nested.rst_preparation_timing);
    }
}

/// Internal mutable collector used only while an immutable result is built.
#[derive(Debug)]
pub(crate) struct EquilibriumTimingCollector {
    report: EquilibriumTimingReport,
}

impl EquilibriumTimingCollector {
    pub(crate) fn new(mode: EquilibriumTimingMode) -> Self {
        Self {
            report: EquilibriumTimingReport {
                enabled: matches!(mode, EquilibriumTimingMode::Enabled),
                ..EquilibriumTimingReport::default()
            },
        }
    }

    pub(crate) fn from_report(report: EquilibriumTimingReport) -> Self {
        Self { report }
    }

    pub(crate) fn measure<T, F>(&mut self, stage: EquilibriumTimingStage, operation: F) -> T
    where
        F: FnOnce() -> T,
    {
        if !self.report.enabled {
            return operation();
        }

        let started = Instant::now();
        let result = operation();
        self.record(stage, started.elapsed());
        result
    }

    pub(crate) fn record(&mut self, stage: EquilibriumTimingStage, duration: Duration) {
        if !self.report.enabled {
            return;
        }

        let target = match stage {
            EquilibriumTimingStage::FacadeResolution => &mut self.report.facade_resolution,
            EquilibriumTimingStage::RequestPreparation => &mut self.report.request_preparation,
            EquilibriumTimingStage::RouteSelection => &mut self.report.route_selection,
            EquilibriumTimingStage::Normalization => &mut self.report.normalization,
            EquilibriumTimingStage::PhFormulationPreparation => {
                &mut self.report.ph_formulation_preparation
            }
            EquilibriumTimingStage::Recovery => &mut self.report.recovery,
            EquilibriumTimingStage::ResultAssembly => &mut self.report.result_assembly,
            EquilibriumTimingStage::RepositoryLookup => &mut self.report.repository_lookup,
            EquilibriumTimingStage::ThermochemistryPreparation => {
                &mut self.report.thermochemistry_preparation
            }
            EquilibriumTimingStage::NumericClosureConstruction => {
                &mut self.report.numeric_closure_construction
            }
            EquilibriumTimingStage::SymbolicConstruction => &mut self.report.symbolic_construction,
            EquilibriumTimingStage::EquationConstruction => &mut self.report.equation_construction,
            EquilibriumTimingStage::NumericalProblemPreparation => {
                &mut self.report.numerical_problem_preparation
            }
            EquilibriumTimingStage::ActiveSetFormulationBuild => {
                &mut self.report.active_set_formulation_build
            }
            EquilibriumTimingStage::RstPreparation => &mut self.report.rst_preparation,
            EquilibriumTimingStage::ProjectionBuild => &mut self.report.projection_build,
            EquilibriumTimingStage::NonlinearSolve => &mut self.report.nonlinear_solve,
            EquilibriumTimingStage::PhaseControl => &mut self.report.phase_control,
            EquilibriumTimingStage::Validation => &mut self.report.validation,
            EquilibriumTimingStage::Postprocessing => &mut self.report.postprocessing,
        };
        *target += duration;
    }

    pub(crate) fn record_rst_preparation_timing(&mut self, timing: RstPreparationTimingReport) {
        self.report.rst_preparation_timing.accumulate(timing);
    }

    /// Adds an already collected nested report to this in-progress report.
    pub(crate) fn accumulate(&mut self, nested: EquilibriumTimingReport) {
        self.report.accumulate(nested);
    }

    pub(crate) fn set_total(&mut self, duration: Duration) {
        if self.report.enabled {
            self.report.total = duration;
        }
    }

    pub(crate) fn finish(self) -> EquilibriumTimingReport {
        self.report
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::thread;

    #[test]
    fn disabled_collection_does_not_publish_measurements() {
        let mut collector = EquilibriumTimingCollector::new(EquilibriumTimingMode::Disabled);
        collector.measure(EquilibriumTimingStage::NonlinearSolve, || {
            thread::yield_now()
        });
        let report = collector.finish();

        assert!(!report.enabled());
        assert_eq!(report.nonlinear_solve(), Duration::ZERO);
    }

    #[test]
    fn enabled_collection_keeps_stage_and_total_evidence() {
        let mut collector = EquilibriumTimingCollector::new(EquilibriumTimingMode::Enabled);
        collector.measure(EquilibriumTimingStage::RouteSelection, || {
            thread::yield_now()
        });
        collector.measure(EquilibriumTimingStage::Normalization, || {
            thread::yield_now()
        });
        collector.measure(EquilibriumTimingStage::SymbolicConstruction, || {
            thread::yield_now()
        });
        collector.set_total(Duration::from_nanos(7));
        let report = collector.finish();

        assert!(report.enabled());
        assert!(report.route_selection() >= Duration::ZERO);
        assert!(report.normalization() >= Duration::ZERO);
        assert!(report.symbolic_construction() >= Duration::ZERO);
        assert_eq!(report.total(), Duration::from_nanos(7));
    }
}
