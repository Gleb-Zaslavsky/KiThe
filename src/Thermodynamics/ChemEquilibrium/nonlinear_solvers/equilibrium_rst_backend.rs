//! RustedSciThe adapters for the canonical log-moles equilibrium formulation.
//!
//! The adapter intentionally hands RST symbolic residual expressions rather
//! than KiThe's hand-written Jacobian. RST owns lambdification and symbolic
//! Jacobian construction, while this crate keeps the thermodynamic model and
//! the validation boundary. The prepared RST problem is wrapped in an opaque
//! adapter type so the domain layer does not depend on RST's concrete problem
//! representation.

use crate::Thermodynamics::ChemEquilibrium::equilibrium_log_moles::{EquilibriumLogMoles, Solvers};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_nonlinear::{
    BackendFailureKind, ReactionExtentError,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_ph_formulation::PreparedPhFormulation;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_problem::PreparedEquilibriumProblem;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_solver_policy::{
    SolverAttemptMetrics, SolverEvaluationTiming, SolverIterationTelemetry, SolverTermination,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_timing::{
    EquilibriumTimingCollector, EquilibriumTimingMode, EquilibriumTimingReport,
    EquilibriumTimingStage, RstPreparationArtifactAction, RstPreparationBackendKind,
    RstPreparationTimingCollector, RstPreparationTimingReport, RstPreparationTimingStage,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_workflows::{
    multiphase_equilibrium_residual_generator_sym,
    multiphase_equilibrium_residual_generator_sym_at_temperature_with_timing,
};
use crate::Thermodynamics::User_substances_error::SubsDataError;
use RustedSciThe::numerical::Nonlinear_systems::error::SolveError as RstSolveError;
use RustedSciThe::numerical::Nonlinear_systems::prelude::{
    BacktrackingLevenbergMarquardtMethod, BoundSymbolicNonlinearProblem, Bounds,
    DampedNewtonMethod, DiagnosticsOptions, LambdifyExecutionPolicy, LevenbergMarquardtMethod,
    LevenbergMarquardtMinpack, NielsenLevenbergMarquardtMethod, NonlinearSolverMethod,
    PowellDoglegMethod, PreparationStage, PreparationTelemetryMode,
    PreparedSymbolicNonlinearProblem, SolveOptions, SymbolicArtifactAction, SymbolicBackendKind,
    SymbolicLambdifyFrontend, SymbolicProblemOptions, TerminationReason, TrustRegionLMMethod,
};
#[cfg(test)]
use RustedSciThe::numerical::Nonlinear_systems::prelude::{
    JacobianProvider, NonlinearProblem, SolveResult,
};
use RustedSciThe::numerical::SolverTelemetryMode;
use RustedSciThe::symbolic::symbolic_engine::Expr;
#[cfg(test)]
use nalgebra::DMatrix;
use nalgebra::DVector;
use std::time::Instant;

/// Symbolic frontend used when preparing an RST nonlinear problem.
///
/// `ExprLegacy` remains the production default until AtomView's numerical
/// robustness envelope is fully characterized.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RstSymbolicFrontend {
    ExprLegacy,
    AtomView,
}

impl From<RstSymbolicFrontend> for SymbolicLambdifyFrontend {
    fn from(frontend: RstSymbolicFrontend) -> Self {
        match frontend {
            RstSymbolicFrontend::ExprLegacy => Self::ExprLegacy,
            RstSymbolicFrontend::AtomView => Self::AtomViewNative,
        }
    }
}

#[cfg(any(test, feature = "equilibrium-benchmark-rst-frontends"))]
use std::cell::RefCell;

#[cfg(any(test, feature = "equilibrium-benchmark-rst-frontends"))]
thread_local! {
    static RST_TEST_FRONTEND_OVERRIDE: RefCell<Option<SymbolicLambdifyFrontend>> = const { RefCell::new(None) };
}

#[cfg(test)]
thread_local! {
    static RST_TEST_HISTORY_CAPTURE_ENABLED: RefCell<bool> = const { RefCell::new(false) };
    static RST_TEST_CAPTURED_RESULTS: RefCell<Vec<SolveResult>> = const { RefCell::new(Vec::new()) };
}

#[cfg(test)]
pub(crate) fn with_test_rst_frontend<T>(
    frontend: SymbolicLambdifyFrontend,
    operation: impl FnOnce() -> T,
) -> T {
    with_rst_frontend_override(frontend, operation)
}

#[cfg(feature = "equilibrium-benchmark-rst-frontends")]
#[doc(hidden)]
pub fn with_benchmark_rst_frontend<T>(
    frontend: SymbolicLambdifyFrontend,
    operation: impl FnOnce() -> T,
) -> T {
    with_rst_frontend_override(frontend, operation)
}

#[cfg(any(test, feature = "equilibrium-benchmark-rst-frontends"))]
fn with_rst_frontend_override<T>(
    frontend: SymbolicLambdifyFrontend,
    operation: impl FnOnce() -> T,
) -> T {
    struct RestoreFrontend(Option<SymbolicLambdifyFrontend>);

    impl Drop for RestoreFrontend {
        fn drop(&mut self) {
            RST_TEST_FRONTEND_OVERRIDE.with(|current| {
                *current.borrow_mut() = self.0.take();
            });
        }
    }

    let previous =
        RST_TEST_FRONTEND_OVERRIDE.with(|current| current.borrow_mut().replace(frontend));
    let _restore = RestoreFrontend(previous);
    operation()
}

#[cfg(test)]
fn rst_test_history_capture_enabled() -> bool {
    RST_TEST_HISTORY_CAPTURE_ENABLED.with(|enabled| *enabled.borrow())
}

/// Immutable symbolic thermochemistry snapshot prepared for the RST adapter.
///
/// The snapshot keeps the species ordering and the associated symbolic
/// standard Gibbs expressions together so the backend bridge can construct its
/// symbolic problem without reopening the lookup pipeline.
pub(crate) struct RstSymbolicThermochemistry {
    /// Symbolic standard Gibbs free energy expressions `G_i^0(T)` for each species.
    /// These are `Expr` values from RustedSciThe's symbolic engine, used to build
    /// the symbolic residual and Jacobian for automatic differentiation.
    standard_gibbs: Vec<Expr>,
}

impl RstSymbolicThermochemistry {
    /// Builds an immutable symbolic thermochemistry snapshot.
    ///
    /// A phase-system bridge may inject phase-qualified symbolic `G0(T)`
    /// expressions directly into `gibbs_sym`. That is the preferred path: it
    /// preserves the already validated component order and does not reopen the
    /// mutable `SubsData` lookup pipeline. The historical facade remains a
    /// fallback for callers that still construct an `EquilibriumLogMoles`
    /// object directly.
    pub(crate) fn from_solver(
        solver: &mut EquilibriumLogMoles,
    ) -> Result<Self, ReactionExtentError> {
        if solver.gibbs_sym.len() == solver.subs_data.substances.len()
            && !solver.gibbs_sym.is_empty()
        {
            return Ok(Self {
                standard_gibbs: solver.gibbs_sym.clone(),
            });
        }
        if !solver.gibbs_sym.is_empty() {
            return Err(ReactionExtentError::DimensionMismatch(format!(
                "symbolic Gibbs snapshot has {} entries for {} species",
                solver.gibbs_sym.len(),
                solver.subs_data.substances.len()
            )));
        }

        solver
            .subs_data
            .calculate_therm_map_of_sym()
            .map_err(ReactionExtentError::SubsDataError)?;
        let gibbs_by_species = solver.subs_data.calculate_dG0_sym_one_phase()?;
        let mut standard_gibbs = Vec::with_capacity(solver.subs_data.substances.len());

        for species in &solver.subs_data.substances {
            let expression = gibbs_by_species.get(species).ok_or_else(|| {
                ReactionExtentError::SubsDataError(SubsDataError::MissingData {
                    field: "symbolic standard Gibbs expression".to_string(),
                    substance: species.clone(),
                })
            })?;
            standard_gibbs.push(expression.clone());
        }

        Ok(Self { standard_gibbs })
    }
}

/// Typed solve contract passed from the equilibrium layer into RST.
///
/// The contract intentionally carries only the numeric semantics we actually
/// want to promise across the boundary: a convergence tolerance and a maximum
/// nonlinear-iteration budget. Feasibility bounds are a separate domain-side
/// concern in the log-mole formulation, so they do not belong in the symbolic
/// RST adapter contract.
#[derive(Debug, Clone, PartialEq)]
pub struct RustedSciTheSolveContract {
    /// Convergence tolerance forwarded to RST.
    pub tolerance: f64,
    /// Maximum nonlinear iterations forwarded to RST.
    pub max_iterations: usize,
    /// RST work-collection level. Ordinary solves keep this disabled.
    pub telemetry_mode: SolverTelemetryMode,
    /// Optional finite box bounds in the backend coordinate order.
    ///
    /// Canonical fixed-`P,T` log-mole problems provide these from conserved
    /// elemental inventory. Compatibility and coupled `P,H` paths leave them
    /// unset until they can construct bounds for every coordinate safely.
    log_mole_bounds: Option<Vec<(f64, f64)>>,
}

impl RustedSciTheSolveContract {
    /// Builds a validated solve contract from the numeric values that the
    /// equilibrium layer wants to pass to RST.
    pub fn new(tolerance: f64, max_iterations: usize) -> Result<Self, ReactionExtentError> {
        if !tolerance.is_finite() || tolerance <= 0.0 {
            return Err(ReactionExtentError::InvalidProblem {
                field: "rst_solve_contract.tolerance",
                message: "tolerance must be finite and strictly positive".to_string(),
            });
        }
        if max_iterations == 0 {
            return Err(ReactionExtentError::InvalidProblem {
                field: "rst_solve_contract.max_iterations",
                message: "max_iterations must be greater than zero".to_string(),
            });
        }
        Ok(Self {
            tolerance,
            max_iterations,
            telemetry_mode: SolverTelemetryMode::Off,
            log_mole_bounds: None,
        })
    }

    /// Selects which native RST statistics are collected for this attempt.
    pub fn with_telemetry_mode(mut self, mode: SolverTelemetryMode) -> Self {
        self.telemetry_mode = mode;
        self
    }

    /// Attaches prevalidated finite bounds for a fixed log-mole formulation.
    pub(crate) fn with_log_mole_bounds(
        mut self,
        bounds: &[(f64, f64)],
    ) -> Result<Self, ReactionExtentError> {
        if bounds.is_empty() {
            return Err(ReactionExtentError::InvalidProblem {
                field: "rst_log_mole_bounds",
                message: "bounds must contain one entry per nonlinear coordinate".to_string(),
            });
        }
        for (index, &(lower, upper)) in bounds.iter().enumerate() {
            if !lower.is_finite() || !upper.is_finite() || lower > upper {
                return Err(ReactionExtentError::InvalidProblem {
                    field: "rst_log_mole_bounds",
                    message: format!("bound {index} must be finite and satisfy lower <= upper"),
                });
            }
        }
        self.log_mole_bounds = Some(bounds.to_vec());
        Ok(self)
    }

    /// Converts the typed contract into the concrete RST options struct.
    pub fn to_options(self) -> Result<SolveOptions, ReactionExtentError> {
        let bounds = self
            .log_mole_bounds
            .map(Bounds::new)
            .transpose()
            .map_err(|error| ReactionExtentError::InvalidProblem {
                field: "rst_log_mole_bounds",
                message: error.to_string(),
            })?;
        Ok(SolveOptions {
            tolerance: self.tolerance,
            max_iterations: self.max_iterations,
            bounds,
            diagnostics: DiagnosticsOptions {
                telemetry_mode: self.telemetry_mode,
                collect_history: false,
                ..DiagnosticsOptions::default()
            },
            ..SolveOptions::default()
        })
    }
}

/// Local prepared symbolic payload owned by the RST adapter layer.
///
/// The immutable RST prepared problem owns symbolic differentiation and
/// lambdification exactly once. Runtime thermochemistry/target changes are
/// represented by a short numeric parameter vector and bound views created at
/// solve time; they never rebuild this symbolic payload.
pub(crate) struct RstSymbolicProblem {
    prepared: PreparedSymbolicNonlinearProblem,
    equation_parameters: Option<Vec<String>>,
    equation_parameter_values: Option<DVector<f64>>,
}

impl RstSymbolicProblem {
    pub(crate) fn from_expressions(
        equations: Vec<Expr>,
        variables: Vec<String>,
        equation_parameters: Option<Vec<String>>,
        equation_parameter_values: Option<DVector<f64>>,
    ) -> Result<Self, ReactionExtentError> {
        Self::from_expressions_with_timing(
            equations,
            variables,
            equation_parameters,
            equation_parameter_values,
            EquilibriumTimingMode::Disabled,
        )
        .map(|(problem, _)| problem)
    }

    pub(crate) fn from_expressions_with_timing(
        equations: Vec<Expr>,
        variables: Vec<String>,
        equation_parameters: Option<Vec<String>>,
        equation_parameter_values: Option<DVector<f64>>,
        timing_mode: EquilibriumTimingMode,
    ) -> Result<(Self, RstPreparationTimingReport), ReactionExtentError> {
        Self::from_expressions_with_timing_and_frontend(
            equations,
            variables,
            equation_parameters,
            equation_parameter_values,
            timing_mode,
            // AtomView is available as an explicit A/B route, but remains
            // opt-in until the real TP-1907 acceptance matrix passes it.
            SymbolicLambdifyFrontend::ExprLegacy,
        )
    }

    pub(crate) fn from_expressions_with_timing_and_frontend(
        equations: Vec<Expr>,
        variables: Vec<String>,
        equation_parameters: Option<Vec<String>>,
        equation_parameter_values: Option<DVector<f64>>,
        timing_mode: EquilibriumTimingMode,
        frontend: SymbolicLambdifyFrontend,
    ) -> Result<(Self, RstPreparationTimingReport), ReactionExtentError> {
        let started = Instant::now();
        let mut timing = RstPreparationTimingCollector::new(timing_mode);
        if equations.is_empty() {
            return Err(ReactionExtentError::InvalidProblem {
                field: "rst_symbolic_equations",
                message: "equation system must not be empty".to_string(),
            });
        }
        if variables.is_empty() {
            return Err(ReactionExtentError::InvalidProblem {
                field: "rst_symbolic_variables",
                message: "symbolic variable list must not be empty".to_string(),
            });
        }
        if variables.len() != equations.len() {
            return Err(ReactionExtentError::DimensionMismatch(format!(
                "RST symbolic variables vs equations: expected {}, got {}",
                equations.len(),
                variables.len()
            )));
        }
        if let Some(values) = equation_parameter_values.as_ref() {
            if equation_parameters.as_ref().map_or(0, Vec::len) != values.len() {
                return Err(ReactionExtentError::DimensionMismatch(format!(
                    "RST symbolic parameter values: expected {}, got {}",
                    equation_parameters.as_ref().map_or(0, Vec::len),
                    values.len()
                )));
            }
        }

        let options = SymbolicProblemOptions::new()
            .with_variables(variables.clone())
            // RST keeps the parallel policy available, but for KiThe's
            // current problem sizes Rayon scheduling costs more than it saves.
            // Sequential is therefore the production default; range/session
            // reuse amortizes this preparation over many equilibrium points.
            .with_lambdify_execution_policy(LambdifyExecutionPolicy::Sequential)
            .with_lambdify_frontend(frontend)
            .with_lambdify_backend()
            .with_preparation_telemetry(match timing_mode {
                EquilibriumTimingMode::Enabled => PreparationTelemetryMode::Collect,
                EquilibriumTimingMode::Disabled => PreparationTelemetryMode::Disabled,
            });
        let options = match equation_parameters.as_ref() {
            Some(parameters) => {
                let options = options.with_equation_parameters(parameters.clone());
                match equation_parameter_values.as_ref() {
                    Some(values) => options.with_equation_parameter_values(values.clone()),
                    None => options,
                }
            }
            None => options,
        };
        let prepared = PreparedSymbolicNonlinearProblem::from_expressions(equations, options)
            .map_err(|error| ReactionExtentError::InvalidProblem {
                field: "rst_symbolic_problem",
                message: format!("failed to prepare RST symbolic problem: {error}"),
            })?;
        if matches!(timing_mode, EquilibriumTimingMode::Enabled) {
            let report = prepared.preparation_report();
            timing.record_backend_metadata(
                match report.effective_backend {
                    SymbolicBackendKind::Lambdify => RstPreparationBackendKind::Lambdify,
                    SymbolicBackendKind::Aot => RstPreparationBackendKind::Aot,
                },
                match report.artifact_action {
                    SymbolicArtifactAction::NotApplicable => {
                        RstPreparationArtifactAction::NotApplicable
                    }
                    SymbolicArtifactAction::Reused => RstPreparationArtifactAction::Reused,
                    SymbolicArtifactAction::Built => RstPreparationArtifactAction::Built,
                    SymbolicArtifactAction::FallbackToLambdify => {
                        RstPreparationArtifactAction::FallbackToLambdify
                    }
                },
                report.build_duration,
                report.generated_residual_jobs,
                report.generated_jacobian_jobs,
            );
            timing.record(
                RstPreparationTimingStage::ProblemMaterialization,
                report.preparation_duration,
            );
            if let Some(detailed) = report.detailed.as_ref() {
                timing.mark_backend_stage_breakdown_available();
                record_rst_preparation_stage(
                    &mut timing,
                    RstPreparationTimingStage::RstInputValidation,
                    detailed,
                    PreparationStage::InputValidation,
                );
                record_rst_preparation_stage(
                    &mut timing,
                    RstPreparationTimingStage::RstExpressionParsing,
                    detailed,
                    PreparationStage::ExpressionParsing,
                );
                record_rst_preparation_stage(
                    &mut timing,
                    RstPreparationTimingStage::RstExpressionGraphMaterialization,
                    detailed,
                    PreparationStage::ExpressionGraphMaterialization,
                );
                record_rst_preparation_stage(
                    &mut timing,
                    RstPreparationTimingStage::SymbolicDifferentiation,
                    detailed,
                    match frontend {
                        SymbolicLambdifyFrontend::ExprLegacy => {
                            PreparationStage::JacobianDifferentiation
                        }
                        SymbolicLambdifyFrontend::AtomViewNative => {
                            PreparationStage::AtomDifferentiation
                        }
                    },
                );
                record_rst_preparation_stage(
                    &mut timing,
                    RstPreparationTimingStage::ResidualLambdification,
                    detailed,
                    PreparationStage::ResidualCallbackPreparation,
                );
                record_rst_preparation_stage(
                    &mut timing,
                    RstPreparationTimingStage::JacobianLambdification,
                    detailed,
                    PreparationStage::JacobianCallbackPreparation,
                );
                record_rst_preparation_stage(
                    &mut timing,
                    RstPreparationTimingStage::RstParameterBinding,
                    detailed,
                    PreparationStage::ParameterBinding,
                );
                record_rst_preparation_stage(
                    &mut timing,
                    RstPreparationTimingStage::RstProblemAssembly,
                    detailed,
                    PreparationStage::PreparedProblemAssembly,
                );
                record_rst_preparation_stage(
                    &mut timing,
                    RstPreparationTimingStage::RstAotMaterialization,
                    detailed,
                    PreparationStage::AotMaterialization,
                );
                record_rst_preparation_stage(
                    &mut timing,
                    RstPreparationTimingStage::RstFinalValidation,
                    detailed,
                    PreparationStage::FinalValidation,
                );
                timing.record_optional_backend_stage(
                    RstPreparationTimingStage::RstUnattributed,
                    Some(detailed.unattributed_wall_time),
                );
            }
        }
        timing.set_total(started.elapsed());
        Ok((
            Self {
                prepared,
                equation_parameters,
                equation_parameter_values,
            },
            timing.finish(),
        ))
    }

    fn set_parameter_values(&mut self, values: DVector<f64>) -> Result<(), ReactionExtentError> {
        if let Some(parameters) = self.equation_parameters.as_ref() {
            if parameters.len() != values.len() {
                return Err(ReactionExtentError::DimensionMismatch(format!(
                    "RST symbolic equation parameters: expected {}, got {}",
                    parameters.len(),
                    values.len()
                )));
            }
        }
        self.equation_parameter_values = Some(values);
        Ok(())
    }

    fn bind(&self) -> Result<BoundSymbolicNonlinearProblem<'_>, ReactionExtentError> {
        let result = match self.equation_parameters.as_ref() {
            Some(_) => self
                .equation_parameter_values
                .as_ref()
                .ok_or_else(|| ReactionExtentError::InvalidProblem {
                    field: "rst_symbolic_parameter_values",
                    message: "parameterized RST problem has no bound numeric values".to_string(),
                })
                .and_then(|values| {
                    self.prepared.bind_values(values.clone()).map_err(|error| {
                        ReactionExtentError::InvalidProblem {
                            field: "rst_symbolic_parameter_values",
                            message: format!("failed to bind RST numeric parameters: {error}"),
                        }
                    })
                }),
            None => self.prepared.bind_without_parameters().map_err(|error| {
                ReactionExtentError::InvalidProblem {
                    field: "rst_symbolic_parameter_values",
                    message: format!("failed to bind non-parameterized RST problem: {error}"),
                }
            }),
        };
        result
    }
}

fn record_rst_preparation_stage(
    timing: &mut RstPreparationTimingCollector,
    target: RstPreparationTimingStage,
    detailed: &RustedSciThe::numerical::Nonlinear_systems::prelude::PreparationTelemetry,
    source: PreparationStage,
) {
    let duration = detailed.stage(source).and_then(|stage| stage.wall_time);
    match target {
        RstPreparationTimingStage::SymbolicDifferentiation
        | RstPreparationTimingStage::ResidualLambdification
        | RstPreparationTimingStage::JacobianLambdification => {
            if let Some(duration) = duration {
                timing.record(target, duration);
            }
        }
        _ => timing.record_optional_backend_stage(target, duration),
    }
}

/// Opaque prepared symbolic problem owned by the RST adapter layer.
///
/// The domain layer may carry this value around as a prepared backend request,
/// but it should never inspect or construct RST internals directly.
pub(crate) struct RstPreparedProblem {
    problem: RstSymbolicProblem,
    log_mole_bounds: Option<Vec<(f64, f64)>>,
    fixed_pt_species_count: Option<usize>,
}

impl RstPreparedProblem {
    pub(crate) fn new(problem: RstSymbolicProblem) -> Self {
        Self {
            problem,
            log_mole_bounds: None,
            fixed_pt_species_count: None,
        }
    }

    /// Test-only export of the prepared symbolic payload for an external RST
    /// reproducer. Production callers continue to see this payload as opaque.
    #[cfg(test)]
    pub(crate) fn symbolic_payload_for_test(
        &self,
    ) -> (Vec<String>, Vec<String>, Vec<f64>, Vec<String>) {
        (
            self.problem
                .prepared
                .as_problem()
                .equations()
                .iter()
                .map(ToString::to_string)
                .collect(),
            self.problem.prepared.variables().to_vec(),
            self.problem
                .equation_parameter_values
                .as_ref()
                .map(|values| values.as_slice().to_vec())
                .unwrap_or_default(),
            self.problem
                .equation_parameters
                .as_ref()
                .cloned()
                .unwrap_or_default(),
        )
    }

    #[cfg(test)]
    pub(crate) fn symbolic_equations_for_test(&self) -> Vec<Expr> {
        self.problem.prepared.as_problem().equations().to_vec()
    }

    /// Test-only export of the symbolic Jacobian retained by RST preparation.
    #[cfg(test)]
    pub(crate) fn symbolic_jacobian_for_test(&self) -> Vec<Vec<String>> {
        self.problem
            .prepared
            .as_problem()
            .symbolic_jacobian()
            .iter()
            .map(|row| row.iter().map(ToString::to_string).collect())
            .collect()
    }

    /// Test-only export of the bounds used by the production RST adapter.
    #[cfg(test)]
    pub(crate) fn log_mole_bounds_for_test(&self) -> Option<Vec<(f64, f64)>> {
        self.log_mole_bounds.as_deref().map(ToOwned::to_owned)
    }

    fn with_fixed_pt_parameters(
        problem: RstSymbolicProblem,
        log_mole_bounds: Vec<(f64, f64)>,
        species_count: usize,
    ) -> Self {
        Self {
            problem,
            log_mole_bounds: Some(log_mole_bounds),
            fixed_pt_species_count: Some(species_count),
        }
    }

    /// Retargets the reusable fixed-`P,T` symbolic graph without rebuilding it.
    ///
    /// `G_i^0(T)` enters the reaction-affinity rows as an equation parameter.
    /// This keeps the residual/Jacobian graph stable across NASA/NIST interval
    /// boundaries: only the already validated numeric thermochemistry changes.
    pub(crate) fn set_fixed_pt_thermochemistry(
        &mut self,
        temperature: f64,
        standard_gibbs: &[f64],
    ) -> Result<(), ReactionExtentError> {
        let species_count =
            self.fixed_pt_species_count
                .ok_or_else(|| {
                    ReactionExtentError::InvalidProblem {
                field: "rst_fixed_pt_parameters",
                message:
                    "this RST problem was not prepared for parameterized fixed-P,T thermochemistry"
                        .to_string(),
            }
                })?;
        if !temperature.is_finite()
            || standard_gibbs.len() != species_count
            || standard_gibbs.iter().any(|value| !value.is_finite())
        {
            return Err(ReactionExtentError::InvalidProblem {
                field: "rst_fixed_pt_parameters",
                message: format!(
                    "temperature must be finite and {species_count} finite standard Gibbs values are required"
                ),
            });
        }
        let mut parameters = Vec::with_capacity(species_count + 1);
        parameters.push(temperature);
        parameters.extend_from_slice(standard_gibbs);
        self.problem
            .set_parameter_values(DVector::from_vec(parameters))
            .map_err(|error| {
                ReactionExtentError::ResidualEvaluation(format!(
                    "failed to update RST fixed-P,T thermochemistry parameters: {error}"
                ))
            })
    }

    /// Evaluates a prepared fixed-`P,T` problem's standard-state closures and
    /// updates the parameterized RST graph in one checked operation.
    pub(crate) fn set_fixed_pt_thermochemistry_from_prepared(
        &mut self,
        prepared: &PreparedEquilibriumProblem,
    ) -> Result<(), ReactionExtentError> {
        let temperature = prepared.problem().conditions().temperature();
        let standard_gibbs = standard_gibbs_at_conditions(prepared)?;
        self.set_fixed_pt_thermochemistry(temperature, &standard_gibbs)
    }

    /// Updates the target and scale parameters of a prepared monolithic P,H
    /// problem without rebuilding symbolic expressions or lambdified code.
    pub(crate) fn set_ph_parameters(
        &mut self,
        target_enthalpy: f64,
        enthalpy_scale: f64,
    ) -> Result<(), ReactionExtentError> {
        if !target_enthalpy.is_finite() || !enthalpy_scale.is_finite() || enthalpy_scale <= 0.0 {
            return Err(ReactionExtentError::InvalidProblem {
                field: "ph_symbolic_parameters",
                message: "P,H symbolic target and scale must be finite; scale must be positive"
                    .to_string(),
            });
        }
        self.problem
            .set_parameter_values(DVector::from_vec(vec![target_enthalpy, enthalpy_scale]))
            .map_err(|error| {
                ReactionExtentError::ResidualEvaluation(format!(
                    "failed to update RST P,H parameters: {error}"
                ))
            })
    }

    fn bound_problem(&self) -> Result<BoundSymbolicNonlinearProblem<'_>, ReactionExtentError> {
        self.problem.bind()
    }

    fn log_mole_bounds(&self) -> Option<&[(f64, f64)]> {
        self.log_mole_bounds.as_deref()
    }

    /// Test-only inspection point for proving that parameterized fixed-`P,T`
    /// graphs remain algebraically identical to their baked-symbolic form.
    #[cfg(test)]
    pub(crate) fn residual_for_test(
        &self,
        log_moles: &[f64],
    ) -> Result<DVector<f64>, ReactionExtentError> {
        let bound = self.problem.bind()?;
        bound
            .residual(&DVector::from_column_slice(log_moles))
            .map_err(|error| {
                ReactionExtentError::ResidualEvaluation(format!(
                    "test evaluation of parameterized RST residual failed: {error}"
                ))
            })
    }

    /// Test-only inspection point paired with [`Self::residual_for_test`].
    #[cfg(test)]
    pub(crate) fn jacobian_for_test(
        &self,
        log_moles: &[f64],
    ) -> Result<DMatrix<f64>, ReactionExtentError> {
        let bound = self.problem.bind()?;
        bound
            .jacobian(&DVector::from_column_slice(log_moles))
            .map_err(|error| {
                ReactionExtentError::ResidualEvaluation(format!(
                    "test evaluation of parameterized RST Jacobian failed: {error}"
                ))
            })
    }

    #[cfg(test)]
    pub(crate) fn residual_components_for_test(
        &self,
        log_moles: &[f64],
    ) -> (Result<(), ReactionExtentError>, Vec<f64>) {
        let bound = match self.problem.bind() {
            Ok(bound) => bound,
            Err(error) => return (Err(error.into()), Vec::new()),
        };
        let mut output = DVector::zeros(bound.dimension());
        let result = bound
            .residual_into(&DVector::from_column_slice(log_moles), &mut output)
            .map_err(|error| {
                ReactionExtentError::ResidualEvaluation(format!(
                    "test evaluation of parameterized RST residual components failed: {error}"
                ))
            });
        (result, output.as_slice().to_vec())
    }
}

/// Candidate and engine diagnostics returned by one RST solve attempt.
#[derive(Debug, Clone)]
pub struct RustedSciTheSolveOutcome {
    /// Final log-mole iterate returned by RST.
    pub solution: Vec<f64>,
    /// Solver-engine counters and termination reason.
    pub metrics: SolverAttemptMetrics,
}

#[cfg(test)]
pub(crate) struct RstSolveFailureTrace {
    pub error: ReactionExtentError,
    pub residual_inputs: Vec<Vec<f64>>,
    pub jacobian_inputs: Vec<Vec<f64>>,
}

/// RST nonlinear methods suitable for the log-moles equilibrium system.
///
/// The ordered production default starts with ordinary LM as a transparent
/// baseline, then tries guarded least-squares/trust-region methods and
/// finally the backtracking LM fallback before handing control to Legacy.
/// Individual methods remain selectable even when they are deliberately
/// omitted from that default.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RustedSciTheSolver {
    /// Classical Levenberg-Marquardt.
    LevenbergMarquardt,
    /// Identity-damped Levenberg-Marquardt with residual-decrease backtracking.
    ///
    /// This is the RST counterpart of the robust legacy LM route and remains
    /// distinct from the canonical and MINPACK LM implementations.
    BacktrackingLevenbergMarquardt,
    /// MINPACK-style Levenberg-Marquardt with trust-region control.
    MinpackLevenbergMarquardt,
    /// Nielsen's adaptive-damping Levenberg-Marquardt method.
    NielsenLevenbergMarquardt,
    /// MINPACK-style trust-region Levenberg-Marquardt.
    TrustRegionLevenbergMarquardt,
    /// Powell dogleg trust-region method.
    ///
    /// This remains available through an explicit solver policy, but is
    /// excluded from [`Self::recommended_cascade`]: the currently supported
    /// RustedSciThe version writes unconditional diagnostics to stdout and its
    /// large real-equilibrium characterization produced rejected candidates.
    PowellDogleg,
    /// Newton method protected by backtracking line search.
    DampedNewton,
}

impl RustedSciTheSolver {
    /// Quiet, robust default order for production equilibrium solves.
    ///
    /// The full backend matrix is still exercised by ignored release tests.
    /// A benchmark-only method must not make a normal library or GUI solve
    /// emit unstructured dependency output.
    pub fn recommended_cascade() -> Vec<Self> {
        vec![
            Self::LevenbergMarquardt,
            Self::MinpackLevenbergMarquardt,
            Self::NielsenLevenbergMarquardt,
            Self::TrustRegionLevenbergMarquardt,
            Self::DampedNewton,
            Self::BacktrackingLevenbergMarquardt,
        ]
    }

    /// Stable identifier used in reports and tests.
    pub fn name(self) -> &'static str {
        match self {
            Self::LevenbergMarquardt => "rst_levenberg_marquardt",
            Self::BacktrackingLevenbergMarquardt => "rst_backtracking_levenberg_marquardt",
            Self::MinpackLevenbergMarquardt => "rst_minpack_levenberg_marquardt",
            Self::NielsenLevenbergMarquardt => "rst_nielsen_levenberg_marquardt",
            Self::TrustRegionLevenbergMarquardt => "rst_trust_region_levenberg_marquardt",
            Self::PowellDogleg => "rst_powell_dogleg",
            Self::DampedNewton => "rst_damped_newton",
        }
    }

    /// Runs one RST method. Residual and Jacobian evaluation remain entirely
    /// inside the prepared RST symbolic problem.
    pub(crate) fn solve(
        self,
        problem: &RstPreparedProblem,
        initial_log_moles: &[f64],
        contract: &RustedSciTheSolveContract,
    ) -> Result<RustedSciTheSolveOutcome, ReactionExtentError> {
        let method = match self {
            Self::LevenbergMarquardt => {
                NonlinearSolverMethod::LevenbergMarquardt(LevenbergMarquardtMethod::default())
            }
            Self::BacktrackingLevenbergMarquardt => {
                NonlinearSolverMethod::BacktrackingLevenbergMarquardt(
                    BacktrackingLevenbergMarquardtMethod::default(),
                )
            }
            Self::MinpackLevenbergMarquardt => NonlinearSolverMethod::LevenbergMarquardtMinpack(
                LevenbergMarquardtMinpack::default(),
            ),
            Self::NielsenLevenbergMarquardt => NonlinearSolverMethod::NielsenLevenbergMarquardt(
                NielsenLevenbergMarquardtMethod::default(),
            ),
            Self::TrustRegionLevenbergMarquardt => {
                NonlinearSolverMethod::TrustRegionLM(TrustRegionLMMethod::default())
            }
            Self::PowellDogleg => {
                NonlinearSolverMethod::PowellDogleg(PowellDoglegMethod::default())
            }
            Self::DampedNewton => {
                NonlinearSolverMethod::DampedNewton(DampedNewtonMethod::default())
            }
        };

        let mut contract = contract.clone();
        let mut initial_log_moles = initial_log_moles.to_vec();
        if let Some(bounds) = problem.log_mole_bounds() {
            if bounds.len() != initial_log_moles.len() {
                return Err(ReactionExtentError::DimensionMismatch(format!(
                    "RST log-mole bounds have {} entries for {} initial coordinates",
                    bounds.len(),
                    initial_log_moles.len()
                )));
            }
            clamp_initial_log_moles_to_roundoff_bounds(&mut initial_log_moles, bounds);
            contract = contract.with_log_mole_bounds(bounds)?;
        }
        let telemetry_mode = contract.telemetry_mode;
        let bound_problem = problem.bound_problem()?;
        #[cfg(test)]
        let options = {
            let mut options = contract.to_options()?;
            if rst_test_history_capture_enabled() {
                options.diagnostics.collect_history = true;
            }
            options
        };
        #[cfg(not(test))]
        let options = contract.to_options()?;
        method
            .solve(
                &bound_problem,
                DVector::from_vec(initial_log_moles),
                options,
            )
            .map(|result| {
                #[cfg(test)]
                if rst_test_history_capture_enabled() {
                    RST_TEST_CAPTURED_RESULTS.with(|captured| {
                        captured.borrow_mut().push(result.clone());
                    });
                }
                let statistics = &result.statistics;
                let counters_collected = statistics.availability.is_collected();
                let timings_collected = statistics.timings_collected;
                let callback_and_linear = statistics
                    .residual_duration
                    .saturating_add(statistics.jacobian_duration)
                    .saturating_add(statistics.linear_solve_duration);
                let iteration_telemetry =
                    (telemetry_mode == SolverTelemetryMode::Detailed).then(|| {
                        statistics
                            .attempts
                            .iter()
                            .map(|attempt| SolverIterationTelemetry {
                                iteration: attempt.iteration,
                                residual_evaluations: attempt.residual_evaluations,
                                state_residual_evaluations: attempt.state_residual_evaluations,
                                trial_residual_evaluations: attempt.trial_residual_evaluations,
                                jacobian_evaluations: attempt.jacobian_evaluations,
                                state_jacobian_evaluations: attempt.state_jacobian_evaluations,
                                trial_jacobian_evaluations: attempt.trial_jacobian_evaluations,
                                linear_factorizations: attempt.linear_factorizations,
                                linear_solves: attempt.linear_solves,
                                accepted_steps: attempt.accepted_steps,
                                rejected_steps: attempt.rejected_steps,
                                residual_micros: attempt.residual_duration.as_micros(),
                                jacobian_micros: attempt.jacobian_duration.as_micros(),
                                linear_solve_micros: attempt.linear_solve_duration.as_micros(),
                                linear_factorization_micros: attempt
                                    .linear_factorization_duration
                                    .as_micros(),
                                linear_system_solve_micros: attempt
                                    .linear_system_solve_duration
                                    .as_micros(),
                            })
                            .collect()
                    });
                RustedSciTheSolveOutcome {
                    solution: result.x.as_slice().to_vec(),
                    metrics: SolverAttemptMetrics {
                        termination: map_termination(result.termination.clone()),
                        backend_converged: matches!(
                            result.termination,
                            TerminationReason::Converged
                        ),
                        iterations: result.iterations,
                        residual_evaluations: statistics.residual_evaluations,
                        jacobian_evaluations: statistics.jacobian_evaluations,
                        linear_solves: statistics.linear_solves,
                        counters_collected,
                        timings_collected,
                        elapsed_millis: if timings_collected {
                            statistics.total_duration.as_millis()
                        } else {
                            0
                        },
                        elapsed_micros: if timings_collected {
                            statistics.total_duration.as_micros()
                        } else {
                            0
                        },
                        evaluation_timing: timings_collected.then(|| SolverEvaluationTiming {
                            residual_evaluation_micros: statistics.residual_duration.as_micros(),
                            jacobian_evaluation_micros: statistics.jacobian_duration.as_micros(),
                            legacy_same_iterate_callback_pairs: None,
                            linear_system_micros: Some(
                                statistics.linear_solve_duration.as_micros(),
                            ),
                            solver_overhead_micros: statistics
                                .total_duration
                                .saturating_sub(callback_and_linear)
                                .as_micros(),
                            matrix_operations_micros: 0,
                            norm_micros: 0,
                            trial_vector_micros: 0,
                            feasibility_micros: 0,
                            globalization_micros: 0,
                            state_residual_evaluations: counters_collected
                                .then_some(statistics.state_residual_evaluations),
                            trial_residual_evaluations: counters_collected
                                .then_some(statistics.trial_residual_evaluations),
                            state_jacobian_evaluations: counters_collected
                                .then_some(statistics.state_jacobian_evaluations),
                            trial_jacobian_evaluations: counters_collected
                                .then_some(statistics.trial_jacobian_evaluations),
                            linear_factorization_micros: Some(
                                statistics.linear_factorization_duration.as_micros(),
                            ),
                            linear_system_solve_micros: Some(
                                statistics.linear_system_solve_duration.as_micros(),
                            ),
                            linear_factorizations: counters_collected
                                .then_some(statistics.linear_factorizations),
                            accepted_steps: counters_collected.then_some(statistics.accepted_steps),
                            rejected_steps: counters_collected.then_some(statistics.rejected_steps),
                            reusable_trial_points: counters_collected
                                .then_some(statistics.reusable_trial_points),
                            iteration_telemetry,
                        }),
                    },
                }
            })
            .map_err(|error| self.map_solve_error(error))
    }

    /// Test-only access to the native RST result with iteration history.
    ///
    /// Production adapters deliberately disable history allocation. The
    /// parity fixture needs the same numerical call with RST's opt-in history
    /// enabled, while retaining the solver's own iteration semantics.
    #[cfg(test)]
    pub(crate) fn begin_history_capture_for_test() {
        RST_TEST_CAPTURED_RESULTS.with(|captured| captured.borrow_mut().clear());
        RST_TEST_HISTORY_CAPTURE_ENABLED.with(|enabled| *enabled.borrow_mut() = true);
    }

    /// Test-only drain of native results captured while the typed facade ran.
    #[cfg(test)]
    pub(crate) fn take_history_capture_for_test() -> Vec<SolveResult> {
        RST_TEST_HISTORY_CAPTURE_ENABLED.with(|enabled| *enabled.borrow_mut() = false);
        RST_TEST_CAPTURED_RESULTS.with(|captured| std::mem::take(&mut *captured.borrow_mut()))
    }

    /// Test-only access to the native RST result with iteration history.
    ///
    /// Production adapters deliberately disable history allocation. The
    /// parity fixture needs the same numerical call with RST's opt-in history
    /// enabled, while retaining the solver's own iteration semantics.
    #[cfg(test)]
    pub(crate) fn solve_with_history_for_test(
        self,
        problem: &RstPreparedProblem,
        initial_log_moles: &[f64],
        contract: &RustedSciTheSolveContract,
    ) -> Result<SolveResult, ReactionExtentError> {
        let method = match self {
            Self::LevenbergMarquardt => {
                NonlinearSolverMethod::LevenbergMarquardt(LevenbergMarquardtMethod::default())
            }
            Self::BacktrackingLevenbergMarquardt => {
                NonlinearSolverMethod::BacktrackingLevenbergMarquardt(
                    BacktrackingLevenbergMarquardtMethod::default(),
                )
            }
            Self::MinpackLevenbergMarquardt => NonlinearSolverMethod::LevenbergMarquardtMinpack(
                LevenbergMarquardtMinpack::default(),
            ),
            Self::NielsenLevenbergMarquardt => NonlinearSolverMethod::NielsenLevenbergMarquardt(
                NielsenLevenbergMarquardtMethod::default(),
            ),
            Self::TrustRegionLevenbergMarquardt => {
                NonlinearSolverMethod::TrustRegionLM(TrustRegionLMMethod::default())
            }
            Self::PowellDogleg => {
                NonlinearSolverMethod::PowellDogleg(PowellDoglegMethod::default())
            }
            Self::DampedNewton => {
                NonlinearSolverMethod::DampedNewton(DampedNewtonMethod::default())
            }
        };

        let mut contract = contract.clone();
        let mut initial_log_moles = initial_log_moles.to_vec();
        if let Some(bounds) = problem.log_mole_bounds() {
            if bounds.len() != initial_log_moles.len() {
                return Err(ReactionExtentError::DimensionMismatch(format!(
                    "RST log-mole bounds have {} entries for {} initial coordinates",
                    bounds.len(),
                    initial_log_moles.len()
                )));
            }
            clamp_initial_log_moles_to_roundoff_bounds(&mut initial_log_moles, bounds);
            contract = contract.with_log_mole_bounds(bounds)?;
        }

        let bound_problem = problem.bound_problem()?;
        let mut options = contract.to_options()?;
        options.diagnostics.collect_history = true;
        method
            .solve(
                &bound_problem,
                DVector::from_vec(initial_log_moles),
                options,
            )
            .map_err(|error| self.map_solve_error(error))
    }

    /// Replays Trust-Region LM while retaining every callback input, including
    /// the trial coordinate whose callback returns an error. Test-only wrapper
    /// preserves the provider's into-path when available.
    #[cfg(test)]
    pub(crate) fn solve_with_failure_trace_for_test(
        self,
        problem: &RstPreparedProblem,
        initial_log_moles: &[f64],
        contract: &RustedSciTheSolveContract,
    ) -> Result<SolveResult, RstSolveFailureTrace> {
        let mut initial_log_moles = initial_log_moles.to_vec();
        let mut contract = contract.clone();
        if let Some(bounds) = problem.log_mole_bounds() {
            if bounds.len() != initial_log_moles.len() {
                return Err(RstSolveFailureTrace {
                    error: ReactionExtentError::DimensionMismatch(format!(
                        "RST log-mole bounds have {} entries for {} initial coordinates",
                        bounds.len(),
                        initial_log_moles.len()
                    )),
                    residual_inputs: Vec::new(),
                    jacobian_inputs: Vec::new(),
                });
            }
            clamp_initial_log_moles_to_roundoff_bounds(&mut initial_log_moles, bounds);
            contract =
                contract
                    .with_log_mole_bounds(bounds)
                    .map_err(|error| RstSolveFailureTrace {
                        error,
                        residual_inputs: Vec::new(),
                        jacobian_inputs: Vec::new(),
                    })?;
        }
        let mut options = contract
            .to_options()
            .map_err(|error| RstSolveFailureTrace {
                error,
                residual_inputs: Vec::new(),
                jacobian_inputs: Vec::new(),
            })?;
        options.diagnostics.collect_history = true;
        let bound_problem = problem
            .bound_problem()
            .map_err(|error| RstSolveFailureTrace {
                error,
                residual_inputs: Vec::new(),
                jacobian_inputs: Vec::new(),
            })?;
        let observed = CallbackInputRecorder::new(bound_problem);
        let method = NonlinearSolverMethod::TrustRegionLM(TrustRegionLMMethod::default());
        match method.solve(&observed, DVector::from_vec(initial_log_moles), options) {
            Ok(result) => Ok(result),
            Err(error) => Err(RstSolveFailureTrace {
                error: self.map_solve_error(error),
                residual_inputs: observed.residual_inputs.into_inner(),
                jacobian_inputs: observed.jacobian_inputs.into_inner(),
            }),
        }
    }

    /// Preserves the distinction between a retryable numerical termination
    /// and a configuration or contract error. The policy layer may try the
    /// next backend only for the former category.
    fn map_solve_error(self, error: RstSolveError) -> ReactionExtentError {
        match error {
            RstSolveError::InvalidConfig(message) => ReactionExtentError::InvalidProblem {
                field: "rst_solve_options",
                message: format!("RustedSciThe {}: {message}", self.name()),
            },
            RstSolveError::DimensionMismatch {
                expected,
                actual,
                context,
            } => ReactionExtentError::DimensionMismatch(format!(
                "RustedSciThe {} {context}: expected {expected}, got {actual}",
                self.name()
            )),
            RstSolveError::InfeasibleInitialGuess {
                index,
                value,
                lower,
                upper,
            } => ReactionExtentError::InvalidProblem {
                field: "initial_log_moles",
                message: format!(
                    "RustedSciThe {}: coordinate {index} = {value} violates bounds [{lower}, {upper}]",
                    self.name()
                ),
            },
            RstSolveError::NonFiniteInitialGuess { index, value } => {
                ReactionExtentError::InvalidProblem {
                    field: "initial_log_moles",
                    message: format!(
                        "RustedSciThe {}: coordinate {index} is non-finite ({value})",
                        self.name()
                    ),
                }
            }
            RstSolveError::InvalidParameterSchema(message)
            | RstSolveError::ParameterSchemaMismatch(message)
            | RstSolveError::InvalidVariableSchema(message) => {
                ReactionExtentError::InvalidProblem {
                    field: "rst_symbolic_schema",
                    message: format!("RustedSciThe {}: {message}", self.name()),
                }
            }
            RstSolveError::UndeclaredSymbol {
                equation_index,
                name,
            } => ReactionExtentError::InvalidProblem {
                field: "rst_symbolic_schema",
                message: format!(
                    "RustedSciThe {}: equation {equation_index} references undeclared symbol {name}",
                    self.name()
                ),
            },
            RstSolveError::NonFiniteParameterValue { index, value } => {
                ReactionExtentError::InvalidProblem {
                    field: "rst_symbolic_parameter_values",
                    message: format!(
                        "RustedSciThe {}: parameter {index} is non-finite ({value})",
                        self.name()
                    ),
                }
            }
            RstSolveError::CompiledAotRuntimeUnavailable(message)
            | RstSolveError::CompiledAotArtifactMissing(message)
            | RstSolveError::CompiledAotArtifactNotBuilt(message)
            | RstSolveError::AotBuildFailed(message) => ReactionExtentError::InvalidProblem {
                field: "rst_backend",
                message: format!("RustedSciThe {}: {message}", self.name()),
            },
            RstSolveError::AotBuildOutputDirMissing => ReactionExtentError::InvalidProblem {
                field: "rst_backend",
                message: format!(
                    "RustedSciThe {} requested an AOT build without an output directory",
                    self.name()
                ),
            },
            RstSolveError::ResidualEvaluation(message) => ReactionExtentError::ResidualEvaluation(
                format!("RustedSciThe {}: {message}", self.name()),
            ),
            RstSolveError::JacobianEvaluation(message) => ReactionExtentError::JacobianEvaluation(
                format!("RustedSciThe {}: {message}", self.name()),
            ),
            RstSolveError::LinearSolveFailure(message) => ReactionExtentError::BackendFailure {
                backend: self.name().to_string(),
                kind: BackendFailureKind::LinearSolve,
                message,
            },
            RstSolveError::SingularJacobian => ReactionExtentError::BackendFailure {
                backend: self.name().to_string(),
                kind: BackendFailureKind::SingularJacobian,
                message: "RustedSciThe reported a singular or ill-conditioned Jacobian".to_string(),
            },
            RstSolveError::NumericalBreakdown(message) => ReactionExtentError::BackendFailure {
                backend: self.name().to_string(),
                kind: BackendFailureKind::NumericalBreakdown,
                message,
            },
        }
    }
}

#[cfg(test)]
struct CallbackInputRecorder<P> {
    inner: P,
    residual_inputs: RefCell<Vec<Vec<f64>>>,
    jacobian_inputs: RefCell<Vec<Vec<f64>>>,
}

#[cfg(test)]
impl<P> CallbackInputRecorder<P> {
    fn new(inner: P) -> Self {
        Self {
            inner,
            residual_inputs: RefCell::new(Vec::new()),
            jacobian_inputs: RefCell::new(Vec::new()),
        }
    }
}

#[cfg(test)]
impl<P: NonlinearProblem> NonlinearProblem for CallbackInputRecorder<P> {
    fn dimension(&self) -> usize {
        self.inner.dimension()
    }

    fn residual(&self, x: &DVector<f64>) -> Result<DVector<f64>, RstSolveError> {
        self.residual_inputs
            .borrow_mut()
            .push(x.as_slice().to_vec());
        self.inner.residual(x)
    }

    fn supports_residual_into(&self) -> bool {
        self.inner.supports_residual_into()
    }

    fn residual_into(&self, x: &DVector<f64>, out: &mut DVector<f64>) -> Result<(), RstSolveError> {
        self.residual_inputs
            .borrow_mut()
            .push(x.as_slice().to_vec());
        self.inner.residual_into(x, out)
    }
}

#[cfg(test)]
impl<P: JacobianProvider> JacobianProvider for CallbackInputRecorder<P> {
    fn jacobian(&self, x: &DVector<f64>) -> Result<DMatrix<f64>, RstSolveError> {
        self.jacobian_inputs
            .borrow_mut()
            .push(x.as_slice().to_vec());
        self.inner.jacobian(x)
    }

    fn supports_jacobian_into(&self) -> bool {
        self.inner.supports_jacobian_into()
    }

    fn jacobian_into(&self, x: &DVector<f64>, out: &mut DMatrix<f64>) -> Result<(), RstSolveError> {
        self.jacobian_inputs
            .borrow_mut()
            .push(x.as_slice().to_vec());
        self.inner.jacobian_into(x, out)
    }
}

#[cfg(test)]
mod error_mapping_tests {
    use super::{RstSolveError, RustedSciTheSolver};
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_nonlinear::{
        BackendFailureKind, ReactionExtentError,
    };

    #[test]
    fn singular_rst_jacobian_preserves_backend_failure_kind() {
        let error =
            RustedSciTheSolver::LevenbergMarquardt.map_solve_error(RstSolveError::SingularJacobian);

        assert!(matches!(
            error,
            ReactionExtentError::BackendFailure {
                kind: BackendFailureKind::SingularJacobian,
                ..
            }
        ));
    }
}

#[cfg(test)]
mod policy_tests {
    use super::RustedSciTheSolver;

    #[test]
    fn production_cascade_keeps_noisy_powell_dogleg_explicit_only() {
        let recommended = RustedSciTheSolver::recommended_cascade();

        assert!(recommended.contains(&RustedSciTheSolver::LevenbergMarquardt));
        assert!(recommended.contains(&RustedSciTheSolver::BacktrackingLevenbergMarquardt));
        assert!(recommended.contains(&RustedSciTheSolver::DampedNewton));
        assert!(!recommended.contains(&RustedSciTheSolver::PowellDogleg));
    }
}

#[cfg(test)]
mod solve_contract_tests {
    use super::{RustedSciTheSolveContract, clamp_initial_log_moles_to_roundoff_bounds};
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_nonlinear::ReactionExtentError;
    use RustedSciThe::numerical::SolverTelemetryMode;

    #[test]
    fn solve_contract_validates_numeric_budget_before_reaching_rst() {
        let error = RustedSciTheSolveContract::new(0.0, 10).unwrap_err();
        assert!(matches!(
            error,
            ReactionExtentError::InvalidProblem {
                field: "rst_solve_contract.tolerance",
                ..
            }
        ));

        let error = RustedSciTheSolveContract::new(1e-8, 0).unwrap_err();
        assert!(matches!(
            error,
            ReactionExtentError::InvalidProblem {
                field: "rst_solve_contract.max_iterations",
                ..
            }
        ));
    }

    #[test]
    fn solve_contract_maps_directly_into_rst_options() {
        let contract = RustedSciTheSolveContract::new(1e-8, 17).unwrap();
        let options = contract.to_options().unwrap();

        assert_eq!(options.tolerance, 1e-8);
        assert_eq!(options.max_iterations, 17);
        assert_eq!(options.diagnostics.telemetry_mode, SolverTelemetryMode::Off);
        assert!(!options.diagnostics.collect_history);
    }

    #[test]
    fn solve_contract_rejects_non_finite_or_reversed_log_mole_bounds() {
        let error = RustedSciTheSolveContract::new(1e-8, 17)
            .unwrap()
            .with_log_mole_bounds(&[(f64::NAN, 1.0)])
            .unwrap_err();
        assert!(matches!(
            error,
            ReactionExtentError::InvalidProblem {
                field: "rst_log_mole_bounds",
                ..
            }
        ));

        let error = RustedSciTheSolveContract::new(1e-8, 17)
            .unwrap()
            .with_log_mole_bounds(&[(2.0, 1.0)])
            .unwrap_err();
        assert!(matches!(
            error,
            ReactionExtentError::InvalidProblem {
                field: "rst_log_mole_bounds",
                ..
            }
        ));
    }

    #[test]
    fn initial_log_moles_clamp_only_removes_roundoff_at_bounds() {
        let lower = -4.0;
        let upper = 4.0;
        let mut near = vec![lower - 4.0 * f64::EPSILON, upper + 4.0 * f64::EPSILON];
        clamp_initial_log_moles_to_roundoff_bounds(&mut near, &[(lower, upper), (lower, upper)]);
        assert_eq!(near, vec![lower, upper]);

        let mut materially_outside = vec![lower - 1.0e-6, upper + 1.0e-6];
        clamp_initial_log_moles_to_roundoff_bounds(
            &mut materially_outside,
            &[(lower, upper), (lower, upper)],
        );
        assert_eq!(materially_outside, vec![lower - 1.0e-6, upper + 1.0e-6]);
    }
}

/// Removes only floating-point noise at a finite log-mole box boundary.
///
/// Recovery can translate a log seed by `ln(s)` and then reconstruct it by
/// `exp`/`ln`; a mathematically boundary-valued coordinate may consequently
/// land a few ULPs outside the capacity bound. The backend must still reject
/// a materially infeasible seed, so this is intentionally a roundoff-sized
/// correction rather than a general projection.
fn clamp_initial_log_moles_to_roundoff_bounds(
    initial_log_moles: &mut [f64],
    bounds: &[(f64, f64)],
) {
    for (value, &(lower, upper)) in initial_log_moles.iter_mut().zip(bounds) {
        let scale = value.abs().max(lower.abs()).max(upper.abs()).max(1.0);
        let tolerance = 8.0 * f64::EPSILON * scale;
        if *value < lower && lower - *value <= tolerance {
            *value = lower;
        } else if *value > upper && *value - upper <= tolerance {
            *value = upper;
        }
    }
}

/// Converts RST's engine-level termination enum into the equilibrium report
/// without stringifying away its machine-readable meaning.
fn map_termination(termination: TerminationReason) -> SolverTermination {
    match termination {
        TerminationReason::Converged => SolverTermination::Converged,
        TerminationReason::MaxIterations => SolverTermination::MaxIterations,
        TerminationReason::StepTooSmall => SolverTermination::StepTooSmall,
        TerminationReason::Stagnation => SolverTermination::Stagnation,
        TerminationReason::RejectedStepLimit => SolverTermination::RejectedStepLimit,
    }
}

/// Builds an RST-owned symbolic residual/Jacobian provider for the current
/// mutable workflow state.
///
/// This is the only bridge allowed to prepare a symbolic nonlinear problem.
/// It prevents the solve policy from mixing manual Jacobians with RST methods.
pub(crate) fn prepare_rst_symbolic_problem(
    solver: &mut EquilibriumLogMoles,
) -> Result<RstPreparedProblem, ReactionExtentError> {
    let thermochemistry = RstSymbolicThermochemistry::from_solver(solver)?;

    let equations = multiphase_equilibrium_residual_generator_sym(
        solver.stoich_matrix.clone(),
        solver.elem_composition.clone(),
        solver.elements_vector.clone(),
        thermochemistry.standard_gibbs,
        solver.phases.clone(),
        solver.P,
        solver.p0,
    )?;
    let variables = Expr::IndexedVars(solver.stoich_matrix.nrows(), "y")
        .0
        .into_iter()
        .map(|expression| expression.to_string())
        .collect();
    RstSymbolicProblem::from_expressions(
        equations,
        variables,
        Some(vec!["T".to_string()]),
        Some(DVector::from_vec(vec![solver.T])),
    )
    .map(RstPreparedProblem::new)
}

#[cfg(test)]
pub(crate) fn prepare_rst_symbolic_problem_from_prepared(
    prepared: &PreparedEquilibriumProblem,
) -> Result<RstPreparedProblem, ReactionExtentError> {
    prepare_rst_symbolic_problem_from_prepared_with_timing(
        prepared,
        EquilibriumTimingMode::Disabled,
    )
    .map(|(problem, _)| problem)
}

#[cfg(test)]
pub(crate) fn prepare_rst_symbolic_problem_from_prepared_with_timing(
    prepared: &PreparedEquilibriumProblem,
    timing_mode: EquilibriumTimingMode,
) -> Result<(RstPreparedProblem, RstPreparationTimingReport), ReactionExtentError> {
    prepare_rst_symbolic_problem_from_prepared_with_options(
        prepared,
        timing_mode,
        RstSymbolicFrontend::ExprLegacy,
    )
}

pub(crate) fn prepare_rst_symbolic_problem_from_prepared_with_options(
    prepared: &PreparedEquilibriumProblem,
    timing_mode: EquilibriumTimingMode,
    frontend: RstSymbolicFrontend,
) -> Result<(RstPreparedProblem, RstPreparationTimingReport), ReactionExtentError> {
    prepare_rst_symbolic_problem_from_prepared_with_frontend(
        prepared,
        timing_mode,
        resolve_requested_frontend(frontend),
    )
}

fn resolve_requested_frontend(frontend: RstSymbolicFrontend) -> SymbolicLambdifyFrontend {
    #[cfg(any(test, feature = "equilibrium-benchmark-rst-frontends"))]
    if let Some(frontend) = RST_TEST_FRONTEND_OVERRIDE.with(|current| *current.borrow()) {
        return frontend;
    }
    frontend.into()
}

pub(crate) fn prepare_rst_symbolic_problem_from_prepared_with_frontend(
    prepared: &PreparedEquilibriumProblem,
    timing_mode: EquilibriumTimingMode,
    frontend: SymbolicLambdifyFrontend,
) -> Result<(RstPreparedProblem, RstPreparationTimingReport), ReactionExtentError> {
    prepare_rst_symbolic_problem_from_prepared_with_timing_and_simplification(
        prepared,
        timing_mode,
        fixed_pt_symbolic_simplification_enabled(),
        frontend,
    )
}

fn fixed_pt_symbolic_simplification_enabled() -> bool {
    // Raw expressions are the production default. Keep the old simplified
    // route available for matched performance and regression measurements.
    if matches!(
        std::env::var("KITHE_RST_SYMBOLIC_SIMPLIFICATION").as_deref(),
        Ok("simplified") | Ok("on")
    ) {
        return true;
    }

    if matches!(
        std::env::var("KITHE_RST_SYMBOLIC_SIMPLIFICATION").as_deref(),
        Ok("raw") | Ok("off")
    ) {
        return false;
    }

    cfg!(feature = "equilibrium-benchmark-simplified-symbolic")
}

fn prepare_rst_symbolic_problem_from_prepared_with_timing_and_simplification(
    prepared: &PreparedEquilibriumProblem,
    timing_mode: EquilibriumTimingMode,
    simplify_residuals: bool,
    frontend: SymbolicLambdifyFrontend,
) -> Result<(RstPreparedProblem, RstPreparationTimingReport), ReactionExtentError> {
    let started = Instant::now();
    let mut timing = RstPreparationTimingCollector::new(timing_mode);
    let species_count = prepared.problem().species().len();
    let parameterized_gibbs = (0..species_count)
        .map(|index| Expr::Var(format!("g0_{index}")))
        .collect::<Vec<_>>();
    let equation_started = Instant::now();
    let equations = multiphase_equilibrium_residual_generator_sym_at_temperature_with_timing(
        prepared.reaction_basis().reactions.clone(),
        prepared.problem().element_composition().clone(),
        prepared.element_totals().to_vec(),
        parameterized_gibbs,
        prepared.problem().phases().to_vec(),
        prepared.problem().conditions().pressure(),
        prepared.problem().conditions().reference_pressure(),
        Expr::Var("T".to_string()),
        simplify_residuals,
        Some(&mut timing),
    )?;
    timing.record(
        RstPreparationTimingStage::EquationGeneration,
        equation_started.elapsed(),
    );
    let variables = timing.measure(RstPreparationTimingStage::VariableMaterialization, || {
        Expr::IndexedVars(species_count, "y")
            .0
            .into_iter()
            .map(|expression| expression.to_string())
            .collect::<Vec<_>>()
    });
    let standard_gibbs_values = timing
        .measure(RstPreparationTimingStage::ParameterEvaluation, || {
            standard_gibbs_at_conditions(prepared)
        })?;
    let mut parameter_names = Vec::with_capacity(species_count + 1);
    parameter_names.push("T".to_string());
    parameter_names.extend((0..species_count).map(|index| format!("g0_{index}")));
    let mut parameter_values = Vec::with_capacity(species_count + 1);
    parameter_values.push(prepared.problem().conditions().temperature());
    parameter_values.extend_from_slice(&standard_gibbs_values);
    let log_mole_bounds = timing.measure(RstPreparationTimingStage::LogMoleBounds, || {
        prepared.finite_log_mole_bounds()
    })?;
    let (problem, nested_timing) = RstSymbolicProblem::from_expressions_with_timing_and_frontend(
        equations,
        variables,
        Some(parameter_names),
        Some(DVector::from_vec(parameter_values)),
        timing_mode,
        frontend,
    )?;
    timing.accumulate(nested_timing);
    timing.set_total(started.elapsed());
    Ok((
        RstPreparedProblem::with_fixed_pt_parameters(problem, log_mole_bounds, species_count),
        timing.finish(),
    ))
}

/// Recreates the pre-parameterization fixed-`P,T` graph for a regression
/// comparison.  Production code must use
/// [`prepare_rst_symbolic_problem_from_prepared`]: rebuilding a baked graph at
/// every coefficient boundary is precisely the cost this adapter removes.
///
/// Keeping this helper test-only lets a real-data test distinguish a solver's
/// inherent convergence limitation from an accidental change introduced by
/// the reusable `g0_i` parameterization.
#[cfg(test)]
pub(crate) fn prepare_baked_rst_symbolic_problem_for_test(
    prepared: &PreparedEquilibriumProblem,
    standard_gibbs: &[Expr],
) -> Result<RstPreparedProblem, ReactionExtentError> {
    if standard_gibbs.len() != prepared.problem().species().len() {
        return Err(ReactionExtentError::DimensionMismatch(format!(
            "symbolic Gibbs snapshot has {} entries for {} species",
            standard_gibbs.len(),
            prepared.problem().species().len()
        )));
    }
    let species_count = prepared.problem().species().len();
    let equations = multiphase_equilibrium_residual_generator_sym(
        prepared.reaction_basis().reactions.clone(),
        prepared.problem().element_composition().clone(),
        prepared.element_totals().to_vec(),
        standard_gibbs.to_vec(),
        prepared.problem().phases().to_vec(),
        prepared.problem().conditions().pressure(),
        prepared.problem().conditions().reference_pressure(),
    )?;
    let variables = Expr::IndexedVars(species_count, "y")
        .0
        .into_iter()
        .map(|expression| expression.to_string())
        .collect::<Vec<_>>();
    let log_mole_bounds = prepared.finite_log_mole_bounds()?;
    RstSymbolicProblem::from_expressions(
        equations,
        variables,
        Some(vec!["T".to_string()]),
        Some(DVector::from_vec(vec![
            prepared.problem().conditions().temperature(),
        ])),
    )
    .map(|problem| RstPreparedProblem {
        problem,
        log_mole_bounds: Some(log_mole_bounds),
        fixed_pt_species_count: None,
    })
}

/// Evaluates standard-state Gibbs closures at the prepared problem's fixed
/// temperature. The RST graph owns only plain finite parameters, never the
/// domain closures themselves.
fn standard_gibbs_at_conditions(
    prepared: &PreparedEquilibriumProblem,
) -> Result<Vec<f64>, ReactionExtentError> {
    let temperature = prepared.problem().conditions().temperature();
    if !temperature.is_finite() {
        return Err(ReactionExtentError::InvalidProblem {
            field: "standard_gibbs_temperature",
            message: "prepared fixed-P,T temperature must be finite".to_string(),
        });
    }
    prepared
        .problem()
        .gibbs()
        .iter()
        .enumerate()
        .map(|(index, function)| {
            let value = function(temperature);
            if value.is_finite() {
                Ok(value)
            } else {
                Err(ReactionExtentError::InvalidProblem {
                    field: "standard_gibbs",
                    message: format!(
                        "standard Gibbs closure {index} returned non-finite value at {temperature} K"
                    ),
                })
            }
        })
        .collect()
}

/// Builds the coupled P,H RST graph and optionally reports its internal
/// preparation split. The symbolic expression and backend materialization
/// boundaries are diagnostic only; they do not alter the graph or solver
/// contract.
pub(crate) fn prepare_rst_symbolic_ph_problem_with_timing(
    formulation: &PreparedPhFormulation,
    timing_mode: EquilibriumTimingMode,
    frontend: RstSymbolicFrontend,
) -> Result<(RstPreparedProblem, EquilibriumTimingReport), ReactionExtentError> {
    let started = Instant::now();
    let mut timing = EquilibriumTimingCollector::new(timing_mode);
    let mut rst_timing = RstPreparationTimingCollector::new(timing_mode);
    let symbolic = timing
        .measure(EquilibriumTimingStage::ThermochemistryPreparation, || {
            formulation
                .thermochemistry()
                .symbolic_for_bounds(formulation.temperature_bounds())
        })
        .map_err(|reason| ReactionExtentError::UnsupportedBackendCapability {
            backend: "RustedSciThe".to_string(),
            capability: "monolithic_p_h_symbolic_single_interval",
            alternatives: format!(
                "legacy LM, NR, or trust-region backends; symbolic P,H requires one native coefficient interval per component ({reason})"
            ),
        })?
        .ok_or_else(|| ReactionExtentError::UnsupportedBackendCapability {
            backend: "RustedSciThe".to_string(),
            capability: "monolithic_p_h_symbolic_thermochemistry",
            alternatives: "legacy LM, NR, or trust-region backends; resolved SubsData thermochemistry".to_string(),
        })?;
    let prepared = formulation.prepared_pt();
    let species_count = prepared.problem().species().len();
    if symbolic.standard_gibbs().len() != species_count
        || symbolic.enthalpy().len() != species_count
    {
        return Err(ReactionExtentError::DimensionMismatch(format!(
            "symbolic P,H thermochemistry has {} Gibbs and {} enthalpy expressions for {species_count} species",
            symbolic.standard_gibbs().len(),
            symbolic.enthalpy().len(),
        )));
    }
    if formulation.pt_row_scale().len() != species_count {
        return Err(ReactionExtentError::DimensionMismatch(format!(
            "symbolic P,H row scale has {} entries for {species_count} fixed-P,T rows",
            formulation.pt_row_scale().len()
        )));
    }

    let equations = timing.measure(EquilibriumTimingStage::SymbolicConstruction, || {
        let temperature = bounded_temperature_expression(formulation.temperature_bounds());
        let (standard_gibbs, enthalpy) =
            rst_timing.measure(RstPreparationTimingStage::TemperatureSubstitution, || {
                let standard_gibbs = symbolic
                    .standard_gibbs()
                    .iter()
                    .cloned()
                    .map(|expression| substitute_temperature(expression, &temperature))
                    .collect::<Vec<_>>();
                let enthalpy = symbolic
                    .enthalpy()
                    .iter()
                    .cloned()
                    .map(|expression| substitute_temperature(expression, &temperature))
                    .collect::<Vec<_>>();
                (standard_gibbs, enthalpy)
            });
        let equation_started = Instant::now();
        let mut equations =
            multiphase_equilibrium_residual_generator_sym_at_temperature_with_timing(
                prepared.reaction_basis().reactions.clone(),
                prepared.problem().element_composition().clone(),
                prepared.element_totals().to_vec(),
                standard_gibbs,
                prepared.problem().phases().to_vec(),
                prepared.problem().conditions().pressure(),
                prepared.problem().conditions().reference_pressure(),
                temperature.clone(),
                false,
                Some(&mut rst_timing),
            )?;
        rst_timing.record(
            RstPreparationTimingStage::EquationGeneration,
            equation_started.elapsed(),
        );
        rst_timing.measure(RstPreparationTimingStage::RowScaling, || {
            for (equation, scale) in equations.iter_mut().zip(formulation.pt_row_scale()) {
                // Row scaling is multiplication by a validated constant. The
                // generator already simplifies each physical residual; running
                // the symbolic simplifier again here only rebuilds the same large
                // expression tree before RST differentiates it.
                *equation = equation.clone() / Expr::Const(*scale);
            }
        });

        rst_timing.measure(RstPreparationTimingStage::EnthalpyEquation, || {
            let y = Expr::IndexedVars(species_count, "y").0;
            let total_enthalpy = y
                .into_iter()
                .zip(enthalpy)
                .fold(Expr::Const(0.0), |total, (log_moles, molar_enthalpy)| {
                    total + log_moles.exp() * molar_enthalpy
                });
            equations.push(
                (total_enthalpy - Expr::Var("target_enthalpy".to_string()))
                    / Expr::Var("enthalpy_scale".to_string()),
            );
        });
        Ok::<_, ReactionExtentError>(equations)
    })?;

    let problem = timing.measure(EquilibriumTimingStage::RstPreparation, || {
        let mut variables =
            rst_timing.measure(RstPreparationTimingStage::VariableMaterialization, || {
                Expr::IndexedVars(species_count, "y")
                    .0
                    .into_iter()
                    .map(|expression| expression.to_string())
                    .collect::<Vec<_>>()
            });
        variables.push("theta_T".to_string());
        let (problem, nested_timing) =
            RstSymbolicProblem::from_expressions_with_timing_and_frontend(
                equations,
                variables,
                Some(vec![
                    "target_enthalpy".to_string(),
                    "enthalpy_scale".to_string(),
                ]),
                Some(DVector::from_vec(vec![
                    formulation.target_enthalpy(),
                    formulation.enthalpy_scale().joules(),
                ])),
                timing_mode,
                resolve_requested_frontend(frontend),
            )?;
        rst_timing.accumulate(nested_timing);
        Ok::<_, ReactionExtentError>(RstPreparedProblem::new(problem))
    })?;
    rst_timing.set_total(started.elapsed());
    timing.record_rst_preparation_timing(rst_timing.finish());
    timing.set_total(started.elapsed());
    Ok((problem, timing.finish()))
}

fn bounded_temperature_expression(
    bounds: crate::Thermodynamics::ChemEquilibrium::equilibrium_constraints::TemperatureBounds,
) -> Expr {
    let theta = Expr::Var("theta_T".to_string());
    let sigmoid = Expr::Const(1.0) / (Expr::Const(1.0) + (Expr::Const(-1.0) * theta).exp());
    Expr::Const(bounds.lower()) + Expr::Const(bounds.upper() - bounds.lower()) * sigmoid
}

/// Replaces the thermochemistry convention `T` with the monolithic bounded
/// temperature expression without string parsing or losing expression shape.
fn substitute_temperature(expression: Expr, temperature: &Expr) -> Expr {
    match expression {
        Expr::Var(name) if name == "T" => temperature.clone(),
        Expr::Var(_) | Expr::Const(_) => expression,
        Expr::Add(left, right) => {
            substitute_temperature(*left, temperature) + substitute_temperature(*right, temperature)
        }
        Expr::Sub(left, right) => {
            substitute_temperature(*left, temperature) - substitute_temperature(*right, temperature)
        }
        Expr::Mul(left, right) => {
            substitute_temperature(*left, temperature) * substitute_temperature(*right, temperature)
        }
        Expr::Div(left, right) => {
            substitute_temperature(*left, temperature) / substitute_temperature(*right, temperature)
        }
        Expr::Pow(left, right) => substitute_temperature(*left, temperature)
            .pow(substitute_temperature(*right, temperature)),
        Expr::Exp(value) => substitute_temperature(*value, temperature).exp(),
        Expr::Ln(value) => substitute_temperature(*value, temperature).ln(),
        Expr::sin(value) => Expr::sin(Box::new(substitute_temperature(*value, temperature))),
        Expr::cos(value) => Expr::cos(Box::new(substitute_temperature(*value, temperature))),
        Expr::tg(value) => Expr::tg(Box::new(substitute_temperature(*value, temperature))),
        Expr::ctg(value) => Expr::ctg(Box::new(substitute_temperature(*value, temperature))),
        Expr::arcsin(value) => Expr::arcsin(Box::new(substitute_temperature(*value, temperature))),
        Expr::arccos(value) => Expr::arccos(Box::new(substitute_temperature(*value, temperature))),
        Expr::arctg(value) => Expr::arctg(Box::new(substitute_temperature(*value, temperature))),
        Expr::arcctg(value) => Expr::arcctg(Box::new(substitute_temperature(*value, temperature))),
    }
}

/// Converts the historical preferred solver into the corresponding explicit
/// legacy backend. Kept here so the policy module has one canonical mapping.
pub fn legacy_backend(solver: Solvers) -> super::equilibrium_solver_policy::SolverBackend {
    super::equilibrium_solver_policy::SolverBackend::Legacy(solver)
}
