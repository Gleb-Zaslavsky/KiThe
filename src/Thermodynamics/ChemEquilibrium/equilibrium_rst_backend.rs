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
    SolverAttemptMetrics, SolverEvaluationTiming, SolverTermination,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_workflows::{
    multiphase_equilibrium_residual_generator_sym,
    multiphase_equilibrium_residual_generator_sym_at_temperature,
};
use crate::Thermodynamics::User_substances_error::SubsDataError;
use RustedSciThe::numerical::Nonlinear_systems::error::SolveError as RstSolveError;
use RustedSciThe::numerical::Nonlinear_systems::prelude::{
    DampedNewtonMethod, LevenbergMarquardtMethod, LevenbergMarquardtMinpack,
    NielsenLevenbergMarquardtMethod, NonlinearSolverMethod, PowellDoglegMethod, SolveOptions,
    SymbolicNonlinearProblem, SymbolicProblemOptions, TerminationReason, TrustRegionLMMethod,
};
use RustedSciThe::numerical::Nonlinear_systems::problem::{
    Bounds, JacobianProvider, NonlinearProblem,
};
use RustedSciThe::symbolic::symbolic_engine::Expr;
use nalgebra::{DMatrix, DVector};
use std::cell::Cell;
use std::time::{Duration, Instant};

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
            log_mole_bounds: None,
        })
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
            ..SolveOptions::default()
        })
    }
}

/// Opaque prepared symbolic problem owned by the RST adapter layer.
///
/// The domain layer may carry this value around as a prepared backend request,
/// but it should never inspect or construct RST internals directly.
pub(crate) struct RstPreparedProblem {
    problem: SymbolicNonlinearProblem,
    log_mole_bounds: Option<Vec<(f64, f64)>>,
    fixed_pt_species_count: Option<usize>,
}

impl RstPreparedProblem {
    pub(crate) fn new(problem: SymbolicNonlinearProblem) -> Self {
        Self {
            problem,
            log_mole_bounds: None,
            fixed_pt_species_count: None,
        }
    }

    fn with_fixed_pt_parameters(
        problem: SymbolicNonlinearProblem,
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

    fn as_problem(&self) -> &SymbolicNonlinearProblem {
        &self.problem
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
        self.problem
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
        self.problem
            .jacobian(&DVector::from_column_slice(log_moles))
            .map_err(|error| {
                ReactionExtentError::ResidualEvaluation(format!(
                    "test evaluation of parameterized RST Jacobian failed: {error}"
                ))
            })
    }
}

/// Borrowed symbolic provider with callback-level timing for release evidence.
///
/// RST owns the numerical loop and its symbolic Jacobian. This wrapper does
/// not reinterpret either; it only measures the trait calls that cross from
/// the solver into the prepared symbolic problem. The remaining backend time
/// therefore covers linear algebra, step control, and engine orchestration.
struct TimedSymbolicProblem<'a> {
    inner: &'a SymbolicNonlinearProblem,
    residual_elapsed: Cell<Duration>,
    jacobian_elapsed: Cell<Duration>,
}

impl<'a> TimedSymbolicProblem<'a> {
    fn new(inner: &'a SymbolicNonlinearProblem) -> Self {
        Self {
            inner,
            residual_elapsed: Cell::new(Duration::ZERO),
            jacobian_elapsed: Cell::new(Duration::ZERO),
        }
    }

    fn timing(&self, total: Duration) -> SolverEvaluationTiming {
        let residual = self.residual_elapsed.get();
        let jacobian = self.jacobian_elapsed.get();
        SolverEvaluationTiming {
            residual_evaluation_micros: residual.as_micros(),
            jacobian_evaluation_micros: jacobian.as_micros(),
            solver_overhead_micros: total
                .saturating_sub(residual)
                .saturating_sub(jacobian)
                .as_micros(),
        }
    }
}

impl NonlinearProblem for TimedSymbolicProblem<'_> {
    fn dimension(&self) -> usize {
        self.inner.dimension()
    }

    fn residual(&self, x: &DVector<f64>) -> Result<DVector<f64>, RstSolveError> {
        let started = Instant::now();
        let result = self.inner.residual(x);
        self.residual_elapsed
            .set(self.residual_elapsed.get() + started.elapsed());
        result
    }
}

impl JacobianProvider for TimedSymbolicProblem<'_> {
    fn jacobian(&self, x: &DVector<f64>) -> Result<DMatrix<f64>, RstSolveError> {
        let started = Instant::now();
        let result = self.inner.jacobian(x);
        self.jacobian_elapsed
            .set(self.jacobian_elapsed.get() + started.elapsed());
        result
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

/// RST nonlinear methods suitable for the log-moles equilibrium system.
///
/// The ordered production default starts with ordinary LM as a transparent
/// baseline, then adds guarded least-squares/trust-region methods before
/// falling back to a line-search Newton method. Individual methods remain
/// selectable even when they are deliberately omitted from that default.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RustedSciTheSolver {
    /// Classical Levenberg-Marquardt.
    LevenbergMarquardt,
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
        ]
    }

    /// Stable identifier used in reports and tests.
    pub fn name(self) -> &'static str {
        match self {
            Self::LevenbergMarquardt => "rst_levenberg_marquardt",
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

        let started = Instant::now();
        let mut contract = contract.clone();
        if let Some(bounds) = problem.log_mole_bounds() {
            if bounds.len() != initial_log_moles.len() {
                return Err(ReactionExtentError::DimensionMismatch(format!(
                    "RST log-mole bounds have {} entries for {} initial coordinates",
                    bounds.len(),
                    initial_log_moles.len()
                )));
            }
            contract = contract.with_log_mole_bounds(bounds)?;
        }
        let timed_problem = TimedSymbolicProblem::new(problem.as_problem());
        method
            .solve(
                &timed_problem,
                DVector::from_vec(initial_log_moles.to_vec()),
                contract.to_options()?,
            )
            .map(|result| {
                let elapsed = started.elapsed();
                RustedSciTheSolveOutcome {
                    solution: result.x.as_slice().to_vec(),
                    metrics: SolverAttemptMetrics {
                        termination: map_termination(result.termination.clone()),
                        backend_converged: matches!(
                            result.termination,
                            TerminationReason::Converged
                        ),
                        iterations: result.iterations,
                        residual_evaluations: result.statistics.residual_evaluations,
                        jacobian_evaluations: result.statistics.jacobian_evaluations,
                        linear_solves: result.statistics.linear_solves,
                        elapsed_millis: elapsed.as_millis(),
                        evaluation_timing: Some(timed_problem.timing(elapsed)),
                    },
                }
            })
            .map_err(|error| self.map_solve_error(error))
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
        assert!(recommended.contains(&RustedSciTheSolver::DampedNewton));
        assert!(!recommended.contains(&RustedSciTheSolver::PowellDogleg));
    }
}

#[cfg(test)]
mod solve_contract_tests {
    use super::RustedSciTheSolveContract;
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_nonlinear::ReactionExtentError;

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
    let options = SymbolicProblemOptions::new()
        .with_variables(variables)
        .with_equation_parameters(vec!["T".to_string()])
        .with_equation_parameter_values(DVector::from_vec(vec![solver.T]))
        .with_lambdify_backend();

    SymbolicNonlinearProblem::from_expressions_with_options(equations, options)
        .map(RstPreparedProblem::new)
        .map_err(|error| {
            ReactionExtentError::ResidualEvaluation(format!(
                "failed to prepare RustedSciThe symbolic equilibrium problem: {error}"
            ))
        })
}

/// Builds the same symbolic RST problem from immutable prepared data.
///
/// The old helper remains for the compatibility solver, but the canonical
/// fixed-`P,T` facade uses this function so symbolic preparation cannot reopen
/// or mutate `SubsData` through `EquilibriumLogMoles`.
pub(crate) fn prepare_rst_symbolic_problem_from_prepared(
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
    let parameterized_gibbs = (0..species_count)
        .map(|index| Expr::Var(format!("g0_{index}")))
        .collect::<Vec<_>>();
    let equations = multiphase_equilibrium_residual_generator_sym(
        prepared.reaction_basis().reactions.clone(),
        prepared.problem().element_composition().clone(),
        prepared.element_totals().to_vec(),
        parameterized_gibbs,
        prepared.problem().phases().to_vec(),
        prepared.problem().conditions().pressure(),
        prepared.problem().conditions().reference_pressure(),
    )?;
    let variables = Expr::IndexedVars(species_count, "y")
        .0
        .into_iter()
        .map(|expression| expression.to_string())
        .collect();
    let standard_gibbs_values = standard_gibbs_at_conditions(prepared)?;
    let mut parameter_names = Vec::with_capacity(species_count + 1);
    parameter_names.push("T".to_string());
    parameter_names.extend((0..species_count).map(|index| format!("g0_{index}")));
    let mut parameter_values = Vec::with_capacity(species_count + 1);
    parameter_values.push(prepared.problem().conditions().temperature());
    parameter_values.extend_from_slice(&standard_gibbs_values);
    let options = SymbolicProblemOptions::new()
        .with_variables(variables)
        .with_equation_parameters(parameter_names)
        .with_equation_parameter_values(DVector::from_vec(parameter_values))
        .with_lambdify_backend();
    let log_mole_bounds = prepared.finite_log_mole_bounds()?;
    SymbolicNonlinearProblem::from_expressions_with_options(equations, options)
        .map(|problem| {
            RstPreparedProblem::with_fixed_pt_parameters(problem, log_mole_bounds, species_count)
        })
        .map_err(|error| {
            ReactionExtentError::ResidualEvaluation(format!(
                "failed to prepare RustedSciThe symbolic equilibrium problem: {error}"
            ))
        })
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
        .collect();
    let options = SymbolicProblemOptions::new()
        .with_variables(variables)
        .with_equation_parameters(vec!["T".to_string()])
        .with_equation_parameter_values(DVector::from_vec(vec![
            prepared.problem().conditions().temperature(),
        ]))
        .with_lambdify_backend();
    let log_mole_bounds = prepared.finite_log_mole_bounds()?;
    SymbolicNonlinearProblem::from_expressions_with_options(equations, options)
        .map(|problem| RstPreparedProblem {
            problem,
            log_mole_bounds: Some(log_mole_bounds),
            fixed_pt_species_count: None,
        })
        .map_err(|error| {
            ReactionExtentError::ResidualEvaluation(format!(
                "failed to prepare baked RustedSciThe symbolic equilibrium problem: {error}"
            ))
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

/// Builds an RST-owned symbolic residual/Jacobian provider for one coupled
/// fixed-phase `P,H` formulation.
///
/// The temperature coordinate is a true unknown, not a mutable equation
/// parameter. The function therefore substitutes the same bounded
/// `T(theta)` transform used by the analytic formulation into both `G0(T)`
/// and `H(T)`, appends the scaled enthalpy row, and lets RST differentiate the
/// complete `(N + 1)` system. No fixed-`P,T` symbolic payload is reused here.
pub(crate) fn prepare_rst_symbolic_ph_problem(
    formulation: &PreparedPhFormulation,
) -> Result<RstPreparedProblem, ReactionExtentError> {
    let symbolic = formulation
        .thermochemistry()
        .symbolic_for_bounds(formulation.temperature_bounds())
        .map_err(|reason| ReactionExtentError::UnsupportedBackendCapability {
            backend: "RustedSciThe".to_string(),
            capability: "monolithic_p_h_symbolic_single_interval",
            alternatives: format!(
                "legacy LM, NR, or trust-region backends; symbolic P,H requires one native coefficient interval per component ({reason})"
            ),
        })?
        .ok_or_else(|| {
        ReactionExtentError::UnsupportedBackendCapability {
            backend: "RustedSciThe".to_string(),
            capability: "monolithic_p_h_symbolic_thermochemistry",
            alternatives: "legacy LM, NR, or trust-region backends; resolved SubsData thermochemistry".to_string(),
        }
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

    let temperature = bounded_temperature_expression(formulation.temperature_bounds());
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
    let mut equations = multiphase_equilibrium_residual_generator_sym_at_temperature(
        prepared.reaction_basis().reactions.clone(),
        prepared.problem().element_composition().clone(),
        prepared.element_totals().to_vec(),
        standard_gibbs,
        prepared.problem().phases().to_vec(),
        prepared.problem().conditions().pressure(),
        prepared.problem().conditions().reference_pressure(),
        temperature,
    )?;
    for (equation, scale) in equations.iter_mut().zip(formulation.pt_row_scale()) {
        *equation = (equation.clone() / Expr::Const(*scale)).simplify();
    }

    let y = Expr::IndexedVars(species_count, "y").0;
    let total_enthalpy = y
        .into_iter()
        .zip(enthalpy)
        .fold(Expr::Const(0.0), |total, (log_moles, molar_enthalpy)| {
            total + log_moles.exp() * molar_enthalpy
        });
    equations.push(
        ((total_enthalpy - Expr::Var("target_enthalpy".to_string()))
            / Expr::Var("enthalpy_scale".to_string()))
        .simplify(),
    );

    let mut variables = Expr::IndexedVars(species_count, "y")
        .0
        .into_iter()
        .map(|expression| expression.to_string())
        .collect::<Vec<_>>();
    variables.push("theta_T".to_string());
    SymbolicNonlinearProblem::from_expressions_with_options(
        equations,
        SymbolicProblemOptions::new()
            .with_variables(variables)
            .with_equation_parameters(vec![
                "target_enthalpy".to_string(),
                "enthalpy_scale".to_string(),
            ])
            .with_equation_parameter_values(DVector::from_vec(vec![
                formulation.target_enthalpy(),
                formulation.enthalpy_scale().joules(),
            ]))
            .with_lambdify_backend(),
    )
    .map(RstPreparedProblem::new)
    .map_err(|error| {
        ReactionExtentError::ResidualEvaluation(format!(
            "failed to prepare RustedSciThe symbolic P,H equilibrium problem: {error}"
        ))
    })
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
