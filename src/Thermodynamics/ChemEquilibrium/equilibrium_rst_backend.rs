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
use crate::Thermodynamics::ChemEquilibrium::equilibrium_solver_policy::{
    SolverAttemptMetrics, SolverTermination,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_workflows::multiphase_equilibrium_residual_generator_sym;
use crate::Thermodynamics::User_substances_error::SubsDataError;
use RustedSciThe::numerical::Nonlinear_systems::error::SolveError as RstSolveError;
use RustedSciThe::numerical::Nonlinear_systems::prelude::{
    DampedNewtonMethod, LevenbergMarquardtMethod, LevenbergMarquardtMinpack,
    NielsenLevenbergMarquardtMethod, NonlinearSolverMethod, PowellDoglegMethod, SolveOptions,
    SymbolicNonlinearProblem, SymbolicProblemOptions, TerminationReason, TrustRegionLMMethod,
};
use RustedSciThe::symbolic::symbolic_engine::Expr;
use nalgebra::DVector;
use std::time::Instant;

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
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RustedSciTheSolveContract {
    /// Convergence tolerance forwarded to RST.
    pub tolerance: f64,
    /// Maximum nonlinear iterations forwarded to RST.
    pub max_iterations: usize,
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
        })
    }

    /// Converts the typed contract into the concrete RST options struct.
    pub fn to_options(self) -> SolveOptions {
        SolveOptions {
            tolerance: self.tolerance,
            max_iterations: self.max_iterations,
            ..SolveOptions::default()
        }
    }
}

/// Opaque prepared symbolic problem owned by the RST adapter layer.
///
/// The domain layer may carry this value around as a prepared backend request,
/// but it should never inspect or construct RST internals directly.
pub(crate) struct RstPreparedProblem(SymbolicNonlinearProblem);

impl RstPreparedProblem {
    pub(crate) fn new(problem: SymbolicNonlinearProblem) -> Self {
        Self(problem)
    }

    fn as_problem(&self) -> &SymbolicNonlinearProblem {
        &self.0
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
/// The ordered default starts with ordinary LM as a transparent baseline,
/// then adds increasingly guarded least-squares/trust-region methods before
/// falling back to a line-search Newton method.
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
    PowellDogleg,
    /// Newton method protected by backtracking line search.
    DampedNewton,
}

impl RustedSciTheSolver {
    /// Default robust order for newly migrated equilibrium solves.
    pub fn recommended_cascade() -> Vec<Self> {
        vec![
            Self::LevenbergMarquardt,
            Self::MinpackLevenbergMarquardt,
            Self::NielsenLevenbergMarquardt,
            Self::TrustRegionLevenbergMarquardt,
            Self::PowellDogleg,
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
        contract: RustedSciTheSolveContract,
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
        method
            .solve(
                problem.as_problem(),
                DVector::from_vec(initial_log_moles.to_vec()),
                contract.to_options(),
            )
            .map(|result| RustedSciTheSolveOutcome {
                solution: result.x.as_slice().to_vec(),
                metrics: SolverAttemptMetrics {
                    termination: map_termination(result.termination.clone()),
                    backend_converged: matches!(result.termination, TerminationReason::Converged),
                    iterations: result.iterations,
                    residual_evaluations: result.statistics.residual_evaluations,
                    jacobian_evaluations: result.statistics.jacobian_evaluations,
                    linear_solves: result.statistics.linear_solves,
                    elapsed_millis: started.elapsed().as_millis(),
                },
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
        let options = contract.to_options();

        assert_eq!(options.tolerance, 1e-8);
        assert_eq!(options.max_iterations, 17);
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

/// Converts the historical preferred solver into the corresponding explicit
/// legacy backend. Kept here so the policy module has one canonical mapping.
pub fn legacy_backend(solver: Solvers) -> super::equilibrium_solver_policy::SolverBackend {
    super::equilibrium_solver_policy::SolverBackend::Legacy(solver)
}
