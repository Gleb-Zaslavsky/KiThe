//! Immutable numerical runner for one prepared fixed-`P,T` equilibrium problem.
//!
//! [`PreparedEquilibriumProblem`](super::equilibrium_problem::PreparedEquilibriumProblem)
//! owns the validated chemistry and all equation-building decisions. This
//! module owns only the numerical policy and backend dispatch. Keeping those
//! responsibilities separate prevents the public fixed-phase workflow from
//! constructing the historical mutable `EquilibriumLogMoles` orchestration
//! object merely to solve one already prepared problem.

use crate::Thermodynamics::ChemEquilibrium::equilibrium_activity::PhaseActivityModel;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_backend_adapter::EquilibriumNonlinearBackend;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_constant_cross_validation::{
    classify_equilibrium_constant_cross_validation, EquilibriumConstantCrossValidationStatus,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_constant_problem::EquilibriumConstantProblem;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_constant_solver::{
    EquilibriumConstantSolver, EquilibriumConstantSolverMode,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_constant_validation::{
    EquilibriumConstantValidationMode, EquilibriumConstantValidationTolerances,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_log_moles::{
    compute_species_moles, EquilibriumLogMoles, EquilibriumSolverSettings,
    LEGACY_MOLE_FEASIBILITY_TOLERANCE,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_nonlinear::ReactionExtentError;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_problem::{
    EquilibriumSolution, LogMolesInitialGuess, PreparedEquilibriumProblem,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_rst_backend::{
    prepare_rst_symbolic_problem_from_prepared, RstPreparedProblem,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_solver_policy::{
    EquilibriumSolveReport, MultiStartAttemptReport, MultiStartSolveReport, SolverBackend,
    SolverCascadeBudget, SolverPolicy,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_validation::{
    validate_equilibrium_candidate, EquilibriumAcceptanceCriteria, EquilibriumCandidateResiduals,
};
use nalgebra::DMatrix;
use std::time::{Duration, Instant};
use RustedSciThe::symbolic::symbolic_engine::Expr;

/// Accepted output of one immutable prepared-problem solve.
#[derive(Debug)]
pub(crate) struct PreparedSolveOutcome {
    /// Accepted physical/log-mole solution.
    pub(crate) solution: EquilibriumSolution,
    /// Backend cascade evidence for the accepted candidate.
    pub(crate) solve_report: EquilibriumSolveReport,
    /// Optional independent K_eq diagnostic, never implicit in the runner.
    pub(crate) keq_validation_status: Option<EquilibriumConstantCrossValidationStatus>,
    /// Time spent in the independent validation stage after backend solving.
    pub(crate) validation_duration: Duration,
    /// Optional explicit multi-start evidence for this accepted result.
    pub(crate) multi_start_report: Option<MultiStartSolveReport>,
}

/// Immutable fixed-formulation numerical runner.
///
/// The runner deliberately has no published mutable solution fields and no
/// phase-transition state. A failed solve simply returns an error; it cannot
/// leave a partially updated solver object behind. Active-set orchestration is
/// a separate concern and remains on its compatibility path until its own
/// immutable transition state is extracted.
#[derive(Clone)]
pub(crate) struct PreparedEquilibriumRunner {
    prepared: PreparedEquilibriumProblem,
    symbolic_standard_gibbs: Vec<Expr>,
    settings: EquilibriumSolverSettings,
}

impl PreparedEquilibriumRunner {
    /// Creates a runner from one validated prepared problem and its optional
    /// symbolic standard-state snapshot.
    pub(crate) fn new(
        prepared: PreparedEquilibriumProblem,
        symbolic_standard_gibbs: Vec<Expr>,
    ) -> Result<Self, ReactionExtentError> {
        let species_count = prepared.problem().species().len();
        if !symbolic_standard_gibbs.is_empty() && symbolic_standard_gibbs.len() != species_count {
            return Err(ReactionExtentError::DimensionMismatch(format!(
                "symbolic Gibbs snapshot has {} entries for {species_count} species",
                symbolic_standard_gibbs.len()
            )));
        }
        Ok(Self {
            prepared,
            symbolic_standard_gibbs,
            settings: EquilibriumSolverSettings::default(),
        })
    }

    /// Mutates numerical policy only; prepared chemistry remains immutable.
    pub(crate) fn configure(&mut self) -> &mut EquilibriumSolverSettings {
        &mut self.settings
    }

    /// Solves the prepared fixed formulation without constructing the mutable
    /// legacy orchestration object.
    pub(crate) fn solve(self) -> Result<PreparedSolveOutcome, ReactionExtentError> {
        let initial_guess = self
            .prepared
            .problem()
            .initial_log_moles()
            .as_slice()
            .to_vec();
        self.solve_with_initial_guess(initial_guess, None)
    }

    /// Solves from a previous accepted point of a continuation sweep.
    pub(crate) fn solve_from_seed(
        self,
        seed: LogMolesInitialGuess,
    ) -> Result<PreparedSolveOutcome, ReactionExtentError> {
        self.solve_with_initial_guess(seed.as_slice().to_vec(), None)
    }

    /// Solves from a continuation seed while borrowing a reusable RST
    /// symbolic problem. The symbolic equations and generated closures stay
    /// alive across points; only their shared temperature parameter changes.
    pub(crate) fn solve_from_seed_with_rst(
        self,
        seed: LogMolesInitialGuess,
        rst_problem: &RstPreparedProblem,
    ) -> Result<PreparedSolveOutcome, ReactionExtentError> {
        self.solve_with_initial_guess(seed.as_slice().to_vec(), Some(rst_problem))
    }

    /// Runs the same prepared formulation from ordered deterministic seeds and
    /// keeps the best accepted candidate. This is an explicit recovery tool,
    /// not an automatic replacement for the backend cascade.
    pub(crate) fn solve_from_initial_guesses(
        self,
        seeds: Vec<LogMolesInitialGuess>,
    ) -> Result<PreparedSolveOutcome, ReactionExtentError> {
        self.solve_from_initial_guesses_with_optional_rst(seeds, None)
    }

    /// Runs explicit multi-start recovery while borrowing a prepared RST
    /// symbolic problem from a higher-level continuation workflow.
    ///
    /// Only the initial log-mole vector changes between starts. Reusing the
    /// caller-owned symbolic problem avoids rebuilding the equation graph for
    /// every `P,H` temperature trial and preserves the same candidate
    /// comparison contract as the standalone method above.
    pub(crate) fn solve_from_initial_guesses_with_rst(
        self,
        seeds: Vec<LogMolesInitialGuess>,
        rst_problem: &RstPreparedProblem,
    ) -> Result<PreparedSolveOutcome, ReactionExtentError> {
        self.solve_from_initial_guesses_with_optional_rst(seeds, Some(rst_problem))
    }

    /// Core multi-start dispatch shared by the public RST and non-RST paths.
    ///
    /// Iterates over the supplied initial seeds, solves each one via
    /// [`solve_with_initial_guess`], and selects the best candidate using
    /// [`select_preferred_candidate_index`]. The optional external RST problem
    /// is forwarded to every seed attempt so the symbolic problem is reused
    /// across starts rather than rebuilt for each one. Returns an error when
    /// all seeds fail, retaining the first error for diagnostics.
    fn solve_from_initial_guesses_with_optional_rst(
        self,
        seeds: Vec<LogMolesInitialGuess>,
        external_rst_problem: Option<&RstPreparedProblem>,
    ) -> Result<PreparedSolveOutcome, ReactionExtentError> {
        if seeds.is_empty() {
            return Err(ReactionExtentError::InvalidProblem {
                field: "multi_start_seeds",
                message: "at least one initial seed is required".to_string(),
            });
        }
        let expected_dimension = self.prepared.problem().initial_moles().len();
        if let Some((seed_index, seed)) = seeds
            .iter()
            .enumerate()
            .find(|(_, seed)| seed.as_slice().len() != expected_dimension)
        {
            return Err(ReactionExtentError::InvalidProblem {
                field: "multi_start_seed_dimension",
                message: format!(
                    "seed {seed_index} has dimension {}, expected {expected_dimension}",
                    seed.as_slice().len()
                ),
            });
        }

        // Prepare one symbolic RST problem for the whole seed batch. A
        // multi-start recovery changes only the initial point; rebuilding the
        // equation graph for every seed would obscure the cost of the actual
        // numerical recovery and defeat prepared-problem reuse.
        self.settings.validate()?;
        let policy = self.settings.solver_policy.clone().unwrap_or_else(|| {
            if self.symbolic_standard_gibbs.is_empty() {
                SolverPolicy::legacy_default(self.settings.solver)
            } else {
                SolverPolicy::production_default(self.settings.solver)
            }
        });
        let rst_problem = if external_rst_problem.is_none()
            && policy
                .ordered_backends()
                .iter()
                .any(|backend| matches!(backend, SolverBackend::RustedSciThe(_)))
        {
            Some(prepare_rst_symbolic_problem_from_prepared(
                &self.prepared,
                &self.symbolic_standard_gibbs,
            )?)
        } else {
            None
        };

        let mut attempts = Vec::with_capacity(seeds.len());
        let mut best: Option<(usize, PreparedSolveOutcome)> = None;
        let mut first_error = None;

        for (start_index, seed) in seeds.into_iter().enumerate() {
            if let Some(control) = &self.settings.execution_control {
                control.check_cancelled()?;
            }
            match self.clone().solve_with_initial_guess(
                seed.as_slice().to_vec(),
                external_rst_problem.or(rst_problem.as_ref()),
            ) {
                Ok(outcome) => {
                    let residual = outcome.solution.validation().residual_l2_norm;
                    attempts.push(MultiStartAttemptReport {
                        start_index,
                        accepted: true,
                        residual_l2_norm: Some(format!("{residual:.17e}")),
                        error: None,
                        started_backend_attempts: outcome.solve_report.started_attempt_count(),
                        nonlinear_iterations: outcome.solve_report.nonlinear_iterations(),
                    });
                    let replace = best.as_ref().is_none_or(|(_, current)| {
                        crate::Thermodynamics::ChemEquilibrium::equilibrium_validation::
                            compare_candidate_reports(
                                outcome.solution.validation(),
                                current.solution.validation(),
                            )
                                == std::cmp::Ordering::Less
                    });
                    if replace {
                        best = Some((start_index, outcome));
                    }
                }
                Err(error) => {
                    if first_error.is_none() {
                        first_error = Some(error.to_string());
                    }
                    attempts.push(MultiStartAttemptReport {
                        start_index,
                        accepted: false,
                        residual_l2_norm: None,
                        error: Some(error.to_string()),
                        started_backend_attempts: started_backend_attempts_from_error(&error),
                        nonlinear_iterations: nonlinear_iterations_from_error(&error),
                    });
                }
            }
        }

        let Some((selected_start, mut outcome)) = best else {
            return Err(ReactionExtentError::InvalidProblem {
                field: "multi_start_solve",
                message: format!(
                    "all {} initial seeds failed; first error: {}",
                    attempts.len(),
                    first_error.unwrap_or_else(|| "unknown multi-start failure".to_string())
                ),
            });
        };
        outcome.multi_start_report = Some(MultiStartSolveReport {
            attempts,
            selected_start,
        });
        Ok(outcome)
    }

    /// Solves one prepared equilibrium problem from a single initial guess.
    ///
    /// Validates solver settings, checks for cancellation, resolves the solver
    /// policy and cascade budget, then dispatches through the common backend
    /// adapter. The optional RST symbolic problem is forwarded when the policy
    /// selects a RustedSciThe backend. The result includes the accepted
    /// solution, validation report, backend cascade evidence, and optional
    /// equilibrium-constant cross-validation.
    fn solve_with_initial_guess(
        self,
        initial_guess: Vec<f64>,
        external_rst_problem: Option<&RstPreparedProblem>,
    ) -> Result<PreparedSolveOutcome, ReactionExtentError> {
        self.settings.validate()?;
        if let Some(control) = &self.settings.execution_control {
            control.check_cancelled()?;
        }
        let policy = self.settings.solver_policy.clone().unwrap_or_else(|| {
            if self.symbolic_standard_gibbs.is_empty() {
                SolverPolicy::legacy_default(self.settings.solver)
            } else {
                SolverPolicy::production_default(self.settings.solver)
            }
        });
        let backends = policy.ordered_backends();
        if backends.is_empty() {
            return Err(ReactionExtentError::InvalidProblem {
                field: "solver_policy",
                message: "a solver policy must contain at least one backend".to_string(),
            });
        }

        let execution_control = self.settings.execution_control.clone();
        let residual = |values: &[f64]| {
            if let Some(control) = &execution_control {
                control.check_cancelled()?;
            }
            self.prepared.residual(values)
        };
        let scale = if self.settings.scaling_flag {
            Some(self.prepared.residual_scaling_contract()?.scale().to_vec())
        } else {
            None
        };
        let scaled_residual = |values: &[f64]| match &scale {
            Some(scale) => self.prepared.scaled_residual(values, scale),
            None => residual(values),
        };
        let uses_legacy = backends
            .iter()
            .any(|backend| matches!(backend, SolverBackend::Legacy(_)));
        let jacobian = if uses_legacy {
            Some(|values: &[f64]| {
                if let Some(scale) = &scale {
                    self.prepared.scaled_jacobian(values, scale)
                } else {
                    self.prepared.jacobian(values)
                }
            })
        } else {
            None
        };
        let owned_rst_problem = if external_rst_problem.is_none()
            && backends
                .iter()
                .any(|backend| matches!(backend, SolverBackend::RustedSciThe(_)))
        {
            Some(prepare_rst_symbolic_problem_from_prepared(
                &self.prepared,
                &self.symbolic_standard_gibbs,
            )?)
        } else {
            None
        };
        let rst_problem = external_rst_problem.or(owned_rst_problem.as_ref());
        let feasible = |candidate: &[f64]| {
            compute_species_moles(candidate)
                .map(|moles| {
                    moles
                        .iter()
                        .all(|&moles_i| moles_i >= -LEGACY_MOLE_FEASIBILITY_TOLERANCE)
                })
                .unwrap_or(false)
        };
        let solver_tolerance = self.settings.solver_params.tol;
        let element_tolerance = 10.0 * solver_tolerance.max(1e-8);
        let validate_candidate = |candidate: &[f64]| {
            let raw_residual = residual(candidate)?;
            let acceptance_residual = scaled_residual(candidate)?;
            validate_equilibrium_candidate(
                EquilibriumCandidateResiduals {
                    log_moles: candidate,
                    raw_residual: &raw_residual,
                    acceptance_residual: &acceptance_residual,
                },
                EquilibriumAcceptanceCriteria::new(
                    solver_tolerance,
                    element_tolerance,
                    solver_tolerance,
                )?
                .with_element_balance_relative_tolerance(element_tolerance)?,
                self.prepared.problem().element_composition(),
                self.prepared.element_totals(),
            )
        };
        let budget = self.settings.solver_budget.unwrap_or_else(|| {
            SolverCascadeBudget::new(
                backends.len(),
                self.settings.solver_params.max_iter,
                self.settings
                    .solver_params
                    .max_iter
                    .saturating_mul(backends.len()),
            )
        });
        let backend_refs: Vec<&dyn EquilibriumNonlinearBackend> = backends
            .iter()
            .map(|backend| backend as &dyn EquilibriumNonlinearBackend)
            .collect();
        let (log_moles, validation, solve_report) =
            EquilibriumLogMoles::solve_backend_cascade_with_control(
                &backend_refs,
                initial_guess,
                &scaled_residual,
                jacobian.as_ref().map(|function| {
                    function as &dyn Fn(&[f64]) -> Result<DMatrix<f64>, ReactionExtentError>
                }),
                &feasible,
                &validate_candidate,
                policy,
                budget,
                &self.settings.solver_params,
                rst_problem,
                execution_control.as_ref(),
            )?;
        if let Some(control) = &execution_control {
            control.check_cancelled()?;
        }
        let solution = self.prepared.accepted_solution(log_moles, validation)?;
        let validation_started = Instant::now();
        let keq_validation_status =
            run_keq_cross_validation(&self.prepared, &solution, &self.settings)?;
        Ok(PreparedSolveOutcome {
            solution,
            solve_report,
            keq_validation_status,
            validation_duration: validation_started.elapsed(),
            multi_start_report: None,
        })
    }
}

fn started_backend_attempts_from_error(error: &ReactionExtentError) -> usize {
    match error {
        ReactionExtentError::AllBackendsFailed { attempts }
        | ReactionExtentError::CascadeAborted { attempts, .. } => attempts
            .iter()
            .filter(|attempt| attempt.is_started())
            .count(),
        _ => 0,
    }
}

fn nonlinear_iterations_from_error(error: &ReactionExtentError) -> usize {
    match error {
        ReactionExtentError::AllBackendsFailed { attempts }
        | ReactionExtentError::CascadeAborted { attempts, .. } => attempts
            .iter()
            .filter(|attempt| attempt.is_started())
            .filter_map(|attempt| attempt.metrics.as_ref())
            .map(|metrics| metrics.iterations)
            .sum(),
        _ => 0,
    }
}

fn run_keq_cross_validation(
    prepared: &PreparedEquilibriumProblem,
    candidate: &EquilibriumSolution,
    settings: &EquilibriumSolverSettings,
) -> Result<Option<EquilibriumConstantCrossValidationStatus>, ReactionExtentError> {
    let mode = settings.keq_validation_mode;
    if mode == EquilibriumConstantValidationMode::Off {
        return Ok(None);
    }
    let applicable = prepared.problem().phases().len() == 1
        && matches!(
            prepared.problem().phases()[0].kind,
            PhaseActivityModel::IdealGas
        );
    if !applicable {
        let message =
            "independent validation currently supports exactly one ideal-gas phase".to_string();
        return match mode {
            EquilibriumConstantValidationMode::Required => {
                Err(ReactionExtentError::ValidationNotApplicable {
                    path: "equilibrium_constant_cross_validation",
                    message,
                })
            }
            EquilibriumConstantValidationMode::WhenApplicable => Ok(Some(
                EquilibriumConstantCrossValidationStatus::ValidatorNotApplicable { message },
            )),
            EquilibriumConstantValidationMode::Off => Ok(None),
        };
    }

    let problem =
        EquilibriumConstantProblem::from_prepared_ideal_gas(prepared, Default::default())?;
    let solver = EquilibriumConstantSolver {
        mode: match mode {
            EquilibriumConstantValidationMode::Off => EquilibriumConstantSolverMode::Off,
            EquilibriumConstantValidationMode::WhenApplicable => {
                EquilibriumConstantSolverMode::WhenApplicable
            }
            EquilibriumConstantValidationMode::Required => EquilibriumConstantSolverMode::Required,
        },
        validation_tolerances: EquilibriumConstantValidationTolerances::default(),
        ..Default::default()
    };
    let validator = solver.solve_if_applicable(&problem);
    let status = classify_equilibrium_constant_cross_validation(
        &problem,
        Ok(candidate.clone()),
        validator,
        settings.keq_validation_tolerances,
    )?;
    if mode == EquilibriumConstantValidationMode::Required
        && !matches!(
            status,
            EquilibriumConstantCrossValidationStatus::Compared(ref report) if report.accepted
        )
    {
        return Err(ReactionExtentError::InvalidCandidate {
            field: "equilibrium_constant_cross_validation",
            message: status.to_string(),
        });
    }
    Ok(Some(status))
}
