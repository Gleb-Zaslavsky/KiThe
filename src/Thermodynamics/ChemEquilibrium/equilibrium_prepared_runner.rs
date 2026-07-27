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
    EquilibriumConstantCrossValidationStatus, classify_equilibrium_constant_cross_validation,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_constant_problem::EquilibriumConstantProblem;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_constant_solver::{
    EquilibriumConstantSolver, EquilibriumConstantSolverMode,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_constant_validation::{
    EquilibriumConstantValidationMode, EquilibriumConstantValidationTolerances,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_log_moles::{
    EquilibriumLogMoles, EquilibriumSolverSettings, LEGACY_MOLE_FEASIBILITY_TOLERANCE,
    compute_species_moles,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_nonlinear::ReactionExtentError;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_problem::{
    EquilibriumSolution, LogMolesInitialGuess, PreparedEquilibriumProblem,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_rst_backend::{
    RstPreparedProblem, prepare_rst_symbolic_problem_from_prepared,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_solver_policy::{
    EquilibriumSolveReport, SolverBackend, SolverCascadeBudget, SolverPolicy,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_validation::{
    EquilibriumAcceptanceCriteria, EquilibriumCandidateResiduals, validate_equilibrium_candidate,
};
use RustedSciThe::symbolic::symbolic_engine::Expr;
use nalgebra::DMatrix;
use std::time::{Duration, Instant};

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
}

/// Immutable fixed-formulation numerical runner.
///
/// The runner deliberately has no published mutable solution fields and no
/// phase-transition state. A failed solve simply returns an error; it cannot
/// leave a partially updated solver object behind. Active-set orchestration is
/// a separate concern and remains on its compatibility path until its own
/// immutable transition state is extracted.
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

    fn solve_with_initial_guess(
        self,
        initial_guess: Vec<f64>,
        external_rst_problem: Option<&RstPreparedProblem>,
    ) -> Result<PreparedSolveOutcome, ReactionExtentError> {
        self.settings.validate()?;
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

        let residual = |values: &[f64]| self.prepared.residual(values);
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
        let (log_moles, validation, solve_report) = EquilibriumLogMoles::solve_backend_cascade(
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
            self.prepared.problem().initial_moles(),
            &self.prepared.reaction_basis().reactions,
            rst_problem,
        )?;
        let solution = self.prepared.accepted_solution(log_moles, validation)?;
        let validation_started = Instant::now();
        let keq_validation_status =
            run_keq_cross_validation(&self.prepared, &solution, &self.settings)?;
        Ok(PreparedSolveOutcome {
            solution,
            solve_report,
            keq_validation_status,
            validation_duration: validation_started.elapsed(),
        })
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
