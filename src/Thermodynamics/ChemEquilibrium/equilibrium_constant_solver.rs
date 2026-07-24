//! Independent small-system solver for equilibrium-constant validation.
//!
//! This path is intentionally narrower than the canonical equilibrium solver.
//! It works in reaction-extent space, supports only tiny one-reaction systems
//! in this first pass, and reports `ln(Q) - ln(K)` evidence independently from
//! the canonical log-moles residual pipeline.

use crate::Thermodynamics::ChemEquilibrium::equilibrium_constant_problem::EquilibriumConstantProblem;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_constant_validation::{
    EquilibriumConstantValidationReport, EquilibriumConstantValidationTolerances,
    validate_equilibrium_constants,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_nonlinear::ReactionExtentError;
use std::fmt;

/// Solver settings for the independent reaction-extent validator.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EquilibriumConstantSolverSettings {
    /// Maximum safeguarded Newton iterations.
    pub max_iterations: usize,
    /// Absolute tolerance on the scalar `ln(Q) - ln(K)` residual.
    pub residual_tolerance: f64,
    /// Initial finite-difference step used for the scalar derivative.
    pub finite_difference_step: f64,
    /// Minimum relative distance from the feasibility boundaries.
    pub feasibility_margin: f64,
    /// Number of samples used to search for a sign-change bracket.
    pub bracket_samples: usize,
}

impl Default for EquilibriumConstantSolverSettings {
    fn default() -> Self {
        Self {
            max_iterations: 64,
            residual_tolerance: 1e-10,
            finite_difference_step: 1e-6,
            feasibility_margin: 1e-12,
            bracket_samples: 9,
        }
    }
}

impl EquilibriumConstantSolverSettings {
    fn validate(self) -> Result<Self, ReactionExtentError> {
        if self.max_iterations == 0 {
            return Err(invalid_solver("max_iterations must be strictly positive"));
        }
        for (field, value) in [
            ("residual_tolerance", self.residual_tolerance),
            ("finite_difference_step", self.finite_difference_step),
            ("feasibility_margin", self.feasibility_margin),
        ] {
            if !value.is_finite() || value <= 0.0 {
                return Err(invalid_solver(format!(
                    "{field} must be finite and strictly positive"
                )));
            }
        }
        if self.bracket_samples < 3 {
            return Err(invalid_solver("bracket_samples must be at least 3"));
        }
        Ok(self)
    }
}

/// Solver policy for the independent validator.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum EquilibriumConstantSolverMode {
    /// Disable the independent extent solver.
    #[default]
    Off,
    /// Run the small-system validator when the problem fits its contract.
    WhenApplicable,
    /// Treat a non-applicable validator as an error.
    Required,
}

/// Summary of one independent extent solve.
#[derive(Debug, Clone, PartialEq)]
pub struct EquilibriumConstantSolveReport {
    /// Whether the extent solver converged.
    pub converged: bool,
    /// Number of safeguarded iterations performed.
    pub iterations: usize,
    /// Reaction extent used to reconstruct the candidate composition.
    pub extent: f64,
    /// Scalar residual `ln(Q) - ln(K)` at the accepted candidate.
    pub log_residual: f64,
    /// Lower feasibility bound implied by mole positivity.
    pub lower_bound: f64,
    /// Upper feasibility bound implied by mole positivity.
    pub upper_bound: f64,
    /// Whether the solve used a sign-change bracket.
    pub bracketed: bool,
}

/// Accepted independent validation result.
#[derive(Debug, Clone, PartialEq)]
pub struct EquilibriumConstantSolveResult {
    /// Accepted composition extent.
    pub extent: f64,
    /// Reconstructed physical mole numbers.
    pub moles: Vec<f64>,
    /// Validation report generated after the extent solver accepted a candidate.
    pub validation: EquilibriumConstantValidationReport,
    /// Solver report with convergence and feasibility evidence.
    pub report: EquilibriumConstantSolveReport,
}

/// One stable summary row for CLI output and snapshot tests.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct EquilibriumConstantSolveRow {
    /// Logical section name.
    pub section: &'static str,
    /// Stable row label.
    pub label: String,
    /// Human-readable row value.
    pub value: String,
}

impl fmt::Display for EquilibriumConstantSolveRow {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "[{}] {} = {}", self.section, self.label, self.value)
    }
}

impl EquilibriumConstantSolveReport {
    /// Stable summary rows for CLI output and snapshot tests.
    pub fn summary_rows(&self) -> Vec<EquilibriumConstantSolveRow> {
        vec![
            EquilibriumConstantSolveRow {
                section: "solve",
                label: "converged".to_string(),
                value: self.converged.to_string(),
            },
            EquilibriumConstantSolveRow {
                section: "solve",
                label: "iterations".to_string(),
                value: self.iterations.to_string(),
            },
            EquilibriumConstantSolveRow {
                section: "solve",
                label: "extent".to_string(),
                value: format!("{:.6e}", self.extent),
            },
            EquilibriumConstantSolveRow {
                section: "solve",
                label: "log_residual".to_string(),
                value: format!("{:.6e}", self.log_residual),
            },
            EquilibriumConstantSolveRow {
                section: "solve",
                label: "lower_bound".to_string(),
                value: format!("{:.6e}", self.lower_bound),
            },
            EquilibriumConstantSolveRow {
                section: "solve",
                label: "upper_bound".to_string(),
                value: format!("{:.6e}", self.upper_bound),
            },
            EquilibriumConstantSolveRow {
                section: "solve",
                label: "bracketed".to_string(),
                value: self.bracketed.to_string(),
            },
        ]
    }
}

impl fmt::Display for EquilibriumConstantSolveReport {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for row in self.summary_rows() {
            writeln!(f, "{row}")?;
        }
        Ok(())
    }
}

/// Independent solver for very small equilibrium-constant problems.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EquilibriumConstantSolver {
    /// Numeric controls for safeguarded Newton iterations.
    pub settings: EquilibriumConstantSolverSettings,
    /// Numerical acceptance gate for the independent validation report.
    pub validation_tolerances: EquilibriumConstantValidationTolerances,
    /// Whether the solver should be run at all when requested through policy.
    pub mode: EquilibriumConstantSolverMode,
}

impl Default for EquilibriumConstantSolver {
    fn default() -> Self {
        Self {
            settings: EquilibriumConstantSolverSettings::default(),
            validation_tolerances: EquilibriumConstantValidationTolerances::default(),
            mode: EquilibriumConstantSolverMode::WhenApplicable,
        }
    }
}

impl EquilibriumConstantSolver {
    /// Solves only when the validator is applicable; otherwise returns `Ok(None)`
    /// for the non-required modes.
    pub fn solve_if_applicable(
        &self,
        problem: &EquilibriumConstantProblem,
    ) -> Result<Option<EquilibriumConstantSolveResult>, ReactionExtentError> {
        if self.mode == EquilibriumConstantSolverMode::Off {
            return Ok(None);
        }
        match self.solve(problem) {
            Ok(result) => Ok(Some(result)),
            Err(ReactionExtentError::ValidationNotApplicable { .. }) => match self.mode {
                EquilibriumConstantSolverMode::Required => {
                    Err(ReactionExtentError::ValidationNotApplicable {
                        path: "equilibrium_constant_solver",
                        message: "validator is required but not applicable".to_string(),
                    })
                }
                EquilibriumConstantSolverMode::Off
                | EquilibriumConstantSolverMode::WhenApplicable => Ok(None),
            },
            Err(error) => Err(error),
        }
    }

    /// Solves a very small equilibrium-constant problem in reaction-extent space.
    pub fn solve(
        &self,
        problem: &EquilibriumConstantProblem,
    ) -> Result<EquilibriumConstantSolveResult, ReactionExtentError> {
        let settings = self.settings.validate()?;
        if self.mode == EquilibriumConstantSolverMode::Off {
            return Err(ReactionExtentError::ValidationNotApplicable {
                path: "equilibrium_constant_solver",
                message: "validation mode is disabled".to_string(),
            });
        }
        if problem.basis().reaction_count() != 1 {
            return self
                .not_applicable("solver currently supports exactly one independent reaction");
        }
        if problem.basis().species().len() != problem.initial_moles().len() {
            return Err(ReactionExtentError::DimensionMismatch(
                "basis species ordering and initial mole vector disagree".to_string(),
            ));
        }

        let reaction = problem.basis().reaction_id(0)?;
        let coefficients = problem.basis().reaction(reaction)?;
        let (lower_bound, upper_bound) =
            feasible_extent_bounds(problem.initial_moles(), &coefficients)?;
        let feasibility_shrink =
            settings.feasibility_margin * (upper_bound - lower_bound).abs().max(1.0);
        let lower_bound = lower_bound + feasibility_shrink;
        let upper_bound = upper_bound - feasibility_shrink;
        if lower_bound >= upper_bound {
            return Err(ReactionExtentError::InvalidProblem {
                field: "equilibrium_constant_solver",
                message: "feasibility margin removed the positive interval".to_string(),
            });
        }

        let samples = sample_points(problem, lower_bound, upper_bound, settings.bracket_samples)?;
        let bracket = find_bracket(&samples)?;
        let bracketed = bracket.is_some();
        let (lower, upper, mut current) = match bracket {
            Some((lo, hi)) => (lo, hi, 0.5 * (lo + hi)),
            None => (
                lower_bound,
                upper_bound,
                best_sample(&samples).ok_or_else(|| {
                    ReactionExtentError::ValidationNotApplicable {
                        path: "equilibrium_constant_solver",
                        message: "could not sample a feasible interior point".to_string(),
                    }
                })?,
            ),
        };
        let mut current_residual = residual_at(problem, current)?;

        for iteration in 0..settings.max_iterations {
            if current_residual.abs() <= settings.residual_tolerance {
                return self.accept(
                    problem,
                    current,
                    current_residual,
                    lower_bound,
                    upper_bound,
                    bracketed,
                    iteration + 1,
                );
            }

            let derivative = finite_difference_derivative(
                problem,
                current,
                lower_bound,
                upper_bound,
                settings.finite_difference_step,
            )?;
            let mut candidate = if derivative.is_finite() && derivative.abs() > 0.0 {
                current - current_residual / derivative
            } else {
                0.5 * (lower + upper)
            };

            if !candidate.is_finite() || candidate <= lower_bound || candidate >= upper_bound {
                candidate = 0.5 * (lower + upper);
            }

            let mut accepted_candidate = None;
            let mut alpha = 1.0_f64;
            while alpha >= 1e-8 {
                let trial = current + alpha * (candidate - current);
                if trial > lower_bound && trial < upper_bound {
                    let trial_residual = residual_at(problem, trial)?;
                    if trial_residual.abs() <= current_residual.abs() {
                        accepted_candidate = Some((trial, trial_residual));
                        break;
                    }
                }
                alpha *= 0.5;
            }

            let (next, next_residual) = accepted_candidate.unwrap_or_else(|| {
                let midpoint = 0.5 * (lower + upper);
                let residual = residual_at(problem, midpoint).unwrap_or(current_residual);
                (midpoint, residual)
            });

            current = next;
            current_residual = next_residual;
        }

        Err(ReactionExtentError::SolveError(
            crate::Thermodynamics::ChemEquilibrium::equilibrium_nonlinear::SolveError::MaxIterations,
        ))
    }

    fn accept(
        &self,
        problem: &EquilibriumConstantProblem,
        extent: f64,
        log_residual: f64,
        lower_bound: f64,
        upper_bound: f64,
        bracketed: bool,
        iterations: usize,
    ) -> Result<EquilibriumConstantSolveResult, ReactionExtentError> {
        let moles = reconstruct_moles(problem, extent)?;
        let validation =
            validate_equilibrium_constants(problem, &moles, self.validation_tolerances)?;
        Ok(EquilibriumConstantSolveResult {
            extent,
            moles,
            validation,
            report: EquilibriumConstantSolveReport {
                converged: true,
                iterations,
                extent,
                log_residual,
                lower_bound,
                upper_bound,
                bracketed,
            },
        })
    }

    fn not_applicable<T>(&self, message: &str) -> Result<T, ReactionExtentError> {
        let message = match self.mode {
            EquilibriumConstantSolverMode::Off => "validation mode is disabled".to_string(),
            EquilibriumConstantSolverMode::WhenApplicable
            | EquilibriumConstantSolverMode::Required => message.to_string(),
        };
        Err(ReactionExtentError::ValidationNotApplicable {
            path: "equilibrium_constant_solver",
            message,
        })
    }
}

fn feasible_extent_bounds(
    initial_moles: &[f64],
    coefficients: &[f64],
) -> Result<(f64, f64), ReactionExtentError> {
    if initial_moles.len() != coefficients.len() {
        return Err(ReactionExtentError::DimensionMismatch(format!(
            "initial moles have {} entries but the reaction basis has {} species",
            initial_moles.len(),
            coefficients.len()
        )));
    }

    let mut lower = f64::NEG_INFINITY;
    let mut upper = f64::INFINITY;
    for (index, (&mole, &nu)) in initial_moles.iter().zip(coefficients.iter()).enumerate() {
        if !mole.is_finite() || mole < 0.0 {
            return Err(ReactionExtentError::InvalidProblem {
                field: "initial_moles",
                message: format!("entry {index} must be finite and non-negative"),
            });
        }
        if nu > 0.0 {
            lower = lower.max(-mole / nu);
        } else if nu < 0.0 {
            upper = upper.min(mole / (-nu));
        }
    }

    if !lower.is_finite() || !upper.is_finite() || lower >= upper {
        return Err(ReactionExtentError::InvalidProblem {
            field: "equilibrium_constant_solver",
            message: "no positive-feasible reaction-extent interval was found".to_string(),
        });
    }

    Ok((lower, upper))
}

fn reconstruct_moles(
    problem: &EquilibriumConstantProblem,
    extent: f64,
) -> Result<Vec<f64>, ReactionExtentError> {
    let reaction = problem.basis().reaction_id(0)?;
    let coefficients = problem.basis().reaction(reaction)?;
    let mut moles = Vec::with_capacity(problem.initial_moles().len());
    for (index, (&mole, &nu)) in problem
        .initial_moles()
        .iter()
        .zip(coefficients.iter())
        .enumerate()
    {
        let value = mole + extent * nu;
        if !value.is_finite() || value <= 0.0 {
            return Err(ReactionExtentError::InvalidCandidate {
                field: "candidate_moles",
                message: format!("species {index} is non-positive at extent {extent}"),
            });
        }
        moles.push(value);
    }
    Ok(moles)
}

fn residual_at(
    problem: &EquilibriumConstantProblem,
    extent: f64,
) -> Result<f64, ReactionExtentError> {
    let moles = reconstruct_moles(problem, extent)?;
    let reaction = problem.basis().reaction_id(0)?;
    Ok(problem.ln_reaction_quotient(reaction, &moles)?
        - problem.ln_equilibrium_constant(reaction)?)
}

fn finite_difference_derivative(
    problem: &EquilibriumConstantProblem,
    extent: f64,
    lower_bound: f64,
    upper_bound: f64,
    step: f64,
) -> Result<f64, ReactionExtentError> {
    let feasible_span = upper_bound - lower_bound;
    let h = (step * extent.abs().max(1.0)).min(0.25 * feasible_span);
    let h = h.max(feasible_span * 1e-8);
    let left = (extent - h).max(lower_bound + 0.5 * h);
    let right = (extent + h).min(upper_bound - 0.5 * h);
    if !(left < extent && extent < right) {
        return Err(ReactionExtentError::ValidationNotApplicable {
            path: "equilibrium_constant_solver",
            message: "could not form a safe finite-difference stencil".to_string(),
        });
    }
    let f_left = residual_at(problem, left)?;
    let f_right = residual_at(problem, right)?;
    let denom = right - left;
    if denom <= 0.0 {
        return Err(ReactionExtentError::ValidationNotApplicable {
            path: "equilibrium_constant_solver",
            message: "finite-difference stencil collapsed".to_string(),
        });
    }
    Ok((f_right - f_left) / denom)
}

fn sample_points(
    problem: &EquilibriumConstantProblem,
    lower_bound: f64,
    upper_bound: f64,
    samples: usize,
) -> Result<Vec<(f64, f64)>, ReactionExtentError> {
    let step = (upper_bound - lower_bound) / (samples as f64 + 1.0);
    let mut points = Vec::with_capacity(samples);
    for index in 0..samples {
        let extent = lower_bound + (index as f64 + 1.0) * step;
        if extent <= lower_bound || extent >= upper_bound {
            continue;
        }
        let residual = residual_at(problem, extent)?;
        points.push((extent, residual));
    }
    Ok(points)
}

fn best_sample(points: &[(f64, f64)]) -> Option<f64> {
    points
        .iter()
        .min_by(|lhs, rhs| lhs.1.abs().total_cmp(&rhs.1.abs()))
        .map(|(extent, _)| *extent)
}

fn find_bracket(points: &[(f64, f64)]) -> Result<Option<(f64, f64)>, ReactionExtentError> {
    for pair in points.windows(2) {
        let (left_extent, left_residual) = pair[0];
        let (right_extent, right_residual) = pair[1];
        if left_residual == 0.0 {
            return Ok(Some((left_extent, left_extent)));
        }
        if left_residual.signum() != right_residual.signum() {
            return Ok(Some((left_extent, right_extent)));
        }
    }

    Ok(None)
}

fn invalid_solver(message: impl Into<String>) -> ReactionExtentError {
    ReactionExtentError::InvalidProblem {
        field: "equilibrium_constant_solver",
        message: message.into(),
    }
}
