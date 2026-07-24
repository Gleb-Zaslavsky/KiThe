//! Validation mode, tolerances, and per-reaction evidence for K_eq cross-validation.
//!
//! # Purpose
//!
//! This module defines the **validation boundary** for the independent
//! equilibrium-constant (K_eq) cross-validation system. It provides:
//!
//! - [`EquilibriumConstantValidationMode`] — whether and how to run validation.
//! - [`EquilibriumConstantValidationTolerances`] — numerical acceptance criteria.
//! - [`EquilibriumConstantValidationEvidence`] — per-reaction evidence from the
//!   independent solver.
//!
//! The actual cross-validation logic (comparing canonical and independent solutions)
//! lives in [`equilibrium_constant_cross_validation`](super::equilibrium_constant_cross_validation).
//! This module provides the **input configuration** and **output evidence** types.
//!
//! # Key Structures
//!
//! | Structure | Role |
//! |-----------|------|
//! | [`EquilibriumConstantValidationMode`] | Enum: Off, WhenApplicable, Required |
//! | [`EquilibriumConstantValidationTolerances`] | Acceptance criteria (max_log_residual) |
//! | [`EquilibriumConstantValidationEvidence`] | Per-reaction validation result |
//!
//! # Dataflow
//!
//! ```text
//!   EquilibriumSolverSettings
//!     ├── keq_validation_mode: EquilibriumConstantValidationMode
//!     │     ├── Off → skip validation entirely
//!     │     ├── WhenApplicable → validate if possible, else accept
//!     │     └── Required → validate, fail if not applicable
//!     │
//!     └── keq_tolerances: EquilibriumConstantValidationTolerances
//!           └── max_log_residual: f64
//!     │
//!     v
//!   EquilibriumLogMoles::run_keq_cross_validation()
//!     ├── Checks mode → Off? Return early
//!     ├── Builds EquilibriumConstantProblem
//!     ├── Runs EquilibriumConstantSolver
//!     ├── Compares solutions via compare_equilibrium_constant_solutions()
//!     └── Returns EquilibriumConstantCrossValidationStatus
//! ```
//!
//! # Examples
//!
//! ```rust, ignore
//! use equilibrium_constant_validation::*;
//!
//! let mode = EquilibriumConstantValidationMode::WhenApplicable;
//! let tolerances = EquilibriumConstantValidationTolerances {
//!     max_log_residual: 1e-8,
//! };
//! assert!(tolerances.validate().is_ok());
//! ```
//!
//! # Non-obvious Details
//!
//! - `EquilibriumConstantValidationMode::Required` is **currently unrealizable**
//!   for many systems (see SourceCraft Diagnostics item B.3 in TODO_ANALYSIS.md).
//!   The independent K_eq solver only supports single-reaction systems.
//! - `EquilibriumConstantValidationTolerances::validate()` rejects non-finite
//!   and non-positive `max_log_residual` values.
//!
//! # Related Modules
//!
//! - [`equilibrium_constant_cross_validation`](super::equilibrium_constant_cross_validation) — cross-validation logic
//! - [`equilibrium_constant_problem`](super::equilibrium_constant_problem) — K_eq problem definition
//! - [`equilibrium_constant_solver`](super::equilibrium_constant_solver) — K_eq solver
//! - [`equilibrium_log_moles`](super::equilibrium_log_moles) — orchestrator that calls validation
//!

use crate::Thermodynamics::ChemEquilibrium::equilibrium_constant_problem::EquilibriumConstantProblem;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_ids::ReactionId;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_nonlinear::ReactionExtentError;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_problem::EquilibriumConditions;
use std::fmt;

/// Controls whether an independent K_eq validation is requested.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum EquilibriumConstantValidationMode {
    /// Do not construct or run the independent validator.
    #[default]
    Off,
    /// Validate supported small systems and report unsupported systems normally.
    WhenApplicable,
    /// Treat a non-applicable validator as a validation failure.
    Required,
}

/// Numerical acceptance criteria owned by the independent validator.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EquilibriumConstantValidationTolerances {
    /// Largest accepted absolute value of `ln(Q) - ln(K)`.
    pub max_log_residual: f64,
}

impl Default for EquilibriumConstantValidationTolerances {
    fn default() -> Self {
        Self {
            max_log_residual: 1e-8,
        }
    }
}

impl EquilibriumConstantValidationTolerances {
    /// Rejects invalid tolerances before candidate evaluation.
    pub fn validate(self) -> Result<Self, ReactionExtentError> {
        if !self.max_log_residual.is_finite() || self.max_log_residual <= 0.0 {
            return Err(ReactionExtentError::InvalidProblem {
                field: "equilibrium_constant_tolerances",
                message: "max_log_residual must be finite and strictly positive".to_string(),
            });
        }
        Ok(self)
    }
}

/// Per-reaction evidence produced independently from the main solver residual.
#[derive(Debug, Clone, PartialEq)]
pub struct EquilibriumConstantReactionReport {
    /// Reaction identity in the validated basis.
    pub reaction: ReactionId,
    /// Natural logarithm of the candidate reaction quotient.
    pub ln_q: f64,
    /// Natural logarithm of the thermodynamic equilibrium constant.
    pub ln_k: f64,
    /// Signed difference `ln(Q) - ln(K)`.
    pub log_residual: f64,
}

/// Independent equilibrium-constant assessment of one candidate composition.
#[derive(Debug, Clone, PartialEq)]
pub struct EquilibriumConstantValidationReport {
    /// Conditions under which both quotients and constants were evaluated.
    pub conditions: EquilibriumConditions,
    /// Ordered report for every independent reaction.
    pub reactions: Vec<EquilibriumConstantReactionReport>,
    /// Largest absolute reaction log residual.
    pub max_abs_log_residual: f64,
    /// Whether the candidate satisfies the validator-owned tolerance.
    pub accepted: bool,
}

/// One stable summary row for CLI output and snapshot tests.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct EquilibriumConstantValidationRow {
    /// Logical section name.
    pub section: &'static str,
    /// Stable row label.
    pub label: String,
    /// Human-readable row value.
    pub value: String,
}

impl fmt::Display for EquilibriumConstantValidationRow {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "[{}] {} = {}", self.section, self.label, self.value)
    }
}

impl EquilibriumConstantValidationReport {
    /// Stable summary rows for CLI output and snapshot tests.
    pub fn summary_rows(&self) -> Vec<EquilibriumConstantValidationRow> {
        let mut rows = vec![
            EquilibriumConstantValidationRow {
                section: "validation",
                label: "temperature".to_string(),
                value: format!("{:.6}", self.conditions.temperature()),
            },
            EquilibriumConstantValidationRow {
                section: "validation",
                label: "pressure".to_string(),
                value: format!("{:.6}", self.conditions.pressure()),
            },
            EquilibriumConstantValidationRow {
                section: "validation",
                label: "reaction_count".to_string(),
                value: self.reactions.len().to_string(),
            },
            EquilibriumConstantValidationRow {
                section: "validation",
                label: "max_abs_log_residual".to_string(),
                value: format!("{:.6e}", self.max_abs_log_residual),
            },
            EquilibriumConstantValidationRow {
                section: "validation",
                label: "accepted".to_string(),
                value: self.accepted.to_string(),
            },
        ];

        for reaction in &self.reactions {
            rows.push(EquilibriumConstantValidationRow {
                section: "reaction",
                label: reaction.reaction.index().to_string(),
                value: format!(
                    "lnQ={:.6e}, lnK={:.6e}, residual={:.6e}",
                    reaction.ln_q, reaction.ln_k, reaction.log_residual
                ),
            });
        }

        rows
    }
}

impl fmt::Display for EquilibriumConstantValidationReport {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for row in self.summary_rows() {
            writeln!(f, "{row}")?;
        }
        Ok(())
    }
}

/// Evaluates a candidate without invoking the canonical equilibrium residual.
pub fn validate_equilibrium_constants(
    problem: &EquilibriumConstantProblem,
    candidate_moles: &[f64],
    tolerances: EquilibriumConstantValidationTolerances,
) -> Result<EquilibriumConstantValidationReport, ReactionExtentError> {
    let tolerances = tolerances.validate()?;
    let mut reactions = Vec::with_capacity(problem.basis().reaction_count());
    let mut max_abs_log_residual = 0.0_f64;

    for index in 0..problem.basis().reaction_count() {
        let reaction = problem.basis().reaction_id(index)?;
        let ln_q = problem.ln_reaction_quotient(reaction, candidate_moles)?;
        let ln_k = problem.ln_equilibrium_constant(reaction)?;
        let log_residual = ln_q - ln_k;
        max_abs_log_residual = max_abs_log_residual.max(log_residual.abs());
        reactions.push(EquilibriumConstantReactionReport {
            reaction,
            ln_q,
            ln_k,
            log_residual,
        });
    }

    Ok(EquilibriumConstantValidationReport {
        conditions: problem.conditions(),
        reactions,
        max_abs_log_residual,
        accepted: max_abs_log_residual <= tolerances.max_log_residual,
    })
}
