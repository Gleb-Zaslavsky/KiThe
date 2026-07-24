//! Typed input boundary, scaling contracts, and solution containers for the canonical solver.
//!
//! # Purpose
//!
//! This module defines the **typed input boundary** for the equilibrium solver.
//! [`EquilibriumProblem`] is a validated, self-contained problem description that
//! is independent of mutable solver state. It ensures that all dimensional and
//! physical invariants are checked **before** any nonlinear backend receives
//! closures, matrices, or an initial iterate.
//!
//! The module also provides:
//! - [`PreparedEquilibriumProblem`] — a ready-to-solve problem with precomputed
//!   reaction basis, scaling, and residual/Jacobian closures.
//! - [`EquilibriumSolution`] — the immutable result of a successful solve.
//! - [`ResidualScalingContract`] and [`VariableScalingContract`] — row and variable
//!   scaling for numerical conditioning.
//! - [`LogMolesInitialGuess`] — typed initial guess with trace-species seeding policies.
//! - Diagnostic structures: [`FormulationDiagnostics`], [`SpeciesCapacityReport`],
//!   [`EquilibriumProblemPreview`].
//!
//! # Physical and Mathematical Background
//!
//! ## Problem Structure
//!
//! An equilibrium problem is defined by:
//!
//! - **Components**: species identities with their thermodynamic phases.
//! - **Initial moles**: `n_i^0` — the starting composition.
//! - **Element composition**: `A_{ij}` — atoms of element `j` in species `i`.
//! - **Gibbs functions**: `g_i(T)` — standard Gibbs free energy of each species.
//! - **Phases**: phase assignments and activity models.
//! - **Conditions**: temperature `T`, pressure `P`, reference pressure `P0`.
//!
//! ## Scaling
//!
//! ### Row Scaling (ResidualScalingContract)
//!
//! Each row of the residual and Jacobian is divided by a characteristic scale factor:
//!
//! ```text
//! F_scaled[i] = F[i] / scale[i]
//! J_scaled[i, :] = J[i, :] / scale[i]
//! ```
//!
//! This prevents species with large mole numbers from dominating the residual norm.
//!
//! ### Variable Scaling (VariableScalingContract) — UNUSED
//!
//! Declared but not currently used in the solve pipeline. Would scale the
//! log-mole variables to improve conditioning of the normal equations.
//!
//! ## Solution
//!
//! [`EquilibriumSolution`] bundles the solved log-moles, physical moles,
//! thermodynamic conditions, and an optional validation report. It is the
//! immutable output of a successful solve.
//!
//! # Key Structures
//!
//! | Structure | Role |
//! |-----------|------|
//! | [`EquilibriumProblem`] | Validated, self-contained problem description |
//! | [`PreparedEquilibriumProblem`] | Ready-to-solve problem with precomputed data |
//! | [`EquilibriumSolution`] | Immutable result of a successful solve |
//! | [`EquilibriumConditions`] | Thermodynamic conditions (T, P, P0) |
//! | [`LogMolesInitialGuess`] | Typed initial guess with trace seeding |
//! | [`ResidualScalingContract`] | Row scaling for residual/Jacobian |
//! | [`VariableScalingContract`] | Variable scaling (currently unused) |
//! | [`FormulationDiagnostics`] | Condition number and capacity analysis |
//! | [`SpeciesCapacityReport`] | Per-species limiting factor analysis |
//! | [`EquilibriumProblemPreview`] | Human-readable problem summary |
//! | [`TraceSpeciesSeedPolicy`] | Policy for seeding trace species |
//!
//! # Dataflow
//!
//! ```text
//!   User input (components, moles, conditions)
//!     │
//!     v
//!   EquilibriumProblem::new()
//!     ├── Validates all dimensions and physical invariants
//!     ├── Stores initial moles, element composition, Gibbs functions
//!     └── Stores phase assignments and activity models
//!     │
//!     v
//!   PreparedEquilibriumProblem::new(problem)
//!     ├── Computes SVD reaction basis
//!     ├── Computes row scaling factors
//!     ├── Builds residual/Jacobian closures
//!     └── Returns PreparedEquilibriumProblem
//!     │
//!     v
//!   EquilibriumLogMoles::from_problem(problem)
//!     ├── Copies data into mutable solver state
//!     └── Ready for solve()
//!     │
//!     v
//!   EquilibriumSolution (after successful solve)
//!     ├── log_moles(), moles() — solution vectors
//!     ├── conditions() — T, P, P0
//!     └── validation() — optional K_eq cross-validation
//! ```
//!
//! # Examples
//!
//! ```rust, ignore
//! use KiThe::Thermodynamics::ChemEquilibrium::equilibrium_problem::*;
//!
//! let problem = EquilibriumProblem::new(
//!     components,           // Vec<EquilibriumComponentDescriptor>
//!     initial_moles,        // Vec<f64>
//!     None,                 // optional log_moles
//!     element_composition,  // DMatrix<f64>
//!     gibbs_functions,      // Vec<GibbsFn>
//!     phases,               // Vec<Phase>
//!     conditions,           // EquilibriumConditions
//! ).unwrap();
//!
//! let prepared = PreparedEquilibriumProblem::new(problem).unwrap();
//! let preview = prepared.preview().unwrap();
//! println!("{}", preview);
//! ```
//!
//! # Non-obvious Details
//!
//! - `EquilibriumProblem::validate()` checks **all** invariants: dimensions,
//!   finite values, positive pressure, non-negative initial moles, phase
//!   assignments, and element composition consistency.
//! - `PreparedEquilibriumProblem` precomputes the SVD reaction basis and scaling
//!   factors, making it efficient to evaluate residuals at multiple points.
//! - `VariableScalingContract` is **declared but unused** in the current pipeline
//!   (see SourceCraft Diagnostics item B.2 in TODO_ANALYSIS.md).
//! - `DEFAULT_TRACE_MOLE_FLOOR = 1e-30` is the default seed for species with
//!   zero initial moles. This is a numerical floor, not a physical cutoff.
//!
//! # Related Modules
//!
//! - [`equilibrium_log_moles`](super::equilibrium_log_moles) — mutable solver that consumes problems
//! - [`equilibrium_component`](super::equilibrium_component) — component descriptor types
//! - [`equilibrium_ids`](super::equilibrium_ids) — typed index wrappers
//! - [`equilibrium_validation`](super::equilibrium_validation) — acceptance gate for solutions
//!

use crate::Thermodynamics::ChemEquilibrium::equilibrium_component::EquilibriumComponentDescriptor;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_ids::{
    ElementId, PhaseIndex, ReactionId, SpeciesId,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_log_moles::{
    GibbsFn, Phase, compute_species_moles, equilibrium_scaling,
    evaluate_equilibrium_logmole_jacobian, evaluate_equilibrium_logmole_residual,
    reaction_phase_stoichiometry, scale_jacobian_rows, scale_residual_rows, species_to_phase_map,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_nonlinear::{
    ReactionBasis, ReactionExtentError, compute_reaction_basis,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_validation::EquilibriumCandidateReport;
use nalgebra::{DMatrix, linalg::SVD};
use std::collections::HashSet;
use std::fmt;

/// Default positive seed for a species absent from the physical initial state.
///
/// This is a numerical coordinate floor, not a physical concentration cutoff.
/// Phase-activation thresholds remain independent settings.
pub const DEFAULT_TRACE_MOLE_FLOOR: f64 = 1e-30;

/// Explicit policy for converting physical mole numbers into log-mole seeds.
///
/// The solver still uses an absolute default floor for backward compatibility,
/// but callers that want a scale-aware seed can now choose that policy
/// explicitly instead of relying on hidden clipping rules.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum TraceSpeciesSeedPolicy {
    /// Use a fixed positive floor for every zero or trace species.
    Absolute { floor: f64 },
    /// Scale the floor to the largest provided mole number, with a minimum cap.
    RelativeToLargestInitialMole {
        /// Fraction of the largest positive mole number used as a floor.
        fraction: f64,
        /// Minimum floor used when the system is entirely zero or the scaled
        /// floor would otherwise become too small.
        minimum_floor: f64,
    },
}

impl TraceSpeciesSeedPolicy {
    /// Returns a validated floor for the supplied initial mole vector.
    pub fn trace_floor_for(self, moles: &[f64]) -> Result<f64, ReactionExtentError> {
        match self {
            Self::Absolute { floor } => {
                if !floor.is_finite() || floor <= 0.0 {
                    return Err(ReactionExtentError::InvalidProblem {
                        field: "trace_floor",
                        message: "trace floor must be finite and strictly positive".to_string(),
                    });
                }
                Ok(floor)
            }
            Self::RelativeToLargestInitialMole {
                fraction,
                minimum_floor,
            } => {
                if !fraction.is_finite() || fraction <= 0.0 || fraction > 1.0 {
                    return Err(ReactionExtentError::InvalidProblem {
                        field: "trace_floor_fraction",
                        message: "fraction must be finite and lie in the interval (0, 1]"
                            .to_string(),
                    });
                }
                if !minimum_floor.is_finite() || minimum_floor <= 0.0 {
                    return Err(ReactionExtentError::InvalidProblem {
                        field: "trace_floor_minimum",
                        message: "minimum floor must be finite and strictly positive".to_string(),
                    });
                }

                let mut largest_initial_mole = 0.0;
                for (index, &value) in moles.iter().enumerate() {
                    if !value.is_finite() || value < 0.0 {
                        return Err(ReactionExtentError::InvalidProblem {
                            field: "initial_moles",
                            message: format!("entry {index} must be finite and non-negative"),
                        });
                    }
                    if value > largest_initial_mole {
                        largest_initial_mole = value;
                    }
                }

                if largest_initial_mole <= 0.0 {
                    return Ok(minimum_floor);
                }

                Ok((largest_initial_mole * fraction).max(minimum_floor))
            }
        }
    }
}

/// Thermodynamic conditions shared by every equilibrium calculation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EquilibriumConditions {
    /// Temperature in Kelvin (K).
    temperature: f64,
    /// System pressure in Pascals (Pa).
    pressure: f64,
    /// Standard-state reference pressure in Pascals (Pa), typically 101325 Pa = 1 atm.
    reference_pressure: f64,
}

impl EquilibriumConditions {
    /// Builds finite, strictly positive thermodynamic conditions.
    pub fn new(
        temperature: f64,
        pressure: f64,
        reference_pressure: f64,
    ) -> Result<Self, ReactionExtentError> {
        for (parameter, value) in [
            ("temperature", temperature),
            ("pressure", pressure),
            ("reference_pressure", reference_pressure),
        ] {
            if !value.is_finite() || value <= 0.0 {
                return Err(ReactionExtentError::InvalidConditions { parameter, value });
            }
        }

        Ok(Self {
            temperature,
            pressure,
            reference_pressure,
        })
    }

    /// Temperature in K.
    pub fn temperature(self) -> f64 {
        self.temperature
    }

    /// System pressure in Pa.
    pub fn pressure(self) -> f64 {
        self.pressure
    }

    /// Standard-state reference pressure in Pa.
    pub fn reference_pressure(self) -> f64 {
        self.reference_pressure
    }
}

/// A solver iterate expressed as `ln(n_i)`, not as physical mole numbers.
///
/// Keeping this distinction typed prevents the historical mistake of passing
/// ordinary mole numbers to a residual that expects log-moles.
#[derive(Debug, Clone, PartialEq)]
pub struct LogMolesInitialGuess(Vec<f64>);

impl LogMolesInitialGuess {
    /// Builds a seed using the canonical trace-species coordinate floor.
    pub fn from_initial_moles(moles: &[f64]) -> Result<Self, ReactionExtentError> {
        Self::from_moles_with_policy(
            moles,
            TraceSpeciesSeedPolicy::Absolute {
                floor: DEFAULT_TRACE_MOLE_FLOOR,
            },
        )
    }

    /// Builds a seed from a typed trace-species policy.
    pub fn from_moles_with_policy(
        moles: &[f64],
        policy: TraceSpeciesSeedPolicy,
    ) -> Result<Self, ReactionExtentError> {
        let trace_floor = policy.trace_floor_for(moles)?;
        Self::from_moles(moles, trace_floor)
    }

    /// Validates a caller-provided log-mole iterate.
    pub fn new(log_moles: Vec<f64>) -> Result<Self, ReactionExtentError> {
        if log_moles.iter().any(|value| !value.is_finite()) {
            return Err(ReactionExtentError::InvalidProblem {
                field: "initial_log_moles",
                message: "every log-mole value must be finite".to_string(),
            });
        }
        Ok(Self(log_moles))
    }

    /// Converts physical mole numbers into a finite log-mole iterate.
    ///
    /// Zero species are represented by the supplied positive trace floor;
    /// negative or non-finite mole numbers are rejected instead of clipped.
    pub fn from_moles(moles: &[f64], trace_floor: f64) -> Result<Self, ReactionExtentError> {
        if !trace_floor.is_finite() || trace_floor <= 0.0 {
            return Err(ReactionExtentError::InvalidProblem {
                field: "trace_floor",
                message: "trace floor must be finite and strictly positive".to_string(),
            });
        }

        let mut log_moles = Vec::with_capacity(moles.len());
        for (index, &moles_i) in moles.iter().enumerate() {
            if !moles_i.is_finite() || moles_i < 0.0 {
                return Err(ReactionExtentError::InvalidProblem {
                    field: "initial_moles",
                    message: format!("entry {index} must be finite and non-negative"),
                });
            }
            log_moles.push(moles_i.max(trace_floor).ln());
        }
        Self::new(log_moles)
    }

    /// Borrow the log-mole values in solver order.
    pub fn as_slice(&self) -> &[f64] {
        &self.0
    }

    pub(crate) fn into_inner(self) -> Vec<f64> {
        self.0
    }
}

/// Fully validated input for one canonical equilibrium solve.
pub struct EquilibriumProblem {
    /// Canonical phase-qualified identity in exact solver order.
    components: Vec<EquilibriumComponentDescriptor>,
    /// Derived display labels retained for legacy numerical/reporting callers.
    labels: Vec<String>,
    /// Initial physical mole numbers `n_i^0` for each species.
    initial_moles: Vec<f64>,
    /// Initial guess in log-mole space `y_i = ln(n_i)`, or a trace-floor seed.
    initial_log_moles: LogMolesInitialGuess,
    /// Element composition matrix `A` of shape `(species × elements)`.
    /// `A[i, j]` = number of atoms of element `j` in species `i`.
    element_composition: DMatrix<f64>,
    /// Standard Gibbs free energy functions `g_i(T)` for each species, in J/mol.
    gibbs: Vec<GibbsFn>,
    /// Phase descriptors: activity model and species indices for each phase.
    phases: Vec<Phase>,
    /// Thermodynamic conditions (T, P, P0) for this problem.
    conditions: EquilibriumConditions,
}

/// Immutable numerical form derived from one validated [`EquilibriumProblem`].
///
/// The prepared problem owns every matrix and ordering decision needed by a
/// nonlinear backend. It is deliberately free of solver progress, caches, and
/// reports, so multiple backends can evaluate exactly the same formulation.
pub struct PreparedEquilibriumProblem {
    /// The original validated problem that was prepared.
    problem: EquilibriumProblem,
    /// SVD-derived reaction basis (stoichiometric nullspace matrix `ν`).
    reaction_basis: ReactionBasis,
    /// Conserved element totals `b_0 = A^T · n0` — total moles of each element.
    element_totals: Vec<f64>,
    /// Maps each species index to its phase index: `species_phase[i]` = phase of species i.
    species_phase: Vec<usize>,
    /// Aggregated stoichiometry per phase: `phase_stoichiometry[p]` = sum of ν over species in phase p.
    phase_stoichiometry: Vec<Vec<f64>>,
}

/// Validated row-scaling metadata for one prepared equilibrium formulation.
///
/// This is intentionally separate from the raw `Vec<f64>` so callers do not
/// silently pass arbitrary arrays around as if they were a stable contract.
#[derive(Debug, Clone, PartialEq)]
pub struct ResidualScalingContract {
    /// Row-scale factors: `scale[i]` multiplies residual/Jacobian row `i`.
    /// Each factor must be finite and strictly positive.
    scale: Vec<f64>,
}

impl ResidualScalingContract {
    /// Builds a validated row-scaling contract from explicit scale factors.
    pub fn new(scale: Vec<f64>) -> Result<Self, ReactionExtentError> {
        if scale.is_empty() {
            return Err(ReactionExtentError::InvalidProblem {
                field: "residual_scale",
                message: "residual scale must contain at least one entry".to_string(),
            });
        }
        for (index, value) in scale.iter().enumerate() {
            if !value.is_finite() || *value <= 0.0 {
                return Err(ReactionExtentError::InvalidProblem {
                    field: "residual_scale",
                    message: format!("scale[{index}] must be finite and positive"),
                });
            }
        }
        Ok(Self { scale })
    }

    /// Borrow the validated scale vector.
    pub fn as_slice(&self) -> &[f64] {
        &self.scale
    }

    /// Applies the scale to one residual vector.
    pub fn apply_residual(&self, residual: Vec<f64>) -> Result<Vec<f64>, ReactionExtentError> {
        scale_residual_rows(residual, &self.scale)
    }

    /// Applies the same scale to one Jacobian matrix.
    pub fn apply_jacobian(
        &self,
        jacobian: DMatrix<f64>,
    ) -> Result<DMatrix<f64>, ReactionExtentError> {
        scale_jacobian_rows(jacobian, &self.scale)
    }
}

/// Validated coordinate scaling for the nonlinear iterate itself.
///
/// This is intentionally separate from [`ResidualScalingContract`]. Row
/// scaling changes the residual/Jacobian rows that a backend sees, while
/// variable scaling changes the coordinate system in which the solver walks.
/// Keeping these two concepts in different types prevents one from being
/// mistaken for the other.
#[derive(Debug, Clone, PartialEq)]
pub struct VariableScalingContract {
    /// Variable-scale factors: `scale[i]` scales coordinate `i` in the solver iterate.
    /// Each factor must be finite and strictly positive.
    scale: Vec<f64>,
}

impl VariableScalingContract {
    /// Builds a validated coordinate-scale contract.
    pub fn new(scale: Vec<f64>) -> Result<Self, ReactionExtentError> {
        if scale.is_empty() {
            return Err(ReactionExtentError::InvalidProblem {
                field: "variable_scale",
                message: "variable scale must contain at least one entry".to_string(),
            });
        }
        for (index, value) in scale.iter().enumerate() {
            if !value.is_finite() || *value <= 0.0 {
                return Err(ReactionExtentError::InvalidProblem {
                    field: "variable_scale",
                    message: format!("scale[{index}] must be finite and positive"),
                });
            }
        }
        Ok(Self { scale })
    }

    /// Borrow the validated coordinate scale.
    pub fn as_slice(&self) -> &[f64] {
        &self.scale
    }

    /// Applies the coordinate scale to one solver iterate.
    ///
    /// This is a standalone transform: it does not alter residual rows or
    /// Jacobian rows.
    pub fn apply_iterate(&self, iterate: &[f64]) -> Result<Vec<f64>, ReactionExtentError> {
        if iterate.len() != self.scale.len() {
            return Err(ReactionExtentError::DimensionMismatch(format!(
                "iterate has {} entries but variable scale has {} entries",
                iterate.len(),
                self.scale.len(),
            )));
        }

        let mut scaled = Vec::with_capacity(iterate.len());
        for (index, (value, factor)) in iterate.iter().zip(self.scale.iter()).enumerate() {
            if !value.is_finite() {
                return Err(ReactionExtentError::InvalidProblem {
                    field: "variable_iterate",
                    message: format!("iterate[{index}] must be finite"),
                });
            }
            scaled.push(value / factor);
        }
        Ok(scaled)
    }

    /// Restores a scaled iterate back to physical coordinates.
    pub fn unscale_iterate(&self, scaled_iterate: &[f64]) -> Result<Vec<f64>, ReactionExtentError> {
        if scaled_iterate.len() != self.scale.len() {
            return Err(ReactionExtentError::DimensionMismatch(format!(
                "scaled iterate has {} entries but variable scale has {} entries",
                scaled_iterate.len(),
                self.scale.len(),
            )));
        }

        let mut iterate = Vec::with_capacity(scaled_iterate.len());
        for (index, (value, factor)) in scaled_iterate.iter().zip(self.scale.iter()).enumerate() {
            if !value.is_finite() {
                return Err(ReactionExtentError::InvalidProblem {
                    field: "variable_iterate",
                    message: format!("scaled iterate[{index}] must be finite"),
                });
            }
            iterate.push(value * factor);
        }
        Ok(iterate)
    }
}

/// Immutable accepted equilibrium state.
///
/// This is a read-only boundary between a successful solve and its callers.
/// It couples log-moles, physical moles, thermodynamic conditions, and the
/// evidence that allowed the candidate to be accepted.
#[derive(Debug, Clone, PartialEq)]
pub struct EquilibriumSolution {
    /// Accepted solution in log-mole space: `y_i = ln(n_i)`.
    log_moles: Vec<f64>,
    /// Reconstructed physical moles: `n_i = exp(y_i)`.
    moles: Vec<f64>,
    /// Thermodynamic conditions under which this solution was obtained.
    conditions: EquilibriumConditions,
    /// Backend-independent acceptance evidence for this solution.
    validation: EquilibriumCandidateReport,
}

impl EquilibriumSolution {
    pub(crate) fn new(
        log_moles: Vec<f64>,
        moles: Vec<f64>,
        conditions: EquilibriumConditions,
        validation: EquilibriumCandidateReport,
    ) -> Result<Self, ReactionExtentError> {
        if log_moles.len() != moles.len() {
            return Err(ReactionExtentError::DimensionMismatch(format!(
                "solution has {} log-moles but {} physical mole values",
                log_moles.len(),
                moles.len()
            )));
        }
        if log_moles.iter().any(|value| !value.is_finite())
            || moles
                .iter()
                .any(|value| !value.is_finite() || *value <= 0.0)
        {
            return Err(ReactionExtentError::InvalidCandidate {
                field: "accepted_solution",
                message: "accepted solution contains non-finite or non-positive values".to_string(),
            });
        }

        Ok(Self {
            log_moles,
            moles,
            conditions,
            validation,
        })
    }

    /// Log-mole nonlinear iterate in deterministic species order.
    pub fn log_moles(&self) -> &[f64] {
        &self.log_moles
    }

    /// Physical species mole numbers in the same deterministic order.
    pub fn moles(&self) -> &[f64] {
        &self.moles
    }

    /// Conditions used during this accepted solve.
    pub fn conditions(&self) -> EquilibriumConditions {
        self.conditions
    }

    /// Backend-independent candidate-validation evidence.
    pub fn validation(&self) -> &EquilibriumCandidateReport {
        &self.validation
    }
}

/// Capacity evidence for one element while forming a single species.
///
/// If a species consumes `a_ie` atoms of an element and the system currently
/// contains `b_e` atoms of that element, the implied species upper bound from
/// that element is `b_e / a_ie`.
#[derive(Debug, Clone, PartialEq)]
pub struct SpeciesCapacityLimit {
    /// Species in the prepared ordering.
    pub species_id: SpeciesId,
    /// Element in the prepared ordering.
    pub element_id: ElementId,
    /// Zero-based element index.
    pub element_index: usize,
    /// Total amount of the limiting element currently available.
    pub element_total_moles: f64,
    /// Stoichiometric coefficient `a_ie` from the element matrix.
    pub stoichiometric_coefficient: f64,
    /// Species upper bound implied by this element.
    pub implied_species_capacity: f64,
}

/// Exact feasibility capacity for a single species.
///
/// This report is intentionally evidence, not solver policy. It can be used
/// for diagnostics, impossible-composition checks, and test fixtures without
/// changing how the nonlinear backend searches.
#[derive(Debug, Clone, PartialEq)]
pub struct SpeciesCapacityReport {
    /// Species in the prepared ordering.
    pub species_id: SpeciesId,
    /// Zero-based species index.
    pub species_index: usize,
    /// Species name in deterministic order.
    pub species_name: String,
    /// Element-by-element capacity evidence, in matrix-column order.
    pub limits: Vec<SpeciesCapacityLimit>,
    /// The smallest capacity evidence among `limits`.
    pub limiting_limit: SpeciesCapacityLimit,
}

impl SpeciesCapacityReport {
    /// Exact upper bound implied by the limiting elemental balance.
    pub fn maximum_moles(&self) -> f64 {
        self.limiting_limit.implied_species_capacity
    }
}

impl fmt::Display for SpeciesCapacityReport {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}: max {:.6} mol limited by element {}",
            self.species_name,
            self.maximum_moles(),
            self.limiting_limit.element_index
        )
    }
}

/// One nearly-null singular direction of the prepared formulation.
#[derive(Debug, Clone, PartialEq)]
pub struct NearNullDirection {
    /// Singular value associated with the direction.
    pub singular_value: f64,
    /// Right-singular vector in species order.
    pub direction: Vec<f64>,
}

/// Optional SVD-based formulation diagnostics for a prepared problem.
///
/// This report is built only when explicitly requested. It is evidence for
/// tests, CLI output, and developer inspection; it does not alter solver
/// policy.
#[derive(Debug, Clone, PartialEq)]
pub struct FormulationDiagnostics {
    /// Numerical rank of the element composition matrix for the requested tolerance.
    pub element_rank: usize,
    /// Number of independent reactions implied by that rank.
    pub reaction_count: usize,
    /// Singular values of `A^T`, in descending order.
    pub singular_values: Vec<f64>,
    /// Largest singular value divided by the smallest singular value above the tolerance.
    pub condition_estimate: Option<f64>,
    /// Tolerance used when classifying nearly-null singular directions.
    pub tolerance: f64,
    /// Singular directions at or below the tolerance.
    pub near_null_directions: Vec<NearNullDirection>,
}

impl fmt::Display for FormulationDiagnostics {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(
            f,
            "element rank = {}, reaction count = {}",
            self.element_rank, self.reaction_count
        )?;
        match self.condition_estimate {
            Some(value) => writeln!(f, "condition estimate = {value:.6e}")?,
            None => writeln!(f, "condition estimate = n/a")?,
        }
        for (index, singular_value) in self.singular_values.iter().enumerate() {
            writeln!(f, "sv[{index}] = {singular_value:.6e}")?;
        }
        for (index, direction) in self.near_null_directions.iter().enumerate() {
            writeln!(
                f,
                "null[{index}] sv={:.6e} dim={}",
                direction.singular_value,
                direction.direction.len()
            )?;
        }
        Ok(())
    }
}

/// Typed summary for CLI, GUI, examples, and regression fixtures.
///
/// This is a read-only snapshot of the validated problem plus the exact
/// capacity evidence derived from it. It does not own any solver state or
/// policy.
#[derive(Debug, Clone, PartialEq)]
pub struct EquilibriumProblemPreview {
    /// Thermodynamic conditions used for the problem.
    pub conditions: EquilibriumConditions,
    /// Species in the deterministic solver ordering.
    pub species: Vec<String>,
    /// Initial physical moles in the same ordering.
    pub initial_moles: Vec<f64>,
    /// Element composition matrix copied from the prepared problem.
    pub element_composition: DMatrix<f64>,
    /// Exact conserved element totals.
    pub element_totals: Vec<f64>,
    /// Rank of the discovered reaction basis.
    pub reaction_basis_rank: usize,
    /// Number of independent reactions in the prepared formulation.
    pub reaction_count: usize,
    /// Species-feasibility evidence in species order.
    pub species_capacity_reports: Vec<SpeciesCapacityReport>,
    /// Optional SVD-based diagnostics requested explicitly by the caller.
    pub formulation_diagnostics: Option<FormulationDiagnostics>,
}

/// One human-readable summary row for a preview report.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct EquilibriumPreviewRow {
    /// Logical section name, e.g. `problem` or `species_capacity`.
    pub section: &'static str,
    /// Stable row label.
    pub label: String,
    /// Human-readable value.
    pub value: String,
}

impl fmt::Display for EquilibriumPreviewRow {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "[{}] {} = {}", self.section, self.label, self.value)
    }
}

impl EquilibriumProblemPreview {
    /// Returns stable summary rows for CLI output and snapshot tests.
    pub fn summary_rows(&self) -> Vec<EquilibriumPreviewRow> {
        let mut rows = vec![
            EquilibriumPreviewRow {
                section: "problem",
                label: "temperature".to_string(),
                value: format!("{:.6}", self.conditions.temperature()),
            },
            EquilibriumPreviewRow {
                section: "problem",
                label: "pressure".to_string(),
                value: format!("{:.6}", self.conditions.pressure()),
            },
            EquilibriumPreviewRow {
                section: "problem",
                label: "species_count".to_string(),
                value: self.species.len().to_string(),
            },
            EquilibriumPreviewRow {
                section: "problem",
                label: "element_count".to_string(),
                value: self.element_composition.ncols().to_string(),
            },
            EquilibriumPreviewRow {
                section: "problem",
                label: "reaction_count".to_string(),
                value: self.reaction_count.to_string(),
            },
        ];

        for report in &self.species_capacity_reports {
            rows.push(EquilibriumPreviewRow {
                section: "species_capacity",
                label: report.species_name.clone(),
                value: report.to_string(),
            });
        }

        if let Some(diagnostics) = &self.formulation_diagnostics {
            rows.push(EquilibriumPreviewRow {
                section: "diagnostics",
                label: "element_rank".to_string(),
                value: diagnostics.to_string(),
            });
        }

        rows
    }
}

impl fmt::Display for EquilibriumProblemPreview {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for row in self.summary_rows() {
            writeln!(f, "{row}")?;
        }
        Ok(())
    }
}

impl PreparedEquilibriumProblem {
    /// Validates and derives the deterministic matrices used by the log-moles formulation.
    pub fn new(problem: EquilibriumProblem) -> Result<Self, ReactionExtentError> {
        problem.validate()?;

        let reaction_basis = compute_reaction_basis(problem.element_composition(), 1e-6)?;
        let species_count = problem.species().len();
        let equation_count =
            reaction_basis.reactions.ncols() + problem.element_composition().ncols();
        if equation_count != species_count {
            return Err(ReactionExtentError::DimensionMismatch(format!(
                "prepared log-moles system has {equation_count} equations for {species_count} species"
            )));
        }

        let species_phase = species_to_phase_map(problem.phases(), species_count)?;
        let phase_stoichiometry =
            reaction_phase_stoichiometry(&reaction_basis.reactions, problem.phases());
        let element_totals = (0..problem.element_composition().ncols())
            .map(|element| {
                problem
                    .initial_moles()
                    .iter()
                    .enumerate()
                    .map(|(species, moles)| {
                        problem.element_composition()[(species, element)] * moles
                    })
                    .sum()
            })
            .collect();

        Ok(Self {
            problem,
            reaction_basis,
            element_totals,
            species_phase,
            phase_stoichiometry,
        })
    }

    /// Evaluates the canonical residual without changing solver state.
    pub fn residual(&self, log_moles: &[f64]) -> Result<Vec<f64>, ReactionExtentError> {
        let conditions = self.problem.conditions();
        evaluate_equilibrium_logmole_residual(
            log_moles,
            &self.reaction_basis.reactions,
            self.problem.element_composition(),
            &self.element_totals,
            self.problem.gibbs(),
            self.problem.phases(),
            conditions.temperature(),
            conditions.pressure(),
            conditions.reference_pressure(),
            &self.species_phase,
            &self.phase_stoichiometry,
        )
    }

    /// Evaluates the analytical Jacobian of the same canonical residual.
    pub fn jacobian(&self, log_moles: &[f64]) -> Result<DMatrix<f64>, ReactionExtentError> {
        evaluate_equilibrium_logmole_jacobian(
            log_moles,
            &self.reaction_basis.reactions,
            self.problem.element_composition(),
            &self.species_phase,
            &self.phase_stoichiometry,
            self.problem.phases().len(),
        )
    }

    /// Reconstructs physical mole numbers in the prepared species ordering.
    ///
    /// This keeps coordinate conversion at the immutable formulation boundary
    /// and rejects an iterate of the wrong size before it reaches a result.
    pub fn reconstruct_moles(&self, log_moles: &[f64]) -> Result<Vec<f64>, ReactionExtentError> {
        if log_moles.len() != self.problem.species().len() {
            return Err(ReactionExtentError::DimensionMismatch(format!(
                "log-mole vector has {} entries for {} prepared species",
                log_moles.len(),
                self.problem.species().len(),
            )));
        }
        compute_species_moles(log_moles)
    }

    /// Packages an accepted iterate into the immutable public result.
    ///
    /// Candidate acceptance remains the responsibility of the common
    /// validation gate; this method prevents result construction from
    /// bypassing the prepared ordering and thermodynamic conditions.
    pub fn accepted_solution(
        &self,
        log_moles: Vec<f64>,
        validation: EquilibriumCandidateReport,
    ) -> Result<EquilibriumSolution, ReactionExtentError> {
        let moles = self.reconstruct_moles(&log_moles)?;
        EquilibriumSolution::new(log_moles, moles, self.problem.conditions(), validation)
    }

    /// Computes deterministic row scales for this prepared formulation.
    pub fn residual_scale(&self) -> Result<Vec<f64>, ReactionExtentError> {
        self.residual_scaling_contract()
            .map(|contract| contract.scale)
    }

    /// Builds the typed row-scaling contract for this formulation.
    pub fn residual_scaling_contract(
        &self,
    ) -> Result<ResidualScalingContract, ReactionExtentError> {
        ResidualScalingContract::new(equilibrium_scaling(
            &self.reaction_basis.reactions,
            self.problem.element_composition(),
            self.problem.gibbs(),
            &self.element_totals,
            self.problem.conditions().temperature(),
        )?)
    }

    /// Evaluates the canonical residual and applies an explicit validated row scale.
    pub fn scaled_residual(
        &self,
        log_moles: &[f64],
        scale: &[f64],
    ) -> Result<Vec<f64>, ReactionExtentError> {
        ResidualScalingContract::new(scale.to_vec())?.apply_residual(self.residual(log_moles)?)
    }

    /// Evaluates the canonical Jacobian and applies the matching row scale.
    pub fn scaled_jacobian(
        &self,
        log_moles: &[f64],
        scale: &[f64],
    ) -> Result<DMatrix<f64>, ReactionExtentError> {
        ResidualScalingContract::new(scale.to_vec())?.apply_jacobian(self.jacobian(log_moles)?)
    }

    /// Computes exact elemental feasibility evidence for one species.
    pub fn species_capacity_report(
        &self,
        species_index: usize,
    ) -> Result<SpeciesCapacityReport, ReactionExtentError> {
        let species_id = self.species_id(species_index)?;
        let mut limits = Vec::new();
        let mut limiting_limit = None::<SpeciesCapacityLimit>;

        for element_index in 0..self.problem.element_composition().ncols() {
            let stoichiometric_coefficient =
                self.problem.element_composition()[(species_index, element_index)];
            if !stoichiometric_coefficient.is_finite() {
                return Err(ReactionExtentError::InvalidProblem {
                    field: "element_composition",
                    message: format!(
                        "species row {species_index} contains a non-finite coefficient"
                    ),
                });
            }
            if stoichiometric_coefficient < 0.0 {
                return Err(ReactionExtentError::InvalidProblem {
                    field: "element_composition",
                    message: format!(
                        "species row {species_index} contains a negative coefficient at element column {element_index}"
                    ),
                });
            }
            if stoichiometric_coefficient == 0.0 {
                continue;
            }

            let element_total_moles = self.element_totals[element_index];
            if !element_total_moles.is_finite() {
                return Err(ReactionExtentError::InvalidProblem {
                    field: "element_totals",
                    message: format!("element total {element_index} is not finite"),
                });
            }
            let implied_species_capacity = element_total_moles / stoichiometric_coefficient;
            if !implied_species_capacity.is_finite() {
                return Err(ReactionExtentError::InvalidProblem {
                    field: "species_capacity",
                    message: format!(
                        "species {species_index} capacity became non-finite for element {element_index}"
                    ),
                });
            }

            let limit = SpeciesCapacityLimit {
                species_id,
                element_id: self.element_id(element_index)?,
                element_index,
                element_total_moles,
                stoichiometric_coefficient,
                implied_species_capacity,
            };

            match &limiting_limit {
                None => limiting_limit = Some(limit.clone()),
                Some(current)
                    if limit.implied_species_capacity < current.implied_species_capacity =>
                {
                    limiting_limit = Some(limit.clone());
                }
                _ => {}
            }
            limits.push(limit);
        }

        let limiting_limit = limiting_limit.ok_or_else(|| ReactionExtentError::InvalidProblem {
            field: "element_composition",
            message: format!(
                "species row {species_index} contains no positive elemental coefficients"
            ),
        })?;

        Ok(SpeciesCapacityReport {
            species_id,
            species_index,
            species_name: self.problem.species()[species_index].clone(),
            limits,
            limiting_limit,
        })
    }

    /// Computes exact capacity evidence for every species in solver order.
    pub fn species_capacity_reports(
        &self,
    ) -> Result<Vec<SpeciesCapacityReport>, ReactionExtentError> {
        (0..self.problem.species().len())
            .map(|index| self.species_capacity_report(index))
            .collect()
    }

    /// Builds a typed, read-only summary of the validated equilibrium problem.
    pub fn preview(&self) -> Result<EquilibriumProblemPreview, ReactionExtentError> {
        Ok(EquilibriumProblemPreview {
            conditions: self.problem.conditions(),
            species: self.problem.species().to_vec(),
            initial_moles: self.problem.initial_moles().to_vec(),
            element_composition: self.problem.element_composition().clone(),
            element_totals: self.element_totals.clone(),
            reaction_basis_rank: self.reaction_basis.rank,
            reaction_count: self.reaction_basis.num_reactions,
            species_capacity_reports: self.species_capacity_reports()?,
            formulation_diagnostics: None,
        })
    }

    /// Builds the preview plus optional SVD diagnostics for the requested tolerance.
    pub fn preview_with_diagnostics(
        &self,
        tolerance: f64,
    ) -> Result<EquilibriumProblemPreview, ReactionExtentError> {
        let mut preview = self.preview()?;
        preview.formulation_diagnostics = Some(self.formulation_diagnostics(tolerance)?);
        Ok(preview)
    }

    /// Computes optional SVD-based formulation diagnostics on demand.
    pub fn formulation_diagnostics(
        &self,
        tolerance: f64,
    ) -> Result<FormulationDiagnostics, ReactionExtentError> {
        if !tolerance.is_finite() || tolerance <= 0.0 {
            return Err(ReactionExtentError::InvalidProblem {
                field: "diagnostic_tolerance",
                message: "diagnostic tolerance must be finite and strictly positive".to_string(),
            });
        }

        let svd = SVD::new(self.problem.element_composition().transpose(), true, true);
        let singular_values: Vec<f64> = svd.singular_values.iter().copied().collect();
        let numerical_rank = singular_values
            .iter()
            .filter(|&&value| value > tolerance)
            .count();

        let max_sv = singular_values.iter().copied().fold(0.0_f64, f64::max);
        let min_positive_sv = singular_values
            .iter()
            .copied()
            .filter(|value| *value > tolerance)
            .reduce(f64::min);
        let condition_estimate = match (max_sv.is_finite(), min_positive_sv) {
            (true, Some(min_sv)) if min_sv > 0.0 => Some(max_sv / min_sv),
            _ => None,
        };

        let mut near_null_directions = Vec::new();
        for reaction in 0..self.reaction_basis.reactions.ncols() {
            near_null_directions.push(NearNullDirection {
                singular_value: 0.0,
                direction: self
                    .reaction_basis
                    .reactions
                    .column(reaction)
                    .iter()
                    .copied()
                    .collect(),
            });
        }

        Ok(FormulationDiagnostics {
            element_rank: numerical_rank,
            reaction_count: self.problem.species().len().saturating_sub(numerical_rank),
            singular_values,
            condition_estimate,
            tolerance,
            near_null_directions,
        })
    }

    /// Borrow the validated domain input used to create this numerical snapshot.
    pub fn problem(&self) -> &EquilibriumProblem {
        &self.problem
    }

    /// Reaction basis in the same species ordering as [`Self::problem`].
    pub fn reaction_basis(&self) -> &ReactionBasis {
        &self.reaction_basis
    }

    /// Conserved elemental totals in the element-matrix column order.
    pub fn element_totals(&self) -> &[f64] {
        &self.element_totals
    }

    /// Phase index for each species in deterministic species order.
    pub fn species_phase(&self) -> &[usize] {
        &self.species_phase
    }

    /// Per-reaction change in total moles for each phase.
    pub fn phase_stoichiometry(&self) -> &[Vec<f64>] {
        &self.phase_stoichiometry
    }

    /// Returns a typed species id in the prepared ordering.
    pub fn species_id(&self, index: usize) -> Result<SpeciesId, ReactionExtentError> {
        SpeciesId::new(index, self.problem.species().len())
    }

    /// Returns a typed element id in the prepared ordering.
    pub fn element_id(&self, index: usize) -> Result<ElementId, ReactionExtentError> {
        ElementId::new(index, self.element_totals.len())
    }

    /// Returns a typed phase id in the prepared ordering.
    pub fn phase_index(&self, index: usize) -> Result<PhaseIndex, ReactionExtentError> {
        PhaseIndex::new(index, self.problem.phases().len())
    }

    /// Returns a typed reaction id in the prepared ordering.
    pub fn reaction_id(&self, index: usize) -> Result<ReactionId, ReactionExtentError> {
        ReactionId::new(index, self.reaction_basis.reactions.ncols())
    }

    pub(crate) fn into_legacy_parts(
        self,
    ) -> (
        Vec<String>,
        Vec<f64>,
        Vec<f64>,
        DMatrix<f64>,
        Vec<GibbsFn>,
        Vec<Phase>,
        EquilibriumConditions,
        ReactionBasis,
        Vec<f64>,
        Vec<usize>,
    ) {
        let (
            species,
            initial_moles,
            initial_log_moles,
            element_composition,
            gibbs,
            phases,
            conditions,
        ) = self.problem.into_parts();
        (
            species,
            initial_moles,
            initial_log_moles,
            element_composition,
            gibbs,
            phases,
            conditions,
            self.reaction_basis,
            self.element_totals,
            self.species_phase,
        )
    }
}

impl EquilibriumProblem {
    /// Validates all input data before it is handed to a nonlinear backend.
    #[allow(clippy::too_many_arguments)]
    pub fn new<C>(
        components: Vec<C>,
        initial_moles: Vec<f64>,
        initial_log_moles: LogMolesInitialGuess,
        element_composition: DMatrix<f64>,
        gibbs: Vec<GibbsFn>,
        phases: Vec<Phase>,
        conditions: EquilibriumConditions,
    ) -> Result<Self, ReactionExtentError>
    where
        C: Into<EquilibriumComponentDescriptor>,
    {
        let components = components
            .into_iter()
            .map(Into::into)
            .collect::<Vec<EquilibriumComponentDescriptor>>();
        let labels = components
            .iter()
            .map(EquilibriumComponentDescriptor::label)
            .collect();
        let problem = Self {
            components,
            labels,
            initial_moles,
            initial_log_moles,
            element_composition,
            gibbs,
            phases,
            conditions,
        };
        problem.validate()?;
        Ok(problem)
    }

    /// Validates structural and numerical invariants without mutating the problem.
    pub fn validate(&self) -> Result<(), ReactionExtentError> {
        let species_count = self.components.len();
        if species_count == 0 {
            return Err(ReactionExtentError::InvalidProblem {
                field: "species",
                message: "at least one species is required".to_string(),
            });
        }
        if self
            .components
            .iter()
            .any(|component| component.substance().trim().is_empty())
        {
            return Err(ReactionExtentError::InvalidProblem {
                field: "components",
                message: "component substance names must not be empty".to_string(),
            });
        }
        let unique_components = self
            .components
            .iter()
            .map(EquilibriumComponentDescriptor::id)
            .collect::<HashSet<_>>();
        if unique_components.len() != species_count {
            return Err(ReactionExtentError::InvalidProblem {
                field: "components",
                message: "phase-qualified component identities must be unique".to_string(),
            });
        }
        if self.initial_moles.len() != species_count
            || self.initial_log_moles.as_slice().len() != species_count
            || self.gibbs.len() != species_count
        {
            return Err(ReactionExtentError::DimensionMismatch(format!(
                "species={species_count}, initial_moles={}, initial_log_moles={}, gibbs={}",
                self.initial_moles.len(),
                self.initial_log_moles.as_slice().len(),
                self.gibbs.len(),
            )));
        }
        if self.element_composition.nrows() != species_count
            || self.element_composition.ncols() == 0
        {
            return Err(ReactionExtentError::DimensionMismatch(format!(
                "element composition must have {species_count} rows and at least one element column"
            )));
        }
        if self
            .initial_moles
            .iter()
            .any(|moles| !moles.is_finite() || *moles < 0.0)
        {
            return Err(ReactionExtentError::InvalidProblem {
                field: "initial_moles",
                message: "initial moles must be finite and non-negative".to_string(),
            });
        }
        if self
            .element_composition
            .iter()
            .any(|coefficient| !coefficient.is_finite() || *coefficient < 0.0)
        {
            return Err(ReactionExtentError::InvalidProblem {
                field: "element_composition",
                message: "elemental coefficients must be finite and non-negative".to_string(),
            });
        }
        for species in 0..species_count {
            if !self
                .element_composition
                .row(species)
                .iter()
                .any(|value| *value > 0.0)
            {
                return Err(ReactionExtentError::InvalidProblem {
                    field: "element_composition",
                    message: format!(
                        "species row {species} must contain at least one positive elemental coefficient"
                    ),
                });
            }
        }
        if self.phases.is_empty() {
            return Err(ReactionExtentError::InvalidProblem {
                field: "phases",
                message: "at least one phase is required".to_string(),
            });
        }
        if self
            .phases
            .iter()
            .flat_map(|phase| phase.species.iter())
            .any(|&index| index >= species_count)
        {
            return Err(ReactionExtentError::InvalidProblem {
                field: "phases",
                message: "a phase references a missing species index".to_string(),
            });
        }
        species_to_phase_map(&self.phases, species_count)?;
        Ok(())
    }

    /// Species names in the solver's deterministic ordering.
    pub fn species(&self) -> &[String] {
        &self.labels
    }

    /// Phase-qualified component descriptors in the deterministic solver ordering.
    pub fn components(&self) -> &[EquilibriumComponentDescriptor] {
        &self.components
    }

    /// Physical initial mole numbers in species order.
    pub fn initial_moles(&self) -> &[f64] {
        &self.initial_moles
    }

    /// Initial iterate in log-mole coordinates.
    pub fn initial_log_moles(&self) -> &LogMolesInitialGuess {
        &self.initial_log_moles
    }

    /// Returns a typed species id in the validated problem ordering.
    pub fn species_id(&self, index: usize) -> Result<SpeciesId, ReactionExtentError> {
        SpeciesId::new(index, self.components.len())
    }

    /// Returns a typed phase id in the validated problem ordering.
    pub fn phase_index(&self, index: usize) -> Result<PhaseIndex, ReactionExtentError> {
        PhaseIndex::new(index, self.phases.len())
    }

    /// Species-by-element composition matrix.
    pub fn element_composition(&self) -> &DMatrix<f64> {
        &self.element_composition
    }

    /// Standard Gibbs-energy functions in species order.
    pub fn gibbs(&self) -> &[GibbsFn] {
        &self.gibbs
    }

    /// Phase layout in species-index order.
    pub fn phases(&self) -> &[Phase] {
        &self.phases
    }

    /// Thermodynamic conditions for this solve.
    pub fn conditions(&self) -> EquilibriumConditions {
        self.conditions
    }

    pub(crate) fn into_parts(
        self,
    ) -> (
        Vec<String>,
        Vec<f64>,
        Vec<f64>,
        DMatrix<f64>,
        Vec<GibbsFn>,
        Vec<Phase>,
        EquilibriumConditions,
    ) {
        (
            self.labels,
            self.initial_moles,
            self.initial_log_moles.into_inner(),
            self.element_composition,
            self.gibbs,
            self.phases,
            self.conditions,
        )
    }
}
