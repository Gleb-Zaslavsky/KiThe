//! Typed reaction-basis contract shared by equilibrium-constant validation.
//!
//! Rows always follow the declared species order and columns represent
//! independent reactions. Construction verifies elemental conservation and
//! applies a deterministic sign and scale convention to every reaction.

use crate::Thermodynamics::ChemEquilibrium::equilibrium_ids::{ReactionId, SpeciesId};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_nonlinear::ReactionExtentError;
use nalgebra::DMatrix;
use std::collections::HashSet;

/// Numerical tolerances used while accepting and normalizing a reaction basis.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ReactionBasisTolerances {
    /// Largest accepted coefficient in `A^T * N`.
    pub conservation: f64,
    /// Coefficients at or below this magnitude are treated as numerical zero.
    pub coefficient_zero: f64,
}

impl Default for ReactionBasisTolerances {
    fn default() -> Self {
        Self {
            conservation: 1e-10,
            coefficient_zero: 1e-12,
        }
    }
}

impl ReactionBasisTolerances {
    fn validate(self) -> Result<Self, ReactionExtentError> {
        if !self.conservation.is_finite() || self.conservation <= 0.0 {
            return Err(invalid_basis(
                "conservation tolerance must be finite and strictly positive",
            ));
        }
        if !self.coefficient_zero.is_finite() || self.coefficient_zero <= 0.0 {
            return Err(invalid_basis(
                "coefficient-zero tolerance must be finite and strictly positive",
            ));
        }
        Ok(self)
    }
}

/// Validated independent reactions aligned to an explicit species ordering.
///
/// Construction verifies that `A^T · N ≈ 0` (elemental conservation) within
/// the configured tolerance and applies a deterministic sign and scale
/// convention to every reaction column.
#[derive(Debug, Clone, PartialEq)]
pub struct ValidatedReactionBasis {
    /// Species names in the exact row order used by the reaction matrix.
    species: Vec<String>,
    /// Number of chemical elements (columns in the composition matrix).
    element_count: usize,
    /// Rank of the element composition matrix (number of independent elements).
    element_rank: usize,
    /// Stoichiometric matrix `N` of shape `(species × reactions)`.
    /// Each column is one independent reaction direction.
    reactions: DMatrix<f64>,
    /// Maximum absolute value in `A^T · N` — the worst elemental conservation error.
    max_conservation_residual: f64,
}

impl ValidatedReactionBasis {
    /// Validates and deterministically normalizes an explicit reaction basis.
    ///
    /// `element_composition` has species in rows and elements in columns;
    /// `reactions` has the same species in rows and reactions in columns.
    pub fn new(
        species: Vec<String>,
        element_composition: &DMatrix<f64>,
        reactions: DMatrix<f64>,
        element_rank: usize,
        tolerances: ReactionBasisTolerances,
    ) -> Result<Self, ReactionExtentError> {
        let tolerances = tolerances.validate()?;
        validate_species(&species)?;

        let species_count = species.len();
        if element_composition.nrows() != species_count || element_composition.ncols() == 0 {
            return Err(invalid_basis(format!(
                "element matrix must have {species_count} rows and at least one column"
            )));
        }
        if reactions.nrows() != species_count {
            return Err(invalid_basis(format!(
                "reaction matrix has {} rows for {species_count} species",
                reactions.nrows()
            )));
        }
        if element_rank > species_count.min(element_composition.ncols()) {
            return Err(invalid_basis(format!(
                "element rank {element_rank} exceeds matrix dimensions"
            )));
        }
        let expected_reactions = species_count - element_rank;
        if reactions.ncols() != expected_reactions {
            return Err(invalid_basis(format!(
                "reaction matrix has {} columns; expected {expected_reactions} from species count and element rank",
                reactions.ncols()
            )));
        }
        if element_composition.iter().any(|value| !value.is_finite())
            || reactions.iter().any(|value| !value.is_finite())
        {
            return Err(invalid_basis(
                "basis matrices must contain only finite values",
            ));
        }

        let reactions = normalize_reactions(reactions, tolerances.coefficient_zero)?;
        let conservation = element_composition.transpose() * &reactions;
        let max_conservation_residual = conservation
            .iter()
            .fold(0.0_f64, |largest, value| largest.max(value.abs()));
        if max_conservation_residual > tolerances.conservation {
            return Err(invalid_basis(format!(
                "reaction basis violates elemental conservation: max |A^T*N| = {max_conservation_residual:e}"
            )));
        }

        Ok(Self {
            species,
            element_count: element_composition.ncols(),
            element_rank,
            reactions,
            max_conservation_residual,
        })
    }

    /// Species names in the exact row order used by the reaction matrix.
    pub fn species(&self) -> &[String] {
        &self.species
    }

    /// Number of element columns used to validate conservation.
    pub fn element_count(&self) -> usize {
        self.element_count
    }

    /// Numerical rank assigned to the element-composition matrix.
    pub fn element_rank(&self) -> usize {
        self.element_rank
    }

    /// Number of independent reaction columns.
    pub fn reaction_count(&self) -> usize {
        self.reactions.ncols()
    }

    /// Deterministically normalized species-by-reaction matrix.
    pub fn reactions(&self) -> &DMatrix<f64> {
        &self.reactions
    }

    /// Largest absolute coefficient in `A^T * N` at construction time.
    pub fn max_conservation_residual(&self) -> f64 {
        self.max_conservation_residual
    }

    /// Returns one reaction column using a typed identifier.
    pub fn reaction(&self, id: ReactionId) -> Result<Vec<f64>, ReactionExtentError> {
        if id.index() >= self.reaction_count() {
            return Err(invalid_basis(format!(
                "reaction index {} is out of bounds for {} reactions",
                id.index(),
                self.reaction_count()
            )));
        }
        Ok(self.reactions.column(id.index()).iter().copied().collect())
    }

    /// Creates a species identifier in this basis ordering.
    pub fn species_id(&self, index: usize) -> Result<SpeciesId, ReactionExtentError> {
        SpeciesId::new(index, self.species.len())
    }

    /// Creates a reaction identifier in this basis ordering.
    pub fn reaction_id(&self, index: usize) -> Result<ReactionId, ReactionExtentError> {
        ReactionId::new(index, self.reaction_count())
    }
}

fn validate_species(species: &[String]) -> Result<(), ReactionExtentError> {
    if species.is_empty() || species.iter().any(|name| name.trim().is_empty()) {
        return Err(invalid_basis("species names must be non-empty"));
    }
    let unique: HashSet<&str> = species.iter().map(String::as_str).collect();
    if unique.len() != species.len() {
        return Err(invalid_basis("species names must be unique"));
    }
    Ok(())
}

fn normalize_reactions(
    mut reactions: DMatrix<f64>,
    zero_tolerance: f64,
) -> Result<DMatrix<f64>, ReactionExtentError> {
    for reaction in 0..reactions.ncols() {
        let pivot = (0..reactions.nrows())
            .find(|&species| reactions[(species, reaction)].abs() > zero_tolerance)
            .ok_or_else(|| {
                invalid_basis(format!("reaction column {reaction} is numerically zero"))
            })?;
        let scale = reactions[(pivot, reaction)].abs();
        let sign = if reactions[(pivot, reaction)] > 0.0 {
            -1.0
        } else {
            1.0
        };
        for species in 0..reactions.nrows() {
            let normalized = sign * reactions[(species, reaction)] / scale;
            reactions[(species, reaction)] = if normalized.abs() <= zero_tolerance {
                0.0
            } else {
                normalized
            };
        }
    }
    Ok(reactions)
}

fn invalid_basis(message: impl Into<String>) -> ReactionExtentError {
    ReactionExtentError::InvalidProblem {
        field: "reaction_basis",
        message: message.into(),
    }
}
