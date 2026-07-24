//! Immutable domain model for independent equilibrium-constant validation.
//!
//! This formulation shares standard-state Gibbs functions with the canonical
//! solver, but computes reaction `ln(K)` and `ln(Q)` directly from a validated
//! reaction basis. It does not call the canonical log-moles residual builder.

use crate::Thermodynamics::ChemEquilibrium::equilibrium_activity::PhaseActivityModel;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_ids::ReactionId;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_log_moles::GibbsFn;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_nonlinear::ReactionExtentError;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_problem::{
    EquilibriumConditions, PreparedEquilibriumProblem,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_reaction_basis::{
    ReactionBasisTolerances, ValidatedReactionBasis,
};

/// Universal gas constant in J/(mol*K).
pub const MOLAR_GAS_CONSTANT: f64 = 8.314_462_618_153_24;

/// Activity models supported by the independent validator.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EquilibriumConstantActivityModel {
    /// Dimensionless ideal-gas activity `a_i = x_i * P / P_ref`.
    IdealGas,
}

/// Immutable input for evaluating reaction equilibrium constants and quotients.
///
/// This problem is intentionally narrower than [`EquilibriumProblem`](super::equilibrium_problem::EquilibriumProblem):
/// it supports only single-reaction systems with one ideal-gas phase. It shares
/// standard-state Gibbs functions with the canonical solver but computes
/// `ln(K)` and `ln(Q)` directly from the reaction basis.
pub struct EquilibriumConstantProblem {
    /// Validated reaction basis with deterministic species ordering.
    basis: ValidatedReactionBasis,
    /// Initial physical mole numbers in basis species order.
    initial_moles: Vec<f64>,
    /// Standard Gibbs free energy functions `g_i(T)` for each species, in J/mol.
    standard_gibbs: Vec<GibbsFn>,
    /// Thermodynamic conditions (T, P, P0) for the evaluation.
    conditions: EquilibriumConditions,
    /// Activity model used for quotient evaluation (currently only IdealGas).
    activity_model: EquilibriumConstantActivityModel,
}

impl EquilibriumConstantProblem {
    /// Derives the independent validator input from the canonical immutable snapshot.
    ///
    /// Only domain data are shared. The canonical residual, Jacobian, scaling,
    /// backend policy, and accepted-candidate report are deliberately absent.
    pub fn from_prepared_ideal_gas(
        prepared: &PreparedEquilibriumProblem,
        basis_tolerances: ReactionBasisTolerances,
    ) -> Result<Self, ReactionExtentError> {
        let source = prepared.problem();
        if source.phases().len() != 1
            || !matches!(source.phases()[0].kind, PhaseActivityModel::IdealGas)
        {
            return Err(invalid_problem(
                "independent validation currently supports exactly one ideal-gas phase",
            ));
        }

        let basis = ValidatedReactionBasis::new(
            source.species().to_vec(),
            source.element_composition(),
            prepared.reaction_basis().reactions.clone(),
            prepared.reaction_basis().rank,
            basis_tolerances,
        )?;
        Self::new(
            basis,
            source.initial_moles().to_vec(),
            source.gibbs().to_vec(),
            source.conditions(),
            EquilibriumConstantActivityModel::IdealGas,
        )
    }

    /// Builds a validated K_eq problem in the reaction basis species order.
    pub fn new(
        basis: ValidatedReactionBasis,
        initial_moles: Vec<f64>,
        standard_gibbs: Vec<GibbsFn>,
        conditions: EquilibriumConditions,
        activity_model: EquilibriumConstantActivityModel,
    ) -> Result<Self, ReactionExtentError> {
        let species_count = basis.species().len();
        if initial_moles.len() != species_count || standard_gibbs.len() != species_count {
            return Err(invalid_problem(format!(
                "basis has {species_count} species, initial moles has {}, and Gibbs functions has {}",
                initial_moles.len(),
                standard_gibbs.len()
            )));
        }
        if initial_moles
            .iter()
            .any(|value| !value.is_finite() || *value < 0.0)
        {
            return Err(invalid_problem(
                "initial mole numbers must be finite and non-negative",
            ));
        }
        if initial_moles.iter().sum::<f64>() <= 0.0 {
            return Err(invalid_problem(
                "at least one initial mole number must be positive",
            ));
        }

        Ok(Self {
            basis,
            initial_moles,
            standard_gibbs,
            conditions,
            activity_model,
        })
    }

    /// Validated reaction basis and deterministic species ordering.
    pub fn basis(&self) -> &ValidatedReactionBasis {
        &self.basis
    }

    /// Initial physical mole numbers in basis order.
    pub fn initial_moles(&self) -> &[f64] {
        &self.initial_moles
    }

    /// Thermodynamic conditions used by both `ln(K)` and `ln(Q)`.
    pub fn conditions(&self) -> EquilibriumConditions {
        self.conditions
    }

    /// Activity model selected for quotient evaluation.
    pub fn activity_model(&self) -> EquilibriumConstantActivityModel {
        self.activity_model
    }

    /// Standard-state Gibbs-energy closures in basis order.
    pub fn standard_gibbs(&self) -> &[GibbsFn] {
        &self.standard_gibbs
    }

    /// Standard reaction Gibbs energy in J/mol for one independent reaction.
    pub fn standard_reaction_gibbs(
        &self,
        reaction: ReactionId,
    ) -> Result<f64, ReactionExtentError> {
        let coefficients = self.basis.reaction(reaction)?;
        let temperature = self.conditions.temperature();
        coefficients
            .iter()
            .zip(&self.standard_gibbs)
            .enumerate()
            .try_fold(0.0, |sum, (species, (nu, gibbs))| {
                let value = gibbs(temperature);
                if !value.is_finite() {
                    return Err(ReactionExtentError::InvalidDG0 {
                        species_index: species,
                        dg0: value,
                        temperature,
                    });
                }
                Ok(sum + nu * value)
            })
    }

    /// Natural logarithm of the dimensionless equilibrium constant.
    pub fn ln_equilibrium_constant(
        &self,
        reaction: ReactionId,
    ) -> Result<f64, ReactionExtentError> {
        Ok(-self.standard_reaction_gibbs(reaction)?
            / (MOLAR_GAS_CONSTANT * self.conditions.temperature()))
    }

    /// Natural logarithms of all independent equilibrium constants.
    pub fn ln_equilibrium_constants(&self) -> Result<Vec<f64>, ReactionExtentError> {
        (0..self.basis.reaction_count())
            .map(|index| {
                let id = self.basis.reaction_id(index)?;
                self.ln_equilibrium_constant(id)
            })
            .collect()
    }

    /// Natural logarithm of the reaction quotient at candidate mole numbers.
    pub fn ln_reaction_quotient(
        &self,
        reaction: ReactionId,
        moles: &[f64],
    ) -> Result<f64, ReactionExtentError> {
        let log_activities = self.log_activities(moles)?;
        let coefficients = self.basis.reaction(reaction)?;
        Ok(coefficients
            .iter()
            .zip(log_activities)
            .map(|(nu, ln_activity)| nu * ln_activity)
            .sum())
    }

    /// Returns `ln(Q) - ln(K)` for every independent reaction.
    pub fn equilibrium_residuals(&self, moles: &[f64]) -> Result<Vec<f64>, ReactionExtentError> {
        (0..self.basis.reaction_count())
            .map(|index| {
                let reaction = self.basis.reaction_id(index)?;
                Ok(self.ln_reaction_quotient(reaction, moles)?
                    - self.ln_equilibrium_constant(reaction)?)
            })
            .collect()
    }

    fn log_activities(&self, moles: &[f64]) -> Result<Vec<f64>, ReactionExtentError> {
        if moles.len() != self.basis.species().len() {
            return Err(invalid_problem(format!(
                "candidate has {} mole numbers for {} species",
                moles.len(),
                self.basis.species().len()
            )));
        }
        if moles
            .iter()
            .any(|value| !value.is_finite() || *value <= 0.0)
        {
            return Err(invalid_problem(
                "reaction quotients require finite, strictly positive candidate moles",
            ));
        }
        let total: f64 = moles.iter().sum();
        if !total.is_finite() || total <= 0.0 {
            return Err(invalid_problem(
                "candidate total moles must be finite and positive",
            ));
        }

        match self.activity_model {
            EquilibriumConstantActivityModel::IdealGas => {
                let pressure_ratio =
                    self.conditions.pressure() / self.conditions.reference_pressure();
                Ok(moles
                    .iter()
                    .map(|moles_i| (moles_i / total * pressure_ratio).ln())
                    .collect())
            }
        }
    }
}

fn invalid_problem(message: impl Into<String>) -> ReactionExtentError {
    ReactionExtentError::InvalidProblem {
        field: "equilibrium_constant_problem",
        message: message.into(),
    }
}
