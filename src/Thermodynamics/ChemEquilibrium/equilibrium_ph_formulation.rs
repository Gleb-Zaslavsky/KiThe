//! Structural contracts for the monolithic fixed-`P,H` formulation.
//!
//! This module does not run a nonlinear solver. It owns the deterministic
//! unknown/residual ordering and the smooth transformation that keeps the
//! solved temperature inside the common thermochemistry interval.
//!
//! The monolithic runner consumes these contracts through the common backend
//! adapter. Formulation tests therefore cover the block dimensions,
//! thermochemical capabilities, and analytic Jacobian independently from the
//! backend cascade.

#![allow(dead_code)]

use std::ops::Range;

use nalgebra::DMatrix;

use crate::Thermodynamics::ChemEquilibrium::equilibrium_activity::phase_activity_models;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_constraints::TemperatureBounds;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_constraints::{
    EnthalpyScale, additive_total_enthalpy,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_log_moles::{
    R, compute_species_moles, evaluate_equilibrium_logmole_residual_with_standard_gibbs,
    scale_jacobian_rows, scale_residual_rows,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_nonlinear::ReactionExtentError;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_ph_thermochemistry::ResolvedThermochemistry;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_problem::PreparedEquilibriumProblem;

/// Column ordering of the monolithic unknown vector.
///
/// Species log-moles retain the canonical fixed-`P,T` ordering. The final
/// scalar coordinate is transformed into the physical temperature.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct PhUnknownLayout {
    species_count: usize,
}

impl PhUnknownLayout {
    /// Creates a non-empty monolithic unknown layout.
    pub(crate) fn new(species_count: usize) -> Result<Self, ReactionExtentError> {
        if species_count == 0 {
            return Err(ReactionExtentError::InvalidProblem {
                field: "ph_unknown_layout",
                message: "at least one species coordinate is required".to_string(),
            });
        }
        Ok(Self { species_count })
    }

    /// Log-mole columns in canonical component order.
    pub(crate) fn log_moles(self) -> Range<usize> {
        0..self.species_count
    }

    /// Column containing the unconstrained temperature coordinate.
    pub(crate) fn temperature(self) -> usize {
        self.species_count
    }

    /// Complete nonlinear unknown dimension.
    pub(crate) fn dimension(self) -> usize {
        self.species_count + 1
    }
}

/// Row ordering of the monolithic residual.
///
/// The existing fixed-`P,T` reaction and element rows remain unchanged. One
/// scaled total-enthalpy row is appended at the end.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct PhResidualLayout {
    reaction_count: usize,
    element_count: usize,
}

impl PhResidualLayout {
    /// Creates a square P,H layout for the supplied species count.
    pub(crate) fn new(
        species_count: usize,
        reaction_count: usize,
        element_count: usize,
    ) -> Result<Self, ReactionExtentError> {
        if reaction_count + element_count != species_count {
            return Err(ReactionExtentError::DimensionMismatch(format!(
                "P,H residual has {reaction_count} reaction rows and {element_count} element rows for {species_count} species"
            )));
        }
        Ok(Self {
            reaction_count,
            element_count,
        })
    }

    /// Reaction-equilibrium rows inherited from fixed P,T.
    pub(crate) fn reactions(self) -> Range<usize> {
        0..self.reaction_count
    }

    /// Element-conservation rows inherited from fixed P,T.
    pub(crate) fn elements(self) -> Range<usize> {
        self.reaction_count..self.reaction_count + self.element_count
    }

    /// Final scaled total-enthalpy row.
    pub(crate) fn enthalpy(self) -> usize {
        self.reaction_count + self.element_count
    }

    /// Complete nonlinear residual dimension.
    pub(crate) fn dimension(self) -> usize {
        self.enthalpy() + 1
    }
}

/// Physical temperature and its chain-rule factor at one solver coordinate.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct TransformedTemperature {
    temperature: f64,
    derivative: f64,
}

impl TransformedTemperature {
    /// Physical temperature in kelvin.
    pub(crate) fn temperature(self) -> f64 {
        self.temperature
    }

    /// `dT/dtheta` for the monolithic Jacobian temperature column.
    pub(crate) fn derivative(self) -> f64 {
        self.derivative
    }
}

/// Smooth logistic map from one unconstrained coordinate into strict bounds.
///
/// Saturated floating-point coordinates are rejected instead of being
/// silently clipped to a thermochemistry boundary. This keeps residual and
/// Jacobian evaluations mathematically consistent.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct BoundedTemperatureTransform {
    bounds: TemperatureBounds,
    span: f64,
}

impl BoundedTemperatureTransform {
    /// Creates a transform for a non-degenerate common temperature interval.
    pub(crate) fn new(bounds: TemperatureBounds) -> Result<Self, ReactionExtentError> {
        let span = bounds.upper() - bounds.lower();
        if !span.is_finite() || span <= 0.0 {
            return Err(ReactionExtentError::InvalidProblem {
                field: "ph_temperature_bounds",
                message: "monolithic P,H requires a non-degenerate temperature interval"
                    .to_string(),
            });
        }
        Ok(Self { bounds, span })
    }

    /// Maps an unconstrained coordinate to temperature and `dT/dtheta`.
    pub(crate) fn evaluate(
        self,
        theta: f64,
    ) -> Result<TransformedTemperature, ReactionExtentError> {
        if !theta.is_finite() {
            return Err(ReactionExtentError::InvalidProblem {
                field: "ph_temperature_coordinate",
                message: "temperature coordinate must be finite".to_string(),
            });
        }

        // Evaluate the logistic without overflowing for large negative theta.
        let fraction = if theta >= 0.0 {
            1.0 / (1.0 + (-theta).exp())
        } else {
            let exponential = theta.exp();
            exponential / (1.0 + exponential)
        };
        let temperature = self.bounds.lower() + self.span * fraction;
        let derivative = self.span * fraction * (1.0 - fraction);
        if !temperature.is_finite()
            || temperature <= self.bounds.lower()
            || temperature >= self.bounds.upper()
            || !derivative.is_finite()
            || derivative <= 0.0
        {
            return Err(ReactionExtentError::InvalidCandidate {
                field: "ph_temperature_coordinate",
                message: format!(
                    "temperature coordinate {theta} saturated the interval [{}, {}] K",
                    self.bounds.lower(),
                    self.bounds.upper()
                ),
            });
        }
        Ok(TransformedTemperature {
            temperature,
            derivative,
        })
    }

    /// Converts a strict interior temperature into the solver coordinate.
    pub(crate) fn coordinate(self, temperature: f64) -> Result<f64, ReactionExtentError> {
        if !temperature.is_finite()
            || temperature <= self.bounds.lower()
            || temperature >= self.bounds.upper()
        {
            return Err(ReactionExtentError::InvalidProblem {
                field: "ph_temperature_seed",
                message: format!(
                    "monolithic P,H seed must lie strictly inside [{}, {}] K",
                    self.bounds.lower(),
                    self.bounds.upper()
                ),
            });
        }
        let fraction = (temperature - self.bounds.lower()) / self.span;
        let coordinate = (fraction / (1.0 - fraction)).ln();
        if !coordinate.is_finite() {
            return Err(ReactionExtentError::InvalidProblem {
                field: "ph_temperature_seed",
                message: "temperature seed produced a non-finite coordinate".to_string(),
            });
        }
        Ok(coordinate)
    }

    /// Converts an inclusive public temperature seed into a strict interior
    /// seed for the logistic coordinate. A boundary value is valid user input
    /// for the P,H request, but has no finite inverse through this transform.
    ///
    /// This normalization happens once while constructing the initial guess;
    /// residual evaluation still rejects saturated coordinates and never clips
    /// a nonlinear iterate to the thermochemistry bounds.
    pub(crate) fn interior_seed(self, temperature: f64) -> Result<f64, ReactionExtentError> {
        if !temperature.is_finite() || !self.bounds.contains(temperature) {
            return Err(ReactionExtentError::InvalidProblem {
                field: "ph_temperature_seed",
                message: format!(
                    "monolithic P,H seed must lie inside [{}, {}] K",
                    self.bounds.lower(),
                    self.bounds.upper()
                ),
            });
        }
        if temperature > self.bounds.lower() && temperature < self.bounds.upper() {
            return Ok(temperature);
        }
        let midpoint = self.bounds.lower() + self.span * 0.5;
        if midpoint <= self.bounds.lower() || midpoint >= self.bounds.upper() {
            return Err(ReactionExtentError::InvalidProblem {
                field: "ph_temperature_seed",
                message: "monolithic P,H bounds cannot produce a strict interior seed".to_string(),
            });
        }
        Ok(midpoint)
    }
}

/// Immutable fixed-active-set P,H equations before backend integration.
///
/// This object is deliberately solver agnostic. It reuses the canonical P,T
/// residual and analytical log-mole Jacobian, then appends the temperature
/// column and total-enthalpy row.
#[derive(Clone)]
pub(crate) struct PreparedPhFormulation {
    prepared_pt: PreparedEquilibriumProblem,
    thermochemistry: ResolvedThermochemistry,
    temperature_transform: BoundedTemperatureTransform,
    unknown_layout: PhUnknownLayout,
    residual_layout: PhResidualLayout,
    pt_row_scale: Vec<f64>,
    target_enthalpy: f64,
    enthalpy_scale: EnthalpyScale,
}

struct PhEvaluationContext {
    log_moles: Vec<f64>,
    moles: Vec<f64>,
    temperature: TransformedTemperature,
    standard_gibbs: Vec<f64>,
    molar_enthalpies: Vec<f64>,
    heat_capacities: Vec<f64>,
}

/// Evaluated monolithic P,H candidate before common acceptance checks.
///
/// The snapshot keeps the unscaled fixed-P,T residual separate from the
/// solver-facing scaled vector. That distinction is essential: numerical row
/// scaling must never weaken the physical conservation and affinity checks.
#[derive(Debug, Clone)]
pub(crate) struct PhCandidateSnapshot {
    pub(crate) log_moles: Vec<f64>,
    pub(crate) moles: Vec<f64>,
    pub(crate) temperature: f64,
    pub(crate) raw_pt_residual: Vec<f64>,
    pub(crate) scaled_pt_residual: Vec<f64>,
    pub(crate) total_enthalpy: f64,
    pub(crate) enthalpy_error: f64,
    pub(crate) scaled_enthalpy_error: f64,
}

impl PreparedPhFormulation {
    /// Prepares one monolithic formulation for a fixed active phase set.
    pub(crate) fn new(
        prepared_pt: PreparedEquilibriumProblem,
        thermochemistry: ResolvedThermochemistry,
        temperature_bounds: TemperatureBounds,
        target_enthalpy: f64,
        enthalpy_scale: EnthalpyScale,
    ) -> Result<Self, ReactionExtentError> {
        if !target_enthalpy.is_finite() {
            return Err(ReactionExtentError::InvalidProblem {
                field: "target_enthalpy",
                message: "monolithic P,H target enthalpy must be finite".to_string(),
            });
        }
        let species_count = prepared_pt.problem().species().len();
        if thermochemistry.len() != species_count {
            return Err(ReactionExtentError::DimensionMismatch(format!(
                "P,H thermochemistry has {} components for {species_count} prepared species",
                thermochemistry.len()
            )));
        }
        let supported = thermochemistry.temperature_bounds();
        if temperature_bounds.lower() < supported.lower()
            || temperature_bounds.upper() > supported.upper()
        {
            return Err(ReactionExtentError::InvalidProblem {
                field: "temperature_bounds",
                message: "monolithic P,H bounds extend beyond thermochemistry capabilities"
                    .to_string(),
            });
        }
        // A monolithic P,H Jacobian contains the explicit temperature
        // derivative of total enthalpy. Validate the complete thermochemical
        // capability before a backend starts, rather than allowing a missing
        // Cp entry to surface as an opaque nonlinear failure at the first
        // iterate. The lower bound is inside the already checked common
        // interval and is deterministic for polynomial-boundary fixtures.
        let preparation_temperature = temperature_bounds.lower();
        thermochemistry.evaluate_gibbs(preparation_temperature)?;
        thermochemistry.evaluate_enthalpy(preparation_temperature)?;
        let heat_capacity = thermochemistry.evaluate_heat_capacity(preparation_temperature)?;
        if let Some(index) = heat_capacity.iter().position(Option::is_none) {
            return Err(ReactionExtentError::InvalidProblem {
                field: "heat_capacity",
                message: format!(
                    "monolithic P,H requires a Cp capability for component {index}"
                ),
            });
        }
        let reaction_count = prepared_pt.reaction_basis().reactions.ncols();
        let element_count = prepared_pt.problem().element_composition().ncols();
        let unknown_layout = PhUnknownLayout::new(species_count)?;
        let residual_layout = PhResidualLayout::new(species_count, reaction_count, element_count)?;
        let pt_row_scale = prepared_pt.residual_scale()?;
        if pt_row_scale.len() != species_count {
            return Err(ReactionExtentError::DimensionMismatch(format!(
                "fixed-P,T scale has {} rows for {species_count} species",
                pt_row_scale.len()
            )));
        }

        Ok(Self {
            prepared_pt,
            thermochemistry,
            temperature_transform: BoundedTemperatureTransform::new(temperature_bounds)?,
            unknown_layout,
            residual_layout,
            pt_row_scale,
            target_enthalpy,
            enthalpy_scale,
        })
    }

    /// Builds `[ln(n), theta_T]` from the prepared P,T seed and a temperature.
    pub(crate) fn initial_unknowns(
        &self,
        temperature: f64,
    ) -> Result<Vec<f64>, ReactionExtentError> {
        self.initial_unknowns_from_log_moles(
            self.prepared_pt.problem().initial_log_moles(),
            temperature,
        )
    }

    /// Builds `[ln(n), theta_T]` from an accepted physical composition and
    /// temperature. This is the only monolithic continuation entry point:
    /// callers must provide a validated accepted state, never a failed
    /// backend iterate.
    pub(crate) fn initial_unknowns_from_log_moles(
        &self,
        log_moles: &crate::Thermodynamics::ChemEquilibrium::equilibrium_problem::
            LogMolesInitialGuess,
        temperature: f64,
    ) -> Result<Vec<f64>, ReactionExtentError> {
        if log_moles.as_slice().len() != self.unknown_layout.species_count {
            return Err(ReactionExtentError::DimensionMismatch(format!(
                "continuation seed has {} entries for {} P,H components",
                log_moles.as_slice().len(),
                self.unknown_layout.species_count
            )));
        }
        let mut unknowns = log_moles.as_slice().to_vec();
        unknowns.push(
            self.temperature_transform
                .coordinate(self.temperature_transform.interior_seed(temperature)?)?,
        );
        Ok(unknowns)
    }

    /// Retargets only the mutable boundary data of a prepared formulation.
    /// Reaction basis, element totals, row scales, layout, and temperature
    /// transform remain shared structural state.
    pub(crate) fn retarget(
        &self,
        prepared_pt: PreparedEquilibriumProblem,
        target_enthalpy: f64,
        enthalpy_scale: EnthalpyScale,
    ) -> Result<Self, ReactionExtentError> {
        if prepared_pt.problem().species().len() != self.unknown_layout.species_count {
            return Err(ReactionExtentError::DimensionMismatch(
                "retargeted P,H preparation changed component count".to_string(),
            ));
        }
        if !target_enthalpy.is_finite() {
            return Err(ReactionExtentError::InvalidProblem {
                field: "target_enthalpy",
                message: "retargeted P,H target enthalpy must be finite".to_string(),
            });
        }
        Ok(Self {
            prepared_pt,
            thermochemistry: self.thermochemistry.clone(),
            temperature_transform: self.temperature_transform,
            unknown_layout: self.unknown_layout,
            residual_layout: self.residual_layout,
            pt_row_scale: self.pt_row_scale.clone(),
            target_enthalpy,
            enthalpy_scale,
        })
    }

    /// Evaluates the scaled monolithic P,H residual.
    pub(crate) fn residual(&self, unknowns: &[f64]) -> Result<Vec<f64>, ReactionExtentError> {
        let snapshot = self.candidate_snapshot(unknowns)?;
        let mut residual = snapshot.scaled_pt_residual;
        residual.push(snapshot.scaled_enthalpy_error);
        debug_assert_eq!(residual.len(), self.residual_layout.dimension());
        Ok(residual)
    }

    /// Evaluates one candidate without selecting or publishing a solution.
    pub(crate) fn candidate_snapshot(
        &self,
        unknowns: &[f64],
    ) -> Result<PhCandidateSnapshot, ReactionExtentError> {
        let context = self.evaluate_context(unknowns)?;
        let problem = self.prepared_pt.problem();
        let conditions = problem.conditions();
        let raw_pt_residual = evaluate_equilibrium_logmole_residual_with_standard_gibbs(
            &context.log_moles,
            &self.prepared_pt.reaction_basis().reactions,
            problem.element_composition(),
            self.prepared_pt.element_totals(),
            &context.standard_gibbs,
            problem.phases(),
            context.temperature.temperature(),
            conditions.pressure(),
            conditions.reference_pressure(),
            self.prepared_pt.species_phase(),
            self.prepared_pt.phase_stoichiometry(),
        )?;
        let scaled_pt_residual = scale_residual_rows(raw_pt_residual.clone(), &self.pt_row_scale)?;
        let total_enthalpy = additive_total_enthalpy(&context.moles, &context.molar_enthalpies)?;
        let enthalpy_error = total_enthalpy - self.target_enthalpy;
        let scaled_enthalpy_error = self.enthalpy_scale.scale_error(enthalpy_error)?;
        Ok(PhCandidateSnapshot {
            log_moles: context.log_moles,
            moles: context.moles,
            temperature: context.temperature.temperature(),
            raw_pt_residual,
            scaled_pt_residual,
            total_enthalpy,
            enthalpy_error,
            scaled_enthalpy_error,
        })
    }

    /// Tests whether an iterate can safely be evaluated by the fixed-active
    /// P,H equations. Detailed diagnostics remain in `candidate_snapshot`.
    pub(crate) fn is_feasible(&self, unknowns: &[f64]) -> bool {
        self.candidate_snapshot(unknowns)
            .map(|snapshot| {
                snapshot
                    .moles
                    .iter()
                    .all(|moles| moles.is_finite() && *moles >= 0.0)
            })
            .unwrap_or(false)
    }

    /// Underlying fixed-P,T problem used by the shared acceptance gate.
    pub(crate) fn prepared_pt(&self) -> &PreparedEquilibriumProblem {
        &self.prepared_pt
    }

    /// Scale used by the total-enthalpy acceptance contract.
    pub(crate) fn enthalpy_scale(&self) -> EnthalpyScale {
        self.enthalpy_scale
    }

    /// Row scales shared by the analytic and RST-symbolic P,H residuals.
    pub(crate) fn pt_row_scale(&self) -> &[f64] {
        &self.pt_row_scale
    }

    /// Prescribed total enthalpy in joules.
    pub(crate) fn target_enthalpy(&self) -> f64 {
        self.target_enthalpy
    }

    /// Resolved thermochemistry aligned with the prepared active system.
    pub(crate) fn thermochemistry(&self) -> &ResolvedThermochemistry {
        &self.thermochemistry
    }

    /// Common thermochemistry interval enforced by the bounded temperature
    /// coordinate.
    pub(crate) fn temperature_bounds(&self) -> TemperatureBounds {
        self.temperature_transform.bounds
    }

    /// Evaluates the full analytical Jacobian with respect to `[ln(n), theta_T]`.
    pub(crate) fn jacobian(&self, unknowns: &[f64]) -> Result<DMatrix<f64>, ReactionExtentError> {
        let context = self.evaluate_context(unknowns)?;
        let problem = self.prepared_pt.problem();
        let conditions = problem.conditions();
        let pt_jacobian = scale_jacobian_rows(
            self.prepared_pt.jacobian(&context.log_moles)?,
            &self.pt_row_scale,
        )?;
        let dimension = self.unknown_layout.dimension();
        let mut jacobian = DMatrix::zeros(dimension, dimension);
        for row in 0..pt_jacobian.nrows() {
            for column in 0..pt_jacobian.ncols() {
                jacobian[(row, column)] = pt_jacobian[(row, column)];
            }
        }

        let temperature = context.temperature.temperature();
        let d_temperature_d_coordinate = context.temperature.derivative();
        let reaction_basis = &self.prepared_pt.reaction_basis().reactions;
        let activity_models = phase_activity_models(problem.phases());
        for reaction in self.residual_layout.reactions() {
            let mut derivative = 0.0;
            for species in self.unknown_layout.log_moles() {
                let coefficient = reaction_basis[(species, reaction)];
                if coefficient == 0.0 {
                    continue;
                }
                derivative -= coefficient * context.molar_enthalpies[species]
                    / (R * temperature * temperature);
                let phase = self.prepared_pt.species_phase()[species];
                derivative += coefficient
                    * activity_models[phase].d_log_activity_d_temperature(
                        conditions.pressure(),
                        conditions.reference_pressure(),
                    )?;
            }
            jacobian[(reaction, self.unknown_layout.temperature())] =
                derivative * d_temperature_d_coordinate / self.pt_row_scale[reaction];
        }

        let enthalpy_row = self.residual_layout.enthalpy();
        for species in self.unknown_layout.log_moles() {
            jacobian[(enthalpy_row, species)] = context.moles[species]
                * context.molar_enthalpies[species]
                / self.enthalpy_scale.joules();
        }
        let partial_temperature_derivative = context
            .moles
            .iter()
            .zip(&context.heat_capacities)
            .map(|(&moles, &heat_capacity)| moles * heat_capacity)
            .sum::<f64>();
        if !partial_temperature_derivative.is_finite() {
            return Err(ReactionExtentError::JacobianEvaluation(
                "P,H enthalpy temperature derivative is not finite".to_string(),
            ));
        }
        jacobian[(enthalpy_row, self.unknown_layout.temperature())] = partial_temperature_derivative
            * d_temperature_d_coordinate
            / self.enthalpy_scale.joules();
        Ok(jacobian)
    }

    /// Decomposes the monolithic `[log-moles, theta_T]` unknown vector into a
    /// structured [`PhEvaluationContext`] containing physical moles, the solved
    /// temperature, and the enthalpy-scale-normalized energy residual.
    ///
    /// The context is consumed by both the analytical residual/Jacobian
    /// formulation and the candidate-snapshot builder. A dimension mismatch
    /// between `unknowns` and the declared layout is caught here rather than
    /// propagating an index-out-of-bounds panic.
    fn evaluate_context(
        &self,
        unknowns: &[f64],
    ) -> Result<PhEvaluationContext, ReactionExtentError> {
        if unknowns.len() != self.unknown_layout.dimension() {
            return Err(ReactionExtentError::DimensionMismatch(format!(
                "P,H iterate has {} entries, expected {}",
                unknowns.len(),
                self.unknown_layout.dimension()
            )));
        }
        let log_moles = unknowns[self.unknown_layout.log_moles()].to_vec();
        let moles = compute_species_moles(&log_moles)?;
        let temperature = self
            .temperature_transform
            .evaluate(unknowns[self.unknown_layout.temperature()])?;
        let standard_gibbs = self
            .thermochemistry
            .evaluate_gibbs(temperature.temperature())?;
        let molar_enthalpies = self
            .thermochemistry
            .evaluate_enthalpy(temperature.temperature())?;
        let heat_capacities = self
            .thermochemistry
            .evaluate_heat_capacity(temperature.temperature())?
            .into_iter()
            .enumerate()
            .map(|(index, value)| {
                value.ok_or_else(|| ReactionExtentError::InvalidProblem {
                    field: "heat_capacity",
                    message: format!("monolithic P,H requires heat-capacity capability {index}"),
                })
            })
            .collect::<Result<Vec<_>, _>>()?;
        Ok(PhEvaluationContext {
            log_moles,
            moles,
            temperature,
            standard_gibbs,
            molar_enthalpies,
            heat_capacities,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::rc::Rc;
    use std::sync::Arc;

    use crate::Thermodynamics::ChemEquilibrium::equilibrium_log_moles::{
        GibbsFn, Phase, PhaseKind,
    };
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_ph_workflow::{
        MolarThermoFunction, ThermochemistryProvenance,
    };
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_problem::{
        EquilibriumConditions, EquilibriumProblem, LogMolesInitialGuess,
    };
    use crate::Thermodynamics::phase_layout::{PhaseComponentId, PhaseId};

    #[test]
    fn ph_layout_appends_one_temperature_column_and_one_enthalpy_row() {
        let unknowns = PhUnknownLayout::new(5).unwrap();
        let residuals = PhResidualLayout::new(5, 3, 2).unwrap();

        assert_eq!(unknowns.log_moles(), 0..5);
        assert_eq!(unknowns.temperature(), 5);
        assert_eq!(unknowns.dimension(), 6);
        assert_eq!(residuals.reactions(), 0..3);
        assert_eq!(residuals.elements(), 3..5);
        assert_eq!(residuals.enthalpy(), 5);
        assert_eq!(residuals.dimension(), 6);
    }

    #[test]
    fn ph_layout_rejects_a_non_square_fixed_pt_block() {
        assert!(PhUnknownLayout::new(0).is_err());
        assert!(PhResidualLayout::new(5, 2, 2).is_err());
    }

    #[test]
    fn bounded_temperature_transform_round_trips_interior_values() {
        let transform =
            BoundedTemperatureTransform::new(TemperatureBounds::new(300.0, 3_000.0).unwrap())
                .unwrap();

        for expected in [301.0, 500.0, 1_650.0, 2_999.0] {
            let coordinate = transform.coordinate(expected).unwrap();
            let state = transform.evaluate(coordinate).unwrap();
            assert!((state.temperature() - expected).abs() < 1.0e-10);
            assert!(state.derivative() > 0.0);
        }
    }

    #[test]
    fn bounded_temperature_transform_has_the_expected_center_derivative() {
        let transform =
            BoundedTemperatureTransform::new(TemperatureBounds::new(300.0, 1_300.0).unwrap())
                .unwrap();
        let state = transform.evaluate(0.0).unwrap();

        assert_eq!(state.temperature(), 800.0);
        assert_eq!(state.derivative(), 250.0);
    }

    #[test]
    fn bounded_temperature_transform_rejects_bounds_and_saturation() {
        assert!(
            BoundedTemperatureTransform::new(TemperatureBounds::new(500.0, 500.0).unwrap())
                .is_err()
        );
        let transform =
            BoundedTemperatureTransform::new(TemperatureBounds::new(300.0, 3_000.0).unwrap())
                .unwrap();
        assert!(transform.coordinate(300.0).is_err());
        assert!(transform.coordinate(3_000.0).is_err());
        assert!(transform.evaluate(1_000.0).is_err());
        assert!(transform.evaluate(f64::NAN).is_err());
    }

    #[test]
    fn boundary_temperature_seed_is_normalized_only_before_the_first_iterate() {
        let transform =
            BoundedTemperatureTransform::new(TemperatureBounds::new(300.0, 1_300.0).unwrap())
                .unwrap();

        assert_eq!(transform.interior_seed(300.0).unwrap(), 800.0);
        assert_eq!(transform.interior_seed(1_300.0).unwrap(), 800.0);
        assert_eq!(transform.interior_seed(700.0).unwrap(), 700.0);
        assert!(transform.interior_seed(299.0).is_err());
    }

    fn prepared_two_species_ph_with_options(
        include_all_cp: bool,
        failing_gibbs: bool,
    ) -> Result<PreparedPhFormulation, ReactionExtentError> {
        let initial_moles = vec![0.75, 0.5];
        let problem = EquilibriumProblem::new(
            vec!["A2".to_string(), "A".to_string()],
            initial_moles.clone(),
            LogMolesInitialGuess::from_initial_moles(&initial_moles).unwrap(),
            DMatrix::from_row_slice(2, 1, &[2.0, 1.0]),
            vec![
                Rc::new(|temperature| 10_000.0 - 20.0 * temperature) as GibbsFn,
                Rc::new(|temperature| -5_000.0 - 5.0 * temperature) as GibbsFn,
            ],
            vec![Phase {
                kind: PhaseKind::IdealGas,
                species: vec![0, 1],
            }],
            EquilibriumConditions::new(1_000.0, 101_325.0, 101_325.0).unwrap(),
        )
        .unwrap();
        let prepared = PreparedEquilibriumProblem::new(problem).unwrap();
        let phase = PhaseId::new(None);
        let provenance = ["A2", "A"]
            .into_iter()
            .map(|substance| {
                ThermochemistryProvenance::new(
                    PhaseComponentId::new(phase.clone(), substance),
                    "synthetic",
                    substance,
                    "gas",
                )
            })
            .collect();
        let gibbs: Vec<MolarThermoFunction> = if failing_gibbs {
            vec![
                Arc::new(|_| {
                    Err(ReactionExtentError::InvalidProblem {
                        field: "fixture_gibbs",
                        message: "synthetic Gibbs capability failed".to_string(),
                    })
                }),
                Arc::new(|temperature| Ok(-5_000.0 - 5.0 * temperature)),
            ]
        } else {
            vec![
                Arc::new(|temperature| Ok(10_000.0 - 20.0 * temperature)),
                Arc::new(|temperature| Ok(-5_000.0 - 5.0 * temperature)),
            ]
        };
        let enthalpy: Vec<MolarThermoFunction> =
            vec![Arc::new(|_| Ok(10_000.0)), Arc::new(|_| Ok(-5_000.0))];
        let heat_capacity: Vec<Option<MolarThermoFunction>> = vec![
            include_all_cp.then(|| Arc::new(|_| Ok(0.0)) as MolarThermoFunction),
            Some(Arc::new(|_| Ok(0.0))),
        ];
        let thermochemistry = ResolvedThermochemistry::from_functions(
            provenance,
            TemperatureBounds::new(300.0, 3_000.0).unwrap(),
            gibbs,
            enthalpy,
            heat_capacity,
        )
        .unwrap();
        PreparedPhFormulation::new(
            prepared,
            thermochemistry,
            TemperatureBounds::new(500.0, 2_500.0).unwrap(),
            5_000.0,
            EnthalpyScale::new(10_000.0).unwrap(),
        )
    }

    fn prepared_two_species_ph_with_cp(include_all_cp: bool) -> PreparedPhFormulation {
        prepared_two_species_ph_with_options(include_all_cp, false).unwrap()
    }

    fn prepared_two_species_ph() -> PreparedPhFormulation {
        prepared_two_species_ph_with_cp(true)
    }

    #[test]
    fn monolithic_ph_residual_and_jacobian_have_the_declared_square_shape() {
        let formulation = prepared_two_species_ph();
        let unknowns = formulation.initial_unknowns(1_000.0).unwrap();

        assert_eq!(formulation.residual(&unknowns).unwrap().len(), 3);
        let jacobian = formulation.jacobian(&unknowns).unwrap();
        assert_eq!(jacobian.shape(), (3, 3));
    }

    #[test]
    fn monolithic_ph_analytic_jacobian_matches_central_differences() {
        let formulation = prepared_two_species_ph();
        let unknowns = formulation.initial_unknowns(1_000.0).unwrap();
        let analytic = formulation.jacobian(&unknowns).unwrap();
        let step = 1.0e-6;

        for column in 0..unknowns.len() {
            let mut plus = unknowns.clone();
            let mut minus = unknowns.clone();
            plus[column] += step;
            minus[column] -= step;
            let residual_plus = formulation.residual(&plus).unwrap();
            let residual_minus = formulation.residual(&minus).unwrap();
            for row in 0..unknowns.len() {
                let finite_difference = (residual_plus[row] - residual_minus[row]) / (2.0 * step);
                assert!(
                    (analytic[(row, column)] - finite_difference).abs() < 1.0e-6,
                    "P,H Jacobian mismatch at ({row}, {column}): analytic={}, finite_difference={finite_difference}",
                    analytic[(row, column)]
                );
            }
        }
    }

    #[test]
    fn monolithic_ph_rejects_missing_heat_capacity_during_preparation() {
        let error = match prepared_two_species_ph_with_options(false, false) {
            Ok(_) => panic!("missing Cp must be rejected before backend execution"),
            Err(error) => error,
        };
        assert!(matches!(
            error,
            ReactionExtentError::InvalidProblem {
                field: "heat_capacity",
                ..
            }
        ));
    }

    #[test]
    fn monolithic_ph_propagates_gibbs_capability_errors_without_nan() {
        let error = match prepared_two_species_ph_with_options(true, true) {
            Ok(_) => panic!("a failing Gibbs capability must reject preparation"),
            Err(error) => error,
        };
        assert!(matches!(
            error,
            ReactionExtentError::InvalidProblem {
                field: "fixture_gibbs",
                ..
            }
        ));
    }
}
