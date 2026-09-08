//! Independent fixed-pressure, fixed-enthalpy validation for one pure phase.
//!
//! This module is the P,H counterpart of `phase_boundary_validation`, but it
//! deliberately does not broaden that fixed-P,T problem type.  It solves a
//! narrow, independently supplied reaction through two scalar equations:
//!
//! ```text
//! ln(Q(xi)) - ln(K(T)) = 0
//! H(xi, T) - H_target = 0
//! ```
//!
//! The implementation owns neither canonical P,H residuals nor a nonlinear
//! backend policy.  Its purpose is validation by a second formulation, not a
//! replacement for the production multiphase workflow.

use nalgebra::DMatrix;

use crate::Thermodynamics::ChemEquilibrium::equilibrium_nonlinear::{
    ReactionExtentError, SolveError,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_ph_thermochemistry::MolarThermoFunction;
use crate::Thermodynamics::ChemEquilibrium::phase_boundary_validation::{
    MOLAR_GAS_CONSTANT, PurePhaseBoundaryElementComposition, PurePhaseBoundaryReactionSpace,
    PurePhaseBoundaryStructuralTolerances, independent_matrix_rank,
};

/// Fixed pressure, target enthalpy, and bounded temperature search interval
/// for the independent P,H problem.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PurePhasePhConditions {
    pressure: f64,
    reference_pressure: f64,
    target_enthalpy: f64,
    lower_temperature: f64,
    upper_temperature: f64,
}

impl PurePhasePhConditions {
    /// Builds finite physical conditions and an ordered positive temperature
    /// search interval. A P,H validator refuses to invent an unbounded search
    /// because that would hide a fixture or modelling error.
    pub fn new(
        pressure: f64,
        reference_pressure: f64,
        target_enthalpy: f64,
        lower_temperature: f64,
        upper_temperature: f64,
    ) -> Result<Self, ReactionExtentError> {
        for (field, value) in [
            ("pressure", pressure),
            ("reference_pressure", reference_pressure),
        ] {
            if !value.is_finite() || value <= 0.0 {
                return Err(ReactionExtentError::InvalidConditions {
                    parameter: field,
                    value,
                });
            }
        }
        if !target_enthalpy.is_finite() {
            return Err(invalid_problem("target enthalpy must be finite"));
        }
        if !lower_temperature.is_finite()
            || !upper_temperature.is_finite()
            || lower_temperature <= 0.0
            || upper_temperature <= lower_temperature
        {
            return Err(invalid_problem(
                "temperature bracket must be finite, positive, and ordered",
            ));
        }
        Ok(Self {
            pressure,
            reference_pressure,
            target_enthalpy,
            lower_temperature,
            upper_temperature,
        })
    }

    pub fn pressure(self) -> f64 {
        self.pressure
    }

    pub fn reference_pressure(self) -> f64 {
        self.reference_pressure
    }

    pub fn target_enthalpy(self) -> f64 {
        self.target_enthalpy
    }

    pub fn lower_temperature(self) -> f64 {
        self.lower_temperature
    }

    pub fn upper_temperature(self) -> f64 {
        self.upper_temperature
    }
}

/// Immutable identity of one independent P,H validation case.
///
/// The identity deliberately contains the physical inputs that distinguish
/// otherwise similarly named problems: inventory, reaction normalization,
/// `P/P0`, target enthalpy, and the declared temperature bracket. It does not
/// attempt to compare function pointers for thermochemistry; fixture
/// provenance owns that separate responsibility.
#[derive(Debug, Clone, PartialEq)]
pub struct PurePhasePhCaseIdentity {
    gas_species: Vec<String>,
    initial_gas_moles: Vec<f64>,
    gas_stoichiometry: Vec<f64>,
    initial_candidate_moles: f64,
    candidate_stoichiometry: f64,
    candidate_name: String,
    conditions: PurePhasePhConditions,
}

/// Narrow independent P,H problem for one phase-forming reaction
///
/// ```text
/// sum_i nu_i A_i(g) + nu_s S(condensed) = 0,  nu_s > 0.
/// ```
///
/// The gas phase is ideal and the pure condensed candidate has unit activity.
/// `initial_candidate_moles` is physical inventory, so this type can later
/// express both an initially absent and an initially present candidate without
/// changing the reaction-coordinate contract.
#[derive(Clone)]
pub struct PurePhasePhProblem {
    gas_species: Vec<String>,
    initial_gas_moles: Vec<f64>,
    gas_stoichiometry: Vec<f64>,
    initial_candidate_moles: f64,
    candidate_stoichiometry: f64,
    candidate_name: String,
    gas_standard_gibbs: Vec<MolarThermoFunction>,
    candidate_standard_gibbs: MolarThermoFunction,
    gas_molar_enthalpies: Vec<MolarThermoFunction>,
    candidate_molar_enthalpy: MolarThermoFunction,
    conditions: PurePhasePhConditions,
    element_composition: Option<PurePhaseBoundaryElementComposition>,
}

impl PurePhasePhProblem {
    /// Constructs a fully typed independent P,H validation problem.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        gas_species: Vec<String>,
        initial_gas_moles: Vec<f64>,
        gas_stoichiometry: Vec<f64>,
        initial_candidate_moles: f64,
        candidate_stoichiometry: f64,
        candidate_name: impl Into<String>,
        gas_standard_gibbs: Vec<MolarThermoFunction>,
        candidate_standard_gibbs: MolarThermoFunction,
        gas_molar_enthalpies: Vec<MolarThermoFunction>,
        candidate_molar_enthalpy: MolarThermoFunction,
        conditions: PurePhasePhConditions,
    ) -> Result<Self, ReactionExtentError> {
        let gas_count = gas_species.len();
        if gas_count == 0 {
            return Err(invalid_problem("at least one gas species is required"));
        }
        if gas_species.iter().any(|name| name.trim().is_empty()) {
            return Err(invalid_problem("gas species identities must be non-empty"));
        }
        let mut names = std::collections::BTreeSet::new();
        if gas_species.iter().any(|name| !names.insert(name.as_str())) {
            return Err(invalid_problem("gas species identities must be unique"));
        }
        if initial_gas_moles.len() != gas_count
            || gas_stoichiometry.len() != gas_count
            || gas_standard_gibbs.len() != gas_count
            || gas_molar_enthalpies.len() != gas_count
        {
            return Err(ReactionExtentError::DimensionMismatch(format!(
                "P,H pure-phase problem has {gas_count} gas species but moles/stoichiometry/Gibbs/enthalpy lengths are {}/{}/{}/{}",
                initial_gas_moles.len(),
                gas_stoichiometry.len(),
                gas_standard_gibbs.len(),
                gas_molar_enthalpies.len(),
            )));
        }
        if initial_gas_moles
            .iter()
            .any(|moles| !moles.is_finite() || *moles <= 0.0)
        {
            return Err(invalid_problem(
                "initial gas mole numbers must be finite and strictly positive",
            ));
        }
        if gas_stoichiometry.iter().any(|nu| !nu.is_finite())
            || gas_stoichiometry.iter().all(|nu| *nu == 0.0)
        {
            return Err(invalid_problem(
                "gas stoichiometry must be finite and include a non-zero coefficient",
            ));
        }
        if !initial_candidate_moles.is_finite() || initial_candidate_moles < 0.0 {
            return Err(invalid_problem(
                "initial candidate moles must be finite and non-negative",
            ));
        }
        if !candidate_stoichiometry.is_finite() || candidate_stoichiometry <= 0.0 {
            return Err(invalid_problem(
                "candidate stoichiometric coefficient must be finite and strictly positive",
            ));
        }
        let candidate_name = candidate_name.into();
        if candidate_name.trim().is_empty() {
            return Err(invalid_problem("candidate identity must be non-empty"));
        }

        Ok(Self {
            gas_species,
            initial_gas_moles,
            gas_stoichiometry,
            initial_candidate_moles,
            candidate_stoichiometry,
            candidate_name,
            gas_standard_gibbs,
            candidate_standard_gibbs,
            gas_molar_enthalpies,
            candidate_molar_enthalpy,
            conditions,
            element_composition: None,
        })
    }

    /// Returns the same independent physical problem with a new enthalpy
    /// target.
    ///
    /// This is intentionally a value-returning operation: a validation sweep
    /// must not mutate a previously accepted reference problem while it is
    /// being compared with a production continuation point.
    pub fn with_target_enthalpy(&self, target_enthalpy: f64) -> Result<Self, ReactionExtentError> {
        let mut retargeted = self.clone();
        retargeted.conditions = PurePhasePhConditions::new(
            self.conditions.pressure,
            self.conditions.reference_pressure,
            target_enthalpy,
            self.conditions.lower_temperature,
            self.conditions.upper_temperature,
        )?;
        Ok(retargeted)
    }

    /// Attaches independent elemental data.  The reaction is rejected rather
    /// than projected when it violates a declared elemental balance.
    pub fn with_element_composition(
        mut self,
        composition: PurePhaseBoundaryElementComposition,
        tolerances: PurePhaseBoundaryStructuralTolerances,
    ) -> Result<Self, ReactionExtentError> {
        validate_structural_tolerances(tolerances)?;
        if composition.gas_species_by_element().nrows() != self.gas_species.len() {
            return Err(invalid_problem(format!(
                "element composition has {} gas rows for {} gas species",
                composition.gas_species_by_element().nrows(),
                self.gas_species.len()
            )));
        }
        let residuals = self.element_balance_residuals_for(&composition)?;
        if residuals
            .iter()
            .any(|residual| residual.abs() > tolerances.max_abs_element_balance)
        {
            return Err(invalid_problem(format!(
                "phase-forming reaction violates elemental conservation: residuals={residuals:?}"
            )));
        }
        self.element_composition = Some(composition);
        Ok(self)
    }

    pub fn gas_species(&self) -> &[String] {
        &self.gas_species
    }

    pub fn candidate_name(&self) -> &str {
        &self.candidate_name
    }

    pub fn conditions(&self) -> PurePhasePhConditions {
        self.conditions
    }

    /// Returns the validated physical identity carried by every independent
    /// result and canonical comparison record for this problem.
    pub fn case_identity(&self) -> PurePhasePhCaseIdentity {
        PurePhasePhCaseIdentity {
            gas_species: self.gas_species.clone(),
            initial_gas_moles: self.initial_gas_moles.clone(),
            gas_stoichiometry: self.gas_stoichiometry.clone(),
            initial_candidate_moles: self.initial_candidate_moles,
            candidate_stoichiometry: self.candidate_stoichiometry,
            candidate_name: self.candidate_name.clone(),
            conditions: self.conditions,
        }
    }

    pub fn gas_stoichiometry(&self) -> &[f64] {
        &self.gas_stoichiometry
    }

    pub fn candidate_stoichiometry(&self) -> f64 {
        self.candidate_stoichiometry
    }

    /// Returns independent reaction-space evidence when composition was
    /// supplied.  This uses only the P8 independent SVD helper.
    pub fn reaction_space(
        &self,
        tolerances: PurePhaseBoundaryStructuralTolerances,
    ) -> Result<PurePhaseBoundaryReactionSpace, ReactionExtentError> {
        validate_structural_tolerances(tolerances)?;
        let composition = self.element_composition.as_ref().ok_or_else(|| {
            ReactionExtentError::ValidationNotApplicable {
                path: "pure_phase_ph_structure",
                message: "strict reaction-space validation requires element composition".into(),
            }
        })?;
        let gas_only_rank =
            independent_matrix_rank(composition.gas_species_by_element(), tolerances)?;
        let mut full = DMatrix::zeros(
            composition.gas_species_by_element().nrows() + 1,
            composition.gas_species_by_element().ncols(),
        );
        full.rows_mut(0, composition.gas_species_by_element().nrows())
            .copy_from(composition.gas_species_by_element());
        for (column, value) in composition.candidate_by_element().iter().enumerate() {
            full[(composition.gas_species_by_element().nrows(), column)] = *value;
        }
        let full_rank = independent_matrix_rank(&full, tolerances)?;
        Ok(PurePhaseBoundaryReactionSpace {
            full_rank,
            gas_only_rank,
            full_reaction_dimension: full.nrows().saturating_sub(full_rank),
            gas_only_reaction_dimension: composition
                .gas_species_by_element()
                .nrows()
                .saturating_sub(gas_only_rank),
            element_balance_residuals: self.element_balance_residuals_for(composition)?,
        })
    }

    /// Enforces the P9 strict family: one full reaction, no gas-only reaction.
    pub fn validate_strict_independent_family(
        &self,
        tolerances: PurePhaseBoundaryStructuralTolerances,
    ) -> Result<PurePhaseBoundaryReactionSpace, ReactionExtentError> {
        validate_structural_tolerances(tolerances)?;
        let evidence = self.reaction_space(tolerances)?;
        if evidence
            .element_balance_residuals
            .iter()
            .any(|residual| residual.abs() > tolerances.max_abs_element_balance)
        {
            return Err(invalid_problem(
                "strict family violates elemental conservation",
            ));
        }
        if evidence.full_reaction_dimension != 1 {
            return Err(invalid_problem(format!(
                "strict family requires one full-system reaction dimension, got {}",
                evidence.full_reaction_dimension
            )));
        }
        if evidence.gas_only_reaction_dimension != 0 {
            return Err(invalid_problem(format!(
                "strict family requires no gas-only reaction dimension, got {}",
                evidence.gas_only_reaction_dimension
            )));
        }
        Ok(evidence)
    }

    /// Solves the inner chemical equation at a prescribed positive temperature.
    pub fn solve_inner_extent_at_temperature(
        &self,
        temperature: f64,
        settings: PurePhasePhSolverSettings,
    ) -> Result<PurePhasePhInnerExtentResult, ReactionExtentError> {
        settings.validate()?;
        self.solve_inner_extent(temperature, settings)
    }

    /// Materializes one feasible independent reaction-coordinate state.
    ///
    /// This deliberately stays below the nested `P,H` solver: callers that
    /// already own an independently accepted fixed-temperature extent can
    /// construct a physical reference state and its additive enthalpy without
    /// borrowing canonical residuals or duplicating closure arithmetic.
    pub fn materialize_state_at_extent(
        &self,
        extent: f64,
        temperature: f64,
    ) -> Result<PurePhasePhMaterializedState, ReactionExtentError> {
        let (gas_moles, candidate_moles, total_enthalpy) =
            self.total_enthalpy_at(extent, temperature)?;
        let chemical_log_residual = self.log_residual_for_gas_moles(&gas_moles, temperature)?;
        let max_abs_element_balance =
            self.accepted_state_max_abs_element_balance(&gas_moles, candidate_moles)?;
        Ok(PurePhasePhMaterializedState {
            extent,
            temperature,
            gas_moles,
            candidate_moles,
            total_enthalpy,
            chemical_log_residual,
            max_abs_element_balance,
        })
    }

    /// Re-evaluates the independent chemical equation for an externally
    /// accepted gas composition. P9.4 uses this for a canonical candidate
    /// rather than assigning physical meaning to heterogeneous solver rows.
    pub(crate) fn chemical_log_residual_for_gas_moles(
        &self,
        gas_moles: &[f64],
        temperature: f64,
    ) -> Result<f64, ReactionExtentError> {
        if gas_moles.len() != self.gas_species.len() {
            return Err(ReactionExtentError::DimensionMismatch(format!(
                "P,H pure-phase chemical check has {} gas values for {} species",
                gas_moles.len(),
                self.gas_species.len()
            )));
        }
        if !temperature.is_finite() || temperature <= 0.0 {
            return Err(ReactionExtentError::InvalidConditions {
                parameter: "temperature",
                value: temperature,
            });
        }
        self.log_residual_for_gas_moles(gas_moles, temperature)
    }

    /// Evaluates the per-element conservation residual of the supplied
    /// elemental composition under this problem's fixed reaction.
    ///
    /// For each element column, the residual is the sum over gas species of
    /// `nu_i * A[i, e]` plus the candidate contribution
    /// `nu_s * A_s[e]`. A non-finite result indicates an inconsistent or
    /// unphysical composition and is rejected as an invalid problem.
    fn element_balance_residuals_for(
        &self,
        composition: &PurePhaseBoundaryElementComposition,
    ) -> Result<Vec<f64>, ReactionExtentError> {
        let columns = composition.gas_species_by_element().ncols();
        let mut residuals = vec![0.0; columns];
        for (gas_index, nu) in self.gas_stoichiometry.iter().enumerate() {
            for element_index in 0..columns {
                residuals[element_index] +=
                    nu * composition.gas_species_by_element()[(gas_index, element_index)];
            }
        }
        for (element_index, amount) in composition.candidate_by_element().iter().enumerate() {
            residuals[element_index] += self.candidate_stoichiometry * amount;
        }
        if residuals.iter().any(|value| !value.is_finite()) {
            return Err(invalid_problem("element-balance residual is non-finite"));
        }
        Ok(residuals)
    }

    /// Computes the accepted-state elemental drift relative to the initial
    /// physical inventory. `None` means that this independent problem was
    /// intentionally created without elemental data and therefore cannot make
    /// a conservation claim.
    fn accepted_state_max_abs_element_balance(
        &self,
        gas_moles: &[f64],
        candidate_moles: f64,
    ) -> Result<Option<f64>, ReactionExtentError> {
        let Some(composition) = self.element_composition.as_ref() else {
            return Ok(None);
        };
        if gas_moles.len() != self.gas_species.len() {
            return Err(ReactionExtentError::DimensionMismatch(format!(
                "accepted P,H state has {} gas moles for {} species",
                gas_moles.len(),
                self.gas_species.len()
            )));
        }
        if gas_moles
            .iter()
            .any(|moles| !moles.is_finite() || *moles <= 0.0)
            || !candidate_moles.is_finite()
            || candidate_moles <= 0.0
        {
            return Err(invalid_candidate(
                "accepted P,H state has non-finite or non-positive physical moles",
            ));
        }

        let mut max_abs_balance = 0.0_f64;
        for element in 0..composition.gas_species_by_element().ncols() {
            let gas_delta = gas_moles
                .iter()
                .zip(&self.initial_gas_moles)
                .enumerate()
                .map(|(species, (accepted, initial))| {
                    (accepted - initial) * composition.gas_species_by_element()[(species, element)]
                })
                .sum::<f64>();
            let candidate_delta = (candidate_moles - self.initial_candidate_moles)
                * composition.candidate_by_element()[element];
            let balance = gas_delta + candidate_delta;
            if !balance.is_finite() {
                return Err(invalid_problem(
                    "accepted P,H element-balance residual is non-finite",
                ));
            }
            max_abs_balance = max_abs_balance.max(balance.abs());
        }
        Ok(Some(max_abs_balance))
    }

    /// Computes the open interval of reaction extents for which every species
    /// (gas and candidate) remains strictly positive.
    ///
    /// Each positive stoichiometric coefficient imposes a lower bound
    /// `-n_i/nu_i` and each negative one an upper bound, intersected across
    /// all species. The interval is then inset by a relative `margin` so the
    /// returned bracketing points are strictly interior and avoid degenerate
    /// zero-mole endpoints. A collapsed or empty interval is reported as
    /// [`ReactionExtentError::ValidationNotApplicable`].
    fn feasible_interior_extent_interval(
        &self,
        margin: f64,
    ) -> Result<(f64, f64), ReactionExtentError> {
        let mut lower = f64::NEG_INFINITY;
        let mut upper = f64::INFINITY;
        for (&moles, &nu) in self
            .initial_gas_moles
            .iter()
            .zip(&self.gas_stoichiometry)
            .chain(std::iter::once((
                &self.initial_candidate_moles,
                &self.candidate_stoichiometry,
            )))
        {
            if nu > 0.0 {
                lower = lower.max(-moles / nu);
            } else if nu < 0.0 {
                upper = upper.min(-moles / nu);
            }
        }
        if !lower.is_finite() || !upper.is_finite() || upper <= lower {
            return Err(ReactionExtentError::ValidationNotApplicable {
                path: "pure_phase_ph_inner_extent",
                message: "reaction has no finite non-empty physical extent interval".into(),
            });
        }
        let inset = margin.max((upper - lower).abs() * margin);
        let left = lower + inset;
        let right = upper - inset;
        if !left.is_finite() || !right.is_finite() || right <= left {
            return Err(ReactionExtentError::ValidationNotApplicable {
                path: "pure_phase_ph_inner_extent",
                message: "physical extent interval collapsed after positivity margin".into(),
            });
        }
        Ok((left, right))
    }

    /// Computes physical gas mole numbers at reaction extent `extent`.
    ///
    /// Each gas species follows `n_i(extent) = n_i(0) + nu_i * extent`. Every
    /// result must be finite and strictly positive; a violation is reported as
    /// an invalid candidate because the log-mole/activity formulation is
    /// undefined for non-positive moles.
    fn gas_moles_at_extent(&self, extent: f64) -> Result<Vec<f64>, ReactionExtentError> {
        if !extent.is_finite() {
            return Err(invalid_candidate("reaction extent must be finite"));
        }
        self.initial_gas_moles
            .iter()
            .zip(&self.gas_stoichiometry)
            .enumerate()
            .map(|(index, (&initial, &nu))| {
                let moles = initial + nu * extent;
                if !moles.is_finite() || moles <= 0.0 {
                    return Err(ReactionExtentError::InvalidCandidate {
                        field: "pure_phase_ph_gas_moles",
                        message: format!("gas species {index} is non-positive at extent {extent}"),
                    });
                }
                Ok(moles)
            })
            .collect()
    }

    /// Computes the physical mole number of the pure candidate phase at extent.
    ///
    /// Follows the same linear extent law as the gas species. The result must
    /// remain finite and strictly positive because the candidate is always
    /// treated as an active, present phase in this fixed-topology formulation.
    fn candidate_moles_at_extent(&self, extent: f64) -> Result<f64, ReactionExtentError> {
        let moles = self.initial_candidate_moles + self.candidate_stoichiometry * extent;
        if !moles.is_finite() || moles <= 0.0 {
            return Err(invalid_candidate(
                "pure candidate moles must be finite and strictly positive in fixed topology",
            ));
        }
        Ok(moles)
    }

    /// Evaluates the chemical `ln(Q) - ln(K)` residual at a reaction extent.
    ///
    /// Convenience wrapper that first materializes the gas moles implied by
    /// `extent`, then delegates to [`log_residual_for_gas_moles`] at the given
    /// temperature.
    fn log_residual_at(&self, extent: f64, temperature: f64) -> Result<f64, ReactionExtentError> {
        let gas_moles = self.gas_moles_at_extent(extent)?;
        self.log_residual_for_gas_moles(&gas_moles, temperature)
    }

    /// Evaluates the independent chemical residual `ln(Q) - ln(K)` for an
    /// explicit gas composition at one temperature.
    ///
    /// The reaction quotient `ln(Q)` is built from ideal-gas activities
    /// `a_i = (n_i/N_gas) * (P/P0)`, and `ln(K) = -DeltaG°(T)/(R*T)` uses the
    /// standard-state Gibbs closures of the gas species and the pure candidate.
    /// A zero residual marks the phase-boundary condition. Non-finite
    /// activities, Gibbs values, or the final residual are rejected with typed
    /// errors rather than propagated as `NaN`.
    fn log_residual_for_gas_moles(
        &self,
        gas_moles: &[f64],
        temperature: f64,
    ) -> Result<f64, ReactionExtentError> {
        if gas_moles
            .iter()
            .any(|moles| !moles.is_finite() || *moles <= 0.0)
        {
            return Err(invalid_candidate(
                "chemical log residual requires finite positive gas moles",
            ));
        }
        let gas_total: f64 = gas_moles.iter().sum();
        if !gas_total.is_finite() || gas_total <= 0.0 {
            return Err(invalid_candidate(
                "total gas moles must be finite and positive",
            ));
        }
        let pressure_ratio = self.conditions.pressure / self.conditions.reference_pressure;
        let mut ln_q = 0.0;
        for (&nu, &moles) in self.gas_stoichiometry.iter().zip(gas_moles) {
            if nu != 0.0 {
                let activity = moles / gas_total * pressure_ratio;
                if !activity.is_finite() || activity <= 0.0 {
                    return Err(invalid_candidate(
                        "ideal-gas activity must be finite and positive",
                    ));
                }
                ln_q += nu * activity.ln();
            }
        }
        let mut delta_g = 0.0;
        for (index, (&nu, function)) in self
            .gas_stoichiometry
            .iter()
            .zip(&self.gas_standard_gibbs)
            .enumerate()
        {
            let value = function(temperature)?;
            if !value.is_finite() {
                return Err(ReactionExtentError::InvalidDG0 {
                    species_index: index,
                    dg0: value,
                    temperature,
                });
            }
            delta_g += nu * value;
        }
        let candidate_gibbs = (self.candidate_standard_gibbs)(temperature)?;
        if !candidate_gibbs.is_finite() {
            return Err(ReactionExtentError::InvalidDG0 {
                species_index: self.gas_species.len(),
                dg0: candidate_gibbs,
                temperature,
            });
        }
        delta_g += self.candidate_stoichiometry * candidate_gibbs;
        let ln_k = -delta_g / (MOLAR_GAS_CONSTANT * temperature);
        let residual = ln_q - ln_k;
        if !residual.is_finite() {
            return Err(ReactionExtentError::ResidualEvaluation(
                "pure-phase P,H chemical log residual is non-finite".into(),
            ));
        }
        Ok(residual)
    }

    /// Evaluates the additive total enthalpy of gas and candidate at a given
    /// extent and temperature.
    ///
    /// Computes `H = sum_i n_i * h_i(T) + n_s * h_s(T)` from the molar-enthalpy
    /// closures. Returns the materialized gas moles, candidate moles, and the
    /// total enthalpy so callers can reconstruct the full state without a
    /// second evaluation. Any non-finite enthalpy contribution aborts with a
    /// typed invalid-candidate error.
    fn total_enthalpy_at(
        &self,
        extent: f64,
        temperature: f64,
    ) -> Result<(Vec<f64>, f64, f64), ReactionExtentError> {
        let gas_moles = self.gas_moles_at_extent(extent)?;
        let candidate_moles = self.candidate_moles_at_extent(extent)?;
        let mut total = 0.0;
        for (index, (&moles, function)) in
            gas_moles.iter().zip(&self.gas_molar_enthalpies).enumerate()
        {
            let enthalpy = function(temperature)?;
            if !enthalpy.is_finite() {
                return Err(ReactionExtentError::InvalidCandidate {
                    field: "pure_phase_ph_enthalpy",
                    message: format!(
                        "gas enthalpy function {index} returned {enthalpy} at {temperature} K"
                    ),
                });
            }
            total += moles * enthalpy;
        }
        let candidate_enthalpy = (self.candidate_molar_enthalpy)(temperature)?;
        if !candidate_enthalpy.is_finite() {
            return Err(ReactionExtentError::InvalidCandidate {
                field: "pure_phase_ph_enthalpy",
                message: format!(
                    "candidate enthalpy function returned {candidate_enthalpy} at {temperature} K"
                ),
            });
        }
        total += candidate_moles * candidate_enthalpy;
        if !total.is_finite() {
            return Err(ReactionExtentError::ResidualEvaluation(
                "pure-phase P,H total enthalpy is non-finite".into(),
            ));
        }
        Ok((gas_moles, candidate_moles, total))
    }

    /// Solves the inner chemical equation `ln(Q) - ln(K) = 0` for the reaction
    /// extent at a fixed temperature using safeguarded bisection.
    ///
    /// The root is searched inside the feasible interior extent interval
    /// derived from positivity constraints. If the endpoint residuals already
    /// satisfy the tolerance the corresponding endpoint is returned with zero
    /// iterations; otherwise a sign change is required to bracket the root and
    /// bisection proceeds up to `max_inner_iterations`. An un-bracketed or
    /// non-converged search is reported as a validation-not-applicable error or
    /// a max-iterations solve error respectively.
    fn solve_inner_extent(
        &self,
        temperature: f64,
        settings: PurePhasePhSolverSettings,
    ) -> Result<PurePhasePhInnerExtentResult, ReactionExtentError> {
        if !temperature.is_finite() || temperature <= 0.0 {
            return Err(ReactionExtentError::InvalidConditions {
                parameter: "temperature",
                value: temperature,
            });
        }
        let (mut left, mut right) =
            self.feasible_interior_extent_interval(settings.feasibility_margin)?;
        let mut f_left = self.log_residual_at(left, temperature)?;
        if f_left.abs() <= settings.max_abs_log_residual {
            return Ok(PurePhasePhInnerExtentResult {
                extent: left,
                log_residual: f_left,
                iterations: 0,
            });
        }
        let f_right = self.log_residual_at(right, temperature)?;
        if f_right.abs() <= settings.max_abs_log_residual {
            return Ok(PurePhasePhInnerExtentResult {
                extent: right,
                log_residual: f_right,
                iterations: 0,
            });
        }
        if f_left.signum() == f_right.signum() {
            return Err(ReactionExtentError::ValidationNotApplicable {
                path: "pure_phase_ph_inner_extent",
                message: format!(
                    "chemical extent root is not bracketed at {temperature} K: endpoint residuals are {f_left:e} and {f_right:e}"
                ),
            });
        }
        for iteration in 1..=settings.max_inner_iterations {
            let mid = 0.5 * (left + right);
            let f_mid = self.log_residual_at(mid, temperature)?;
            if f_mid.abs() <= settings.max_abs_log_residual {
                return Ok(PurePhasePhInnerExtentResult {
                    extent: mid,
                    log_residual: f_mid,
                    iterations: iteration,
                });
            }
            if f_left.signum() != f_mid.signum() {
                right = mid;
            } else {
                left = mid;
                f_left = f_mid;
            }
        }
        Err(ReactionExtentError::SolveError(SolveError::MaxIterations))
    }
}

/// Physical state reconstructed by the independent extent formulation.
///
/// It is intentionally a small value object: this is enough to derive a
/// reference `H_target` from independent local thermochemistry while keeping
/// solver iterations, phase-control state, and canonical caches out of the
/// validation boundary.
#[derive(Debug, Clone, PartialEq)]
pub struct PurePhasePhMaterializedState {
    /// Reaction coordinate used to materialize this state.
    pub extent: f64,
    /// Physical temperature in K.
    pub temperature: f64,
    /// Ordered physical gas mole amounts.
    pub gas_moles: Vec<f64>,
    /// Physical amount of the pure condensed candidate.
    pub candidate_moles: f64,
    /// Additive total enthalpy in J.
    pub total_enthalpy: f64,
    /// Independent `ln(Q)-ln(K)` residual at this physical state.
    pub chemical_log_residual: f64,
    /// Independently recomputed elemental-balance error when the problem has
    /// explicit elemental composition.
    pub max_abs_element_balance: Option<f64>,
}

/// Immutable canonical fixed-topology evidence supplied to the P9 comparator.
///
/// The canonical P,H runner owns its own acceptance report. This compact
/// record holds only physical data needed for an explicit comparison; it does
/// not expose mutable solver internals.
#[derive(Debug, Clone, PartialEq)]
pub struct PurePhasePhCanonicalEvidence {
    /// Physical P,H case this evidence claims to describe.
    pub case_identity: PurePhasePhCaseIdentity,
    /// Temperature of the canonical accepted state, in K.
    pub temperature: f64,
    /// Physical gas mole numbers at the accepted state.
    pub gas_moles: Vec<f64>,
    /// Physical mole number of the pure candidate at the accepted state.
    pub candidate_moles: f64,
    /// Additive total enthalpy of the accepted state, in J.
    pub total_enthalpy: f64,
    /// Independent `ln(Q)-ln(K)` evaluated at the canonical accepted state.
    pub chemical_log_residual: Option<f64>,
    /// Explicit elemental conservation error from canonical physical moles.
    pub max_abs_element_balance: Option<f64>,
}

impl PurePhasePhCanonicalEvidence {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        case_identity: PurePhasePhCaseIdentity,
        temperature: f64,
        gas_moles: Vec<f64>,
        candidate_moles: f64,
        total_enthalpy: f64,
        chemical_log_residual: Option<f64>,
        max_abs_element_balance: Option<f64>,
    ) -> Result<Self, ReactionExtentError> {
        if case_identity.gas_species.len() != gas_moles.len() {
            return Err(ReactionExtentError::DimensionMismatch(format!(
                "canonical P,H evidence has {} gas identities but {} gas moles",
                case_identity.gas_species.len(),
                gas_moles.len()
            )));
        }
        if !temperature.is_finite()
            || temperature <= 0.0
            || gas_moles
                .iter()
                .any(|moles| !moles.is_finite() || *moles <= 0.0)
            || !candidate_moles.is_finite()
            || candidate_moles <= 0.0
            || !total_enthalpy.is_finite()
        {
            return Err(invalid_problem(
                "canonical P,H evidence must carry finite positive physical state values",
            ));
        }
        for (field, value) in [
            ("canonical chemical log residual", chemical_log_residual),
            ("canonical element balance", max_abs_element_balance),
        ] {
            if value.is_some_and(|value| !value.is_finite() || value < 0.0) {
                return Err(invalid_problem(format!(
                    "{field} must be finite and non-negative"
                )));
            }
        }
        Ok(Self {
            case_identity,
            temperature,
            gas_moles,
            candidate_moles,
            total_enthalpy,
            chemical_log_residual,
            max_abs_element_balance,
        })
    }
}

/// Explicit tolerances for P9 fixed-topology comparison.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PurePhasePhCrossValidationTolerances {
    /// Maximum accepted absolute temperature difference, in K.
    pub max_abs_temperature_delta: f64,
    /// Relative temperature comparison tolerance.
    pub max_relative_temperature_delta: f64,
    /// Maximum accepted absolute mole-number difference per species.
    pub max_abs_mole_delta: f64,
    /// Relative mole-number comparison tolerance.
    pub max_relative_mole_delta: f64,
    /// Maximum accepted absolute total-enthalpy difference, in J.
    pub max_abs_enthalpy_delta: f64,
    /// Relative total-enthalpy comparison tolerance.
    pub max_relative_enthalpy_delta: f64,
    /// Maximum accepted absolute chemical log-residual magnitude.
    pub max_abs_chemical_log_residual: f64,
    /// Maximum accepted absolute elemental-balance error.
    pub max_abs_element_balance: f64,
}

impl Default for PurePhasePhCrossValidationTolerances {
    fn default() -> Self {
        Self {
            max_abs_temperature_delta: 1e-5,
            max_relative_temperature_delta: 1e-10,
            max_abs_mole_delta: 1e-7,
            max_relative_mole_delta: 1e-8,
            max_abs_enthalpy_delta: 1e-4,
            max_relative_enthalpy_delta: 1e-10,
            max_abs_chemical_log_residual: 1e-8,
            max_abs_element_balance: 1e-8,
        }
    }
}

impl PurePhasePhCrossValidationTolerances {
    /// Rejects any non-finite or non-positive comparison tolerance.
    ///
    /// Iterates over every typed axis and returns a [`ReactionExtentError`]
    /// naming the offending field, so cross-validation cannot silently compare
    /// two results with a disabled or degenerate tolerance.
    fn validate(self) -> Result<(), ReactionExtentError> {
        for (name, value) in [
            ("max_abs_temperature_delta", self.max_abs_temperature_delta),
            ("max_abs_mole_delta", self.max_abs_mole_delta),
            ("max_abs_enthalpy_delta", self.max_abs_enthalpy_delta),
            (
                "max_abs_chemical_log_residual",
                self.max_abs_chemical_log_residual,
            ),
            ("max_abs_element_balance", self.max_abs_element_balance),
        ] {
            if !value.is_finite() || value <= 0.0 {
                return Err(invalid_problem(format!(
                    "{name} must be finite and positive"
                )));
            }
        }
        for (name, value) in [
            (
                "max_relative_temperature_delta",
                self.max_relative_temperature_delta,
            ),
            ("max_relative_mole_delta", self.max_relative_mole_delta),
            (
                "max_relative_enthalpy_delta",
                self.max_relative_enthalpy_delta,
            ),
        ] {
            if !value.is_finite() || value < 0.0 {
                return Err(invalid_problem(format!(
                    "{name} must be finite and non-negative"
                )));
            }
        }
        Ok(())
    }
}

/// Axis-by-axis fixed-topology P,H comparison.
///
/// `None` means an axis is genuinely unavailable, never implicit success.
#[derive(Debug, Clone, PartialEq)]
pub struct PurePhasePhCrossValidationReport {
    /// Whether the canonical and independent problems identify the same case.
    pub identity_agreement: Option<bool>,
    /// Whether solved temperatures agree within tolerance.
    pub temperature_agreement: Option<bool>,
    /// Whether the chemical log residuals agree within tolerance.
    pub thermodynamic_agreement: Option<bool>,
    /// Whether gas/candidate mole numbers agree within tolerance.
    pub composition_agreement: Option<bool>,
    /// Whether total enthalpies agree within tolerance.
    pub enthalpy_agreement: Option<bool>,
    /// Whether element-balance errors agree within tolerance.
    pub conservation_agreement: Option<bool>,
    /// Independent accepted-state elemental balance, where element data exist.
    pub independent_max_abs_element_balance: Option<f64>,
    /// Canonical accepted-state elemental balance, when it was published.
    pub canonical_max_abs_element_balance: Option<f64>,
}

impl PurePhasePhCrossValidationReport {
    /// All P9.4 comparison axes were observed and agree.
    pub fn is_complete_match(&self) -> bool {
        [
            self.identity_agreement,
            self.temperature_agreement,
            self.thermodynamic_agreement,
            self.composition_agreement,
            self.enthalpy_agreement,
            self.conservation_agreement,
        ]
        .into_iter()
        .all(|axis| axis == Some(true))
    }
}

/// Compares independent P,H evidence with an accepted fixed-topology canonical
/// candidate. It deliberately has no phase-control or nonlinear-backend API.
pub fn compare_pure_phase_ph_validation(
    problem: &PurePhasePhProblem,
    independent: &PurePhasePhEquilibriumResult,
    canonical: &PurePhasePhCanonicalEvidence,
    tolerances: PurePhasePhCrossValidationTolerances,
) -> Result<PurePhasePhCrossValidationReport, ReactionExtentError> {
    tolerances.validate()?;
    let expected_identity = problem.case_identity();
    if independent.case_identity != expected_identity {
        return Err(invalid_problem(
            "independent P,H result belongs to a different validation case",
        ));
    }
    if canonical.case_identity != expected_identity {
        return Err(invalid_problem(
            "canonical P,H evidence belongs to a different validation case",
        ));
    }
    let temperature_agreement = Some(within_abs_relative(
        independent.temperature,
        canonical.temperature,
        tolerances.max_abs_temperature_delta,
        tolerances.max_relative_temperature_delta,
    ));
    let composition_agreement = Some(
        independent.gas_moles.len() == canonical.gas_moles.len()
            && independent
                .gas_moles
                .iter()
                .zip(&canonical.gas_moles)
                .all(|(left, right)| {
                    within_abs_relative(
                        *left,
                        *right,
                        tolerances.max_abs_mole_delta,
                        tolerances.max_relative_mole_delta,
                    )
                })
            && within_abs_relative(
                independent.candidate_moles,
                canonical.candidate_moles,
                tolerances.max_abs_mole_delta,
                tolerances.max_relative_mole_delta,
            ),
    );
    let enthalpy_agreement = Some(
        within_abs_relative(
            independent.total_enthalpy,
            canonical.total_enthalpy,
            tolerances.max_abs_enthalpy_delta,
            tolerances.max_relative_enthalpy_delta,
        ) && within_abs_relative(
            independent.enthalpy_residual,
            0.0,
            tolerances.max_abs_enthalpy_delta,
            tolerances.max_relative_enthalpy_delta,
        ),
    );
    let thermodynamic_agreement = canonical.chemical_log_residual.map(|residual| {
        independent.chemical_log_residual.abs() <= tolerances.max_abs_chemical_log_residual
            && residual.abs() <= tolerances.max_abs_chemical_log_residual
    });
    let conservation_agreement = match (
        independent.max_abs_element_balance,
        canonical.max_abs_element_balance,
    ) {
        (Some(independent_balance), Some(canonical_balance)) => Some(
            independent_balance <= tolerances.max_abs_element_balance
                && canonical_balance <= tolerances.max_abs_element_balance,
        ),
        _ => None,
    };
    Ok(PurePhasePhCrossValidationReport {
        identity_agreement: Some(true),
        temperature_agreement,
        thermodynamic_agreement,
        composition_agreement,
        enthalpy_agreement,
        conservation_agreement,
        independent_max_abs_element_balance: independent.max_abs_element_balance,
        canonical_max_abs_element_balance: canonical.max_abs_element_balance,
    })
}

/// Accepts a comparison when its absolute error is small on an absolute scale
/// or on the magnitude scale of the two physical values.
fn within_abs_relative(left: f64, right: f64, absolute: f64, relative: f64) -> bool {
    (left - right).abs() <= absolute + relative * left.abs().max(right.abs())
}

/// Scalar controls for the independent nested P,H route.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PurePhasePhSolverSettings {
    /// Bisection budget for the inner extent solve at a fixed temperature.
    pub max_inner_iterations: usize,
    /// Bisection budget for the outer temperature solve.
    pub max_outer_iterations: usize,
    /// Absolute tolerance on the inner chemical log residual.
    ///
    /// The default is intentionally tighter than the outer enthalpy tolerance:
    /// the outer residual is evaluated on this chemical root, so a loose inner
    /// acceptance can otherwise impose a numerical noise floor on P,H bisection.
    pub max_abs_log_residual: f64,
    /// Absolute tolerance on the outer enthalpy residual, in J.
    pub max_abs_enthalpy_residual: f64,
    /// Relative inset margin used to keep extent brackets strictly interior.
    pub feasibility_margin: f64,
    /// Number of equal temperature subintervals used only when an interior
    /// condensed-phase extent does not exist at one or both global bounds.
    ///
    /// This is a branch-discovery control, not an outer-solver iteration
    /// budget. It lets an independent P,H validator find a bounded interior
    /// branch that terminates at a physical phase boundary instead of treating
    /// the absence of a positive extent at the global endpoint as a failure.
    pub interior_branch_scan_subdivisions: usize,
}

impl Default for PurePhasePhSolverSettings {
    fn default() -> Self {
        Self {
            max_inner_iterations: 128,
            max_outer_iterations: 128,
            max_abs_log_residual: 1e-12,
            max_abs_enthalpy_residual: 1e-7,
            feasibility_margin: 1e-12,
            interior_branch_scan_subdivisions: 64,
        }
    }
}

impl PurePhasePhSolverSettings {
    /// Validates the inner/outer iteration budgets and residual tolerances.
    ///
    /// Both iteration counts must be strictly positive and every tolerance must
    /// be finite and positive. A violation returns an invalid-problem error so
    /// the nested solver never runs with a degenerate budget.
    fn validate(self) -> Result<(), ReactionExtentError> {
        if self.max_inner_iterations == 0 || self.max_outer_iterations == 0 {
            return Err(invalid_problem(
                "inner and outer iteration budgets must be positive",
            ));
        }
        if self.interior_branch_scan_subdivisions < 2 {
            return Err(invalid_problem(
                "interior branch scan must contain at least two subintervals",
            ));
        }
        for (name, value) in [
            ("max_abs_log_residual", self.max_abs_log_residual),
            ("max_abs_enthalpy_residual", self.max_abs_enthalpy_residual),
        ] {
            if !value.is_finite() || value <= 0.0 {
                return Err(invalid_problem(format!(
                    "{name} must be finite and positive"
                )));
            }
        }
        if !self.feasibility_margin.is_finite()
            || self.feasibility_margin <= 0.0
            || self.feasibility_margin >= 0.25
        {
            return Err(invalid_problem(
                "feasibility margin must lie strictly between zero and 0.25",
            ));
        }
        Ok(())
    }
}

/// Accepted inner chemical root at one temperature.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PurePhasePhInnerExtentResult {
    /// Converged reaction extent.
    pub extent: f64,
    /// Chemical `ln(Q)-ln(K)` residual at the converged extent.
    pub log_residual: f64,
    /// Number of inner bisection iterations used.
    pub iterations: usize,
}

/// Compact I1-I2 evidence from the independent nested P,H route.
#[derive(Debug, Clone, PartialEq)]
pub struct PurePhasePhEquilibriumResult {
    /// Physical problem that produced this accepted scalar result.
    pub case_identity: PurePhasePhCaseIdentity,
    /// Solved equilibrium temperature, in K.
    pub temperature: f64,
    /// Converged reaction extent.
    pub extent: f64,
    /// Physical gas mole numbers at the solved state.
    pub gas_moles: Vec<f64>,
    /// Physical mole number of the pure candidate at the solved state.
    pub candidate_moles: f64,
    /// Additive total enthalpy of the solved state, in J.
    pub total_enthalpy: f64,
    /// Chemical `ln(Q)-ln(K)` residual at the solved state.
    pub chemical_log_residual: f64,
    /// `H(T) - H_target` residual of the solved state, in J.
    pub enthalpy_residual: f64,
    /// Independent accepted-state elemental balance, when explicit element
    /// composition was supplied to the validation problem.
    pub max_abs_element_balance: Option<f64>,
    /// Number of outer temperature bisection iterations.
    pub outer_iterations: usize,
    /// Total number of inner extent solves across all temperature trials.
    pub total_inner_solves: usize,
    /// Total number of inner bisection iterations across all trials.
    pub total_inner_iterations: usize,
}

/// Independent nested scalar P,H validator.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PurePhasePhValidator {
    pub settings: PurePhasePhSolverSettings,
}

impl Default for PurePhasePhValidator {
    fn default() -> Self {
        Self {
            settings: PurePhasePhSolverSettings::default(),
        }
    }
}

impl PurePhasePhValidator {
    /// Solves the nested independent extent/temperature problem transactionally.
    pub fn solve(
        &self,
        problem: &PurePhasePhProblem,
    ) -> Result<PurePhasePhEquilibriumResult, ReactionExtentError> {
        self.settings.validate()?;
        let lower =
            self.optional_interior_evaluation(problem, problem.conditions.lower_temperature)?;
        let upper =
            self.optional_interior_evaluation(problem, problem.conditions.upper_temperature)?;
        if let (Some(lower), Some(upper)) = (lower, upper) {
            let total_inner_iterations = lower.inner.iterations + upper.inner.iterations;
            let total_inner_solves = 2;
            if lower.enthalpy_residual.abs() <= self.settings.max_abs_enthalpy_residual {
                let inner_iterations = lower.inner.iterations;
                return Ok(lower.into_result(problem.case_identity(), 0, 1, inner_iterations));
            }
            if upper.enthalpy_residual.abs() <= self.settings.max_abs_enthalpy_residual {
                return Ok(upper.into_result(
                    problem.case_identity(),
                    0,
                    total_inner_solves,
                    total_inner_iterations,
                ));
            }
            if lower.enthalpy_residual.signum() != upper.enthalpy_residual.signum() {
                return self.solve_bracketed_temperature(
                    problem,
                    lower,
                    upper,
                    total_inner_solves,
                    total_inner_iterations,
                );
            }
        }

        self.solve_across_interior_temperature_branches(problem)
    }

    /// Evaluates a chemical interior state, treating only the deliberate
    /// absence of a positive extent root as a branch boundary.
    ///
    /// Thermochemistry, composition, and arithmetic failures remain hard
    /// errors. This distinction prevents branch discovery from masking broken
    /// local data while still modelling the ordinary physical disappearance of
    /// a pure condensed phase at a temperature boundary.
    fn optional_interior_evaluation(
        &self,
        problem: &PurePhasePhProblem,
        temperature: f64,
    ) -> Result<Option<TemperatureEvaluation>, ReactionExtentError> {
        match self.evaluate_temperature(problem, temperature) {
            Ok(evaluation) => Ok(Some(evaluation)),
            Err(ReactionExtentError::ValidationNotApplicable {
                path: "pure_phase_ph_inner_extent",
                ..
            }) => Ok(None),
            Err(error) => Err(error),
        }
    }

    /// Discovers a continuous temperature segment on which a positive
    /// condensed-phase extent exists and brackets the enthalpy root there.
    ///
    /// The scan is used only after the fast global-endpoint bracket was not
    /// applicable. Adjacent valid samples define a continuous interior branch;
    /// a missing inner root breaks that branch rather than allowing a bracket
    /// to jump across a phase boundary. If several branches could satisfy a
    /// target, the first one in the caller's increasing temperature interval
    /// is selected deterministically and reported through its resulting state.
    fn solve_across_interior_temperature_branches(
        &self,
        problem: &PurePhasePhProblem,
    ) -> Result<PurePhasePhEquilibriumResult, ReactionExtentError> {
        let lower = problem.conditions.lower_temperature;
        let upper = problem.conditions.upper_temperature;
        let subdivisions = self.settings.interior_branch_scan_subdivisions;
        let step = (upper - lower) / subdivisions as f64;
        let mut previous: Option<TemperatureEvaluation> = None;
        let mut valid_samples = 0usize;
        let mut total_inner_solves = 0usize;
        let mut total_inner_iterations = 0usize;
        for index in 0..=subdivisions {
            let temperature = lower + step * index as f64;
            let Some(current) = self.optional_interior_evaluation(problem, temperature)? else {
                previous = None;
                continue;
            };
            valid_samples += 1;
            total_inner_solves += 1;
            total_inner_iterations += current.inner.iterations;
            if current.enthalpy_residual.abs() <= self.settings.max_abs_enthalpy_residual {
                return Ok(current.into_result(
                    problem.case_identity(),
                    0,
                    total_inner_solves,
                    total_inner_iterations,
                ));
            }
            if let Some(left) = previous.as_ref()
                && left.enthalpy_residual.signum() != current.enthalpy_residual.signum()
            {
                return self.solve_bracketed_temperature(
                    problem,
                    left.clone(),
                    current,
                    total_inner_solves,
                    total_inner_iterations,
                );
            }
            previous = Some(current);
        }
        let message = if valid_samples == 0 {
            format!("no positive condensed-phase extent root exists inside {lower}..{upper} K")
        } else {
            format!(
                "target enthalpy is not bracketed on any continuous interior branch inside {lower}..{upper} K ({valid_samples} valid scan samples)"
            )
        };
        Err(ReactionExtentError::ValidationNotApplicable {
            path: "pure_phase_ph_outer_temperature",
            message,
        })
    }

    /// Bisects an already validated interior temperature branch.
    fn solve_bracketed_temperature(
        &self,
        problem: &PurePhasePhProblem,
        mut left: TemperatureEvaluation,
        mut right: TemperatureEvaluation,
        mut total_inner_solves: usize,
        mut total_inner_iterations: usize,
    ) -> Result<PurePhasePhEquilibriumResult, ReactionExtentError> {
        for iteration in 1..=self.settings.max_outer_iterations {
            let temperature = 0.5 * (left.temperature + right.temperature);
            let mid = self
                .optional_interior_evaluation(problem, temperature)?
                .ok_or_else(|| ReactionExtentError::ValidationNotApplicable {
                    path: "pure_phase_ph_outer_temperature",
                    message: format!(
                        "interior extent branch disappeared inside an already bracketed temperature interval {}..{} K",
                        left.temperature, right.temperature
                    ),
                })?;
            total_inner_solves += 1;
            total_inner_iterations += mid.inner.iterations;
            if mid.enthalpy_residual.abs() <= self.settings.max_abs_enthalpy_residual {
                return Ok(mid.into_result(
                    problem.case_identity(),
                    iteration,
                    total_inner_solves,
                    total_inner_iterations,
                ));
            }
            if left.enthalpy_residual.signum() != mid.enthalpy_residual.signum() {
                right = mid;
            } else {
                left = mid;
            }
        }
        Err(ReactionExtentError::SolveError(SolveError::MaxIterations))
    }

    /// Evaluates one outer temperature trial for the inner fixed-extent solve.
    ///
    /// Solves the chemical equation at the supplied temperature, materializes
    /// the resulting gas and candidate moles, computes the additive total
    /// enthalpy, and packages everything (plus the enthalpy residual against
    /// the target) into a [`TemperatureEvaluation`] used by the scalar
    /// bracketing loop.
    fn evaluate_temperature(
        &self,
        problem: &PurePhasePhProblem,
        temperature: f64,
    ) -> Result<TemperatureEvaluation, ReactionExtentError> {
        let inner = problem.solve_inner_extent(temperature, self.settings)?;
        let state = problem.materialize_state_at_extent(inner.extent, temperature)?;
        Ok(TemperatureEvaluation {
            temperature,
            inner,
            gas_moles: state.gas_moles,
            candidate_moles: state.candidate_moles,
            total_enthalpy: state.total_enthalpy,
            enthalpy_residual: state.total_enthalpy - problem.conditions.target_enthalpy,
            max_abs_element_balance: state.max_abs_element_balance,
        })
    }
}

#[derive(Debug, Clone)]
struct TemperatureEvaluation {
    temperature: f64,
    inner: PurePhasePhInnerExtentResult,
    gas_moles: Vec<f64>,
    candidate_moles: f64,
    total_enthalpy: f64,
    enthalpy_residual: f64,
    max_abs_element_balance: Option<f64>,
}

impl TemperatureEvaluation {
    /// Converts a completed temperature evaluation into the public equilibrium
    /// result, attaching the outer-loop and accumulated inner-solve counters.
    ///
    /// The conversion is a pure projection: it merges the solved temperature,
    /// extent, moles, enthalpies, and residuals with the iteration bookkeeping
    /// collected by the calling bracketing loop.
    fn into_result(
        self,
        case_identity: PurePhasePhCaseIdentity,
        outer_iterations: usize,
        total_inner_solves: usize,
        total_inner_iterations: usize,
    ) -> PurePhasePhEquilibriumResult {
        PurePhasePhEquilibriumResult {
            case_identity,
            temperature: self.temperature,
            extent: self.inner.extent,
            gas_moles: self.gas_moles,
            candidate_moles: self.candidate_moles,
            total_enthalpy: self.total_enthalpy,
            chemical_log_residual: self.inner.log_residual,
            enthalpy_residual: self.enthalpy_residual,
            max_abs_element_balance: self.max_abs_element_balance,
            outer_iterations,
            total_inner_solves,
            total_inner_iterations,
        }
    }
}

/// Rejects non-finite or non-positive structural element-balance tolerances.
///
/// Kept as a free helper so both the P,H validator and any caller that
/// re-uses the shared boundary composition type validate the same contract
/// with one code path.
fn validate_structural_tolerances(
    tolerances: PurePhaseBoundaryStructuralTolerances,
) -> Result<(), ReactionExtentError> {
    for (name, value) in [
        (
            "max_abs_element_balance",
            tolerances.max_abs_element_balance,
        ),
        (
            "rank_absolute_tolerance",
            tolerances.rank_absolute_tolerance,
        ),
        (
            "rank_relative_tolerance",
            tolerances.rank_relative_tolerance,
        ),
    ] {
        if !value.is_finite() || value <= 0.0 {
            return Err(invalid_problem(format!(
                "{name} must be finite and positive"
            )));
        }
    }
    Ok(())
}

fn invalid_problem(message: impl Into<String>) -> ReactionExtentError {
    ReactionExtentError::InvalidProblem {
        field: "pure_phase_ph_validator",
        message: message.into(),
    }
}

fn invalid_candidate(message: impl Into<String>) -> ReactionExtentError {
    ReactionExtentError::InvalidCandidate {
        field: "pure_phase_ph_validator",
        message: message.into(),
    }
}
