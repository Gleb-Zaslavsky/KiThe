//! Independent phase-boundary validation through an equilibrium constant.
//!
//! # Purpose
//!
//! This module provides an independent validation path for the simplest
//! phase-appearance problems:
//!
//! - one ideal-gas phase;
//! - one pure one-component condensed candidate phase;
//! - one independent phase-forming reaction;
//! - fixed P,T.
//!
//! The validator deliberately does NOT use:
//!
//! - the canonical log-moles residual;
//! - the analytical Jacobian;
//! - elemental-potential reconstruction;
//! - TPD;
//! - PhaseManager;
//! - ActiveSet / outer-loop logic.
//!
//! Instead, it uses the ordinary reaction relation
//!
//!     Δ_r G = R T [ln(Q) - ln(K)]
//!
//! and exploits the fact that the activity of a pure condensed phase is one.
//!
//! If the reaction is oriented so that the pure candidate phase has a
//! positive stoichiometric coefficient, then at zero candidate amount
//!
//!     ln(Q) - ln(K) < 0
//!
//! means that creating the phase lowers Gibbs energy.
//!
//! This makes the module useful as an independent cross-validator for the
//! canonical TPD + ActiveSet phase-selection pipeline.
//!
//! # Important scope restriction
//!
//! This is intentionally NOT a general multiphase equilibrium solver.
//! The full system must have exactly one independent reaction direction.
//! The module exists to validate simple phase boundaries and simple
//! two-phase equilibria by "a second mathematics".

use crate::Thermodynamics::ChemEquilibrium::equilibrium_log_moles::GibbsFn;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_nonlinear::ReactionExtentError;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_phase_stability::{
    PhaseStabilityReport, PhaseStabilityStatus,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_problem::EquilibriumConditions;

use nalgebra::DMatrix;
use std::fmt;

/// Universal gas constant in J/(mol*K).
///
/// If KiThe already exposes one canonical constant, import that one instead.
pub const MOLAR_GAS_CONSTANT: f64 = 8.314_462_618_153_24;

/// Independent elemental composition for a pure-phase boundary problem.
///
/// Rows of `gas_species_by_element` follow `PurePhaseBoundaryProblem` gas
/// species order; columns follow `elements`. The candidate composition is one
/// additional row in the same element order. This data belongs to the
/// independent validator rather than the canonical reaction-basis builder.
#[derive(Debug, Clone, PartialEq)]
pub struct PurePhaseBoundaryElementComposition {
    elements: Vec<String>,
    gas_species_by_element: DMatrix<f64>,
    candidate_by_element: Vec<f64>,
}

impl PurePhaseBoundaryElementComposition {
    /// Creates finite non-negative elemental composition data. The gas-row
    /// count is checked when it is attached to a concrete boundary problem.
    pub fn new(
        elements: Vec<String>,
        gas_species_by_element: DMatrix<f64>,
        candidate_by_element: Vec<f64>,
    ) -> Result<Self, ReactionExtentError> {
        if elements.is_empty() {
            return Err(invalid_problem(
                "element composition requires at least one element",
            ));
        }
        if elements.iter().any(|element| element.trim().is_empty()) {
            return Err(invalid_problem("element names must be non-empty"));
        }
        let mut unique_elements = std::collections::BTreeSet::new();
        if elements
            .iter()
            .any(|element| !unique_elements.insert(element.clone()))
        {
            return Err(invalid_problem("element names must be unique"));
        }
        if gas_species_by_element.ncols() != elements.len()
            || candidate_by_element.len() != elements.len()
        {
            return Err(invalid_problem(format!(
                "element composition has {} columns and {} candidate entries for {} elements",
                gas_species_by_element.ncols(),
                candidate_by_element.len(),
                elements.len()
            )));
        }
        if gas_species_by_element.nrows() == 0 {
            return Err(invalid_problem(
                "element composition requires at least one gas-species row",
            ));
        }
        if gas_species_by_element
            .iter()
            .chain(candidate_by_element.iter())
            .any(|value| !value.is_finite() || *value < 0.0)
        {
            return Err(invalid_problem(
                "elemental composition entries must be finite and non-negative",
            ));
        }
        if candidate_by_element.iter().all(|value| *value == 0.0) {
            return Err(invalid_problem(
                "candidate phase must contain at least one element",
            ));
        }

        Ok(Self {
            elements,
            gas_species_by_element,
            candidate_by_element,
        })
    }

    pub fn elements(&self) -> &[String] {
        &self.elements
    }

    pub fn gas_species_by_element(&self) -> &DMatrix<f64> {
        &self.gas_species_by_element
    }

    pub fn candidate_by_element(&self) -> &[f64] {
        &self.candidate_by_element
    }
}

/// Explicit tolerances for independent structural validation.
///
/// These values govern reaction conservation and numerical matrix rank only.
/// They are intentionally unrelated to residual acceptance or production
/// phase hysteresis thresholds.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PurePhaseBoundaryStructuralTolerances {
    /// Maximum accepted absolute elemental-balance residual.
    pub max_abs_element_balance: f64,
    /// Absolute tolerance used to decide numerical matrix rank.
    pub rank_absolute_tolerance: f64,
    /// Relative tolerance used to decide numerical matrix rank.
    pub rank_relative_tolerance: f64,
}

impl Default for PurePhaseBoundaryStructuralTolerances {
    fn default() -> Self {
        Self {
            max_abs_element_balance: 1e-12,
            rank_absolute_tolerance: 1e-12,
            rank_relative_tolerance: 1e-10,
        }
    }
}

impl PurePhaseBoundaryStructuralTolerances {
    /// Validates that every structural tolerance is finite and positive.
    ///
    /// Returns the validated tolerances so callers can chain the check into
    /// problem construction. Any non-finite or non-positive value is rejected
    /// as an invalid problem rather than being silently used downstream.
    fn validate(self) -> Result<Self, ReactionExtentError> {
        for (name, value) in [
            ("max_abs_element_balance", self.max_abs_element_balance),
            ("rank_absolute_tolerance", self.rank_absolute_tolerance),
            ("rank_relative_tolerance", self.rank_relative_tolerance),
        ] {
            if !value.is_finite() || value <= 0.0 {
                return Err(invalid_problem(format!(
                    "{name} must be finite and strictly positive"
                )));
            }
        }
        Ok(self)
    }
}

/// Independent rank evidence for the supplied gas-plus-pure-phase family.
#[derive(Debug, Clone, PartialEq)]
pub struct PurePhaseBoundaryReactionSpace {
    /// Rank of the combined gas-plus-candidate element matrix.
    pub full_rank: usize,
    /// Rank of the gas-only element matrix.
    pub gas_only_rank: usize,
    /// Dimension of the full-system reaction space.
    pub full_reaction_dimension: usize,
    /// Dimension of the gas-only reaction space.
    pub gas_only_reaction_dimension: usize,
    /// Per-element conservation residual of the supplied reaction.
    pub element_balance_residuals: Vec<f64>,
}

// ============================================================================
// Problem definition
// ============================================================================

/// Independent equilibrium-constant problem for appearance of one pure phase.
///
/// The reaction is represented as
///
///     Σ_i ν_i A_i(g) + ν_s S(condensed) = 0
///
/// and MUST be oriented so that
///
///     ν_s > 0.
///
/// Therefore positive reaction extent means creation of the candidate
/// condensed phase.
///
/// `gas_moles_at_absence` is the physical gas composition at ξ = 0, i.e. at
/// exactly zero amount of the candidate phase. In a cross-validation workflow
/// this will normally be the accepted gas-only equilibrium composition
/// obtained before the canonical outer loop activates the phase.
///
/// The pure candidate does not occur in the reaction quotient because
///
///     a_s = 1,
///     ln(a_s) = 0.
///
/// It DOES occur in ΔG° and therefore in K.
#[derive(Clone)]
pub struct PurePhaseBoundaryProblem {
    /// Names of gas species, used only for diagnostics/reporting.
    gas_species: Vec<String>,

    /// Physical gas mole numbers at zero candidate-phase amount.
    ///
    /// These values must be strictly positive because `ln(Q)` is evaluated
    /// directly from ideal-gas activities.
    gas_moles_at_absence: Vec<f64>,

    /// Stoichiometric coefficients of gas species.
    ///
    /// Sign convention follows
    ///
    ///     n_i(ξ) = n_i(0) + ν_i ξ.
    gas_stoichiometry: Vec<f64>,

    /// Stoichiometric coefficient of the pure candidate phase.
    ///
    /// It is required to be positive, so `ξ > 0` creates the candidate.
    candidate_stoichiometry: f64,

    /// Standard-state Gibbs functions for gas species, in the same order as
    /// `gas_species`.
    gas_standard_gibbs: Vec<GibbsFn>,

    /// Standard-state Gibbs function of the pure condensed candidate phase.
    candidate_standard_gibbs: GibbsFn,

    /// Thermodynamic conditions for the boundary calculation.
    conditions: EquilibriumConditions,

    /// Optional human-readable name of the candidate phase.
    candidate_name: String,

    /// Optional independent elemental composition. Generic K_eq validation
    /// remains available without it, while strict family validation requires
    /// this additional structural evidence.
    element_composition: Option<PurePhaseBoundaryElementComposition>,
}

impl PurePhaseBoundaryProblem {
    /// Constructs and validates the independent phase-boundary problem.
    ///
    /// This constructor intentionally accepts the reaction directly rather
    /// than deriving it from the canonical reaction basis. For a truly
    /// independent validator that is useful: a test fixture can provide a
    /// hand-validated reaction without depending on the same basis builder
    /// used by the canonical equilibrium path.
    pub fn new(
        gas_species: Vec<String>,
        gas_moles_at_absence: Vec<f64>,
        gas_stoichiometry: Vec<f64>,
        candidate_stoichiometry: f64,
        gas_standard_gibbs: Vec<GibbsFn>,
        candidate_standard_gibbs: GibbsFn,
        conditions: EquilibriumConditions,
        candidate_name: impl Into<String>,
    ) -> Result<Self, ReactionExtentError> {
        let n = gas_species.len();

        if n == 0 {
            return Err(invalid_problem("at least one gas species is required"));
        }

        if gas_moles_at_absence.len() != n
            || gas_stoichiometry.len() != n
            || gas_standard_gibbs.len() != n
        {
            return Err(invalid_problem(format!(
                "gas species count is {n}, but moles/stoichiometry/Gibbs lengths are {}/{}/{}",
                gas_moles_at_absence.len(),
                gas_stoichiometry.len(),
                gas_standard_gibbs.len(),
            )));
        }

        if gas_moles_at_absence
            .iter()
            .any(|n_i| !n_i.is_finite() || *n_i <= 0.0)
        {
            return Err(invalid_problem(
                "gas mole numbers at phase absence must be finite and strictly positive",
            ));
        }

        if gas_stoichiometry.iter().any(|nu| !nu.is_finite()) {
            return Err(invalid_problem(
                "gas stoichiometric coefficients must be finite",
            ));
        }

        if !candidate_stoichiometry.is_finite() || candidate_stoichiometry <= 0.0 {
            return Err(invalid_problem(
                "candidate stoichiometric coefficient must be finite and strictly positive",
            ));
        }

        // A reaction with no gas contribution is not useful for the present
        // validator: Q would be identically one.
        if gas_stoichiometry.iter().all(|nu| *nu == 0.0) {
            return Err(invalid_problem(
                "at least one gas stoichiometric coefficient must be non-zero",
            ));
        }

        Ok(Self {
            gas_species,
            gas_moles_at_absence,
            gas_stoichiometry,
            candidate_stoichiometry,
            gas_standard_gibbs,
            candidate_standard_gibbs,
            conditions,
            candidate_name: candidate_name.into(),
            element_composition: None,
        })
    }

    /// Attaches independent elemental composition and rejects a supplied
    /// phase-forming reaction that does not conserve every declared element.
    ///
    /// This is deliberately a fallible builder instead of a hidden repair
    /// step: invalid stoichiometry is a fixture/input error, never something
    /// the validator may project into a different reaction.
    pub fn with_element_composition(
        mut self,
        composition: PurePhaseBoundaryElementComposition,
        tolerances: PurePhaseBoundaryStructuralTolerances,
    ) -> Result<Self, ReactionExtentError> {
        let tolerances = tolerances.validate()?;
        if composition.gas_species_by_element.nrows() != self.gas_species.len() {
            return Err(invalid_problem(format!(
                "element composition has {} gas rows for {} gas species",
                composition.gas_species_by_element.nrows(),
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

    pub fn gas_moles_at_absence(&self) -> &[f64] {
        &self.gas_moles_at_absence
    }

    pub fn gas_stoichiometry(&self) -> &[f64] {
        &self.gas_stoichiometry
    }

    pub fn candidate_stoichiometry(&self) -> f64 {
        self.candidate_stoichiometry
    }

    pub fn conditions(&self) -> EquilibriumConditions {
        self.conditions
    }

    /// Returns a stable identity for the structural and state data supplied to
    /// the independent validation problem.
    ///
    /// Gibbs closures cannot be compared by identity in Rust. Callers that
    /// consume an independent result therefore use this layout/state identity
    /// and re-evaluate its boundary report before accepting it as evidence.
    pub fn case_identity(&self) -> PurePhaseBoundaryCaseIdentity {
        PurePhaseBoundaryCaseIdentity {
            gas_species: self.gas_species.clone(),
            candidate_name: self.candidate_name.clone(),
            gas_mole_bits: self
                .gas_moles_at_absence
                .iter()
                .map(|value| value.to_bits())
                .collect(),
            gas_stoichiometry_bits: self
                .gas_stoichiometry
                .iter()
                .map(|value| value.to_bits())
                .collect(),
            candidate_stoichiometry_bits: self.candidate_stoichiometry.to_bits(),
            temperature_bits: self.conditions.temperature().to_bits(),
            pressure_bits: self.conditions.pressure().to_bits(),
            reference_pressure_bits: self.conditions.reference_pressure().to_bits(),
        }
    }

    /// Optional independent elemental composition attached to this problem.
    pub fn element_composition(&self) -> Option<&PurePhaseBoundaryElementComposition> {
        self.element_composition.as_ref()
    }

    /// Computes rank and reaction-space dimensions from independently supplied
    /// element matrices. This does not use the canonical reaction basis.
    pub fn reaction_space(
        &self,
        tolerances: PurePhaseBoundaryStructuralTolerances,
    ) -> Result<PurePhaseBoundaryReactionSpace, ReactionExtentError> {
        let tolerances = tolerances.validate()?;
        let composition = self.element_composition.as_ref().ok_or_else(|| {
            ReactionExtentError::ValidationNotApplicable {
                path: "pure_phase_boundary_structure",
                message: "strict reaction-space validation requires element composition".into(),
            }
        })?;
        let gas_only_rank =
            independent_matrix_rank(composition.gas_species_by_element(), tolerances)?;
        let mut full = DMatrix::zeros(
            composition.gas_species_by_element.nrows() + 1,
            composition.gas_species_by_element.ncols(),
        );
        full.rows_mut(0, composition.gas_species_by_element.nrows())
            .copy_from(composition.gas_species_by_element());
        for (column, value) in composition.candidate_by_element().iter().enumerate() {
            full[(composition.gas_species_by_element.nrows(), column)] = *value;
        }
        let full_rank = independent_matrix_rank(&full, tolerances)?;
        let element_balance_residuals = self.element_balance_residuals_for(composition)?;

        Ok(PurePhaseBoundaryReactionSpace {
            full_rank,
            gas_only_rank,
            full_reaction_dimension: full.nrows().saturating_sub(full_rank),
            gas_only_reaction_dimension: composition
                .gas_species_by_element
                .nrows()
                .saturating_sub(gas_only_rank),
            element_balance_residuals,
        })
    }

    /// Validates the narrow independent family used by the synthetic
    /// cross-validation suite: exactly one full-system reaction direction and
    /// no residual gas-only reaction direction when the candidate is absent.
    ///
    /// For one appended pure candidate, a conserved non-zero phase-forming
    /// direction makes `full_dimension = gas_only_dimension + 1`. The gas-only
    /// condition is therefore retained as explicit defensive evidence and a
    /// clearer diagnostic, not presented as an unrelated mathematical premise.
    pub fn validate_strict_independent_family(
        &self,
        tolerances: PurePhaseBoundaryStructuralTolerances,
    ) -> Result<PurePhaseBoundaryReactionSpace, ReactionExtentError> {
        let tolerances = tolerances.validate()?;
        let reaction_space = self.reaction_space(tolerances)?;
        if reaction_space
            .element_balance_residuals
            .iter()
            .any(|residual| residual.abs() > tolerances.max_abs_element_balance)
        {
            return Err(invalid_problem(format!(
                "strict family reaction violates elemental conservation: residuals={:?}",
                reaction_space.element_balance_residuals
            )));
        }
        if reaction_space.full_reaction_dimension != 1 {
            return Err(invalid_problem(format!(
                "strict family requires one full-system reaction dimension, got {}",
                reaction_space.full_reaction_dimension
            )));
        }
        if reaction_space.gas_only_reaction_dimension != 0 {
            return Err(invalid_problem(format!(
                "strict family requires no gas-only reaction dimension, got {}",
                reaction_space.gas_only_reaction_dimension
            )));
        }
        Ok(reaction_space)
    }

    /// Evaluates the per-element conservation residual of the supplied
    /// elemental composition under this problem's fixed reaction.
    ///
    /// Requires the composition's gas-species rows to match this problem's
    /// gas-stoichiometry ordering. The residual for each element column is the
    /// sum of `nu_i * A[i,e]` over gas species plus the candidate contribution
    /// `nu_s * A_s[e]`. A non-finite result is rejected as an invalid problem.
    fn element_balance_residuals_for(
        &self,
        composition: &PurePhaseBoundaryElementComposition,
    ) -> Result<Vec<f64>, ReactionExtentError> {
        if composition.gas_species_by_element.nrows() != self.gas_stoichiometry.len() {
            return Err(invalid_problem(
                "cannot evaluate element balance with incompatible gas composition rows",
            ));
        }
        let mut residuals = vec![0.0; composition.gas_species_by_element.ncols()];
        for (gas_index, coefficient) in self.gas_stoichiometry.iter().enumerate() {
            for element_index in 0..residuals.len() {
                residuals[element_index] +=
                    coefficient * composition.gas_species_by_element[(gas_index, element_index)];
            }
        }
        for (element_index, candidate) in composition.candidate_by_element.iter().enumerate() {
            residuals[element_index] += self.candidate_stoichiometry * candidate;
        }
        if residuals.iter().any(|value| !value.is_finite()) {
            return Err(invalid_problem(
                "reaction elemental-balance residual is non-finite",
            ));
        }
        Ok(residuals)
    }

    /// Standard Gibbs energy of the phase-forming reaction:
    ///
    ///     ΔG°_r =
    ///         Σ_gas ν_i g_i°(T)
    ///         + ν_s g_s°(T).
    pub fn standard_reaction_gibbs(&self) -> Result<f64, ReactionExtentError> {
        let t = self.conditions.temperature();

        let mut delta_g = 0.0;

        for (index, (nu, gibbs)) in self
            .gas_stoichiometry
            .iter()
            .zip(&self.gas_standard_gibbs)
            .enumerate()
        {
            let value = gibbs(t);

            if !value.is_finite() {
                return Err(ReactionExtentError::InvalidDG0 {
                    species_index: index,
                    dg0: value,
                    temperature: t,
                });
            }

            delta_g += nu * value;
        }

        let candidate_gibbs = (self.candidate_standard_gibbs)(t);

        if !candidate_gibbs.is_finite() {
            // The independent boundary layout consists of all gas components
            // followed by the one pure condensed candidate. Reusing the
            // typed thermochemical error keeps candidate and gas closures
            // equally diagnosable without a string-only special case.
            return Err(ReactionExtentError::InvalidDG0 {
                species_index: self.gas_species.len(),
                dg0: candidate_gibbs,
                temperature: t,
            });
        }

        delta_g += self.candidate_stoichiometry * candidate_gibbs;

        Ok(delta_g)
    }

    /// Natural logarithm of the dimensionless equilibrium constant.
    ///
    ///     ln K = -ΔG°_r / RT.
    pub fn ln_equilibrium_constant(&self) -> Result<f64, ReactionExtentError> {
        Ok(-self.standard_reaction_gibbs()? / (MOLAR_GAS_CONSTANT * self.conditions.temperature()))
    }

    /// Reconstructs gas mole numbers at reaction extent ξ:
    ///
    ///     n_i(ξ) = n_i(0) + ν_i ξ.
    ///
    /// The pure condensed candidate has
    ///
    ///     n_s(ξ) = ν_s ξ,
    ///
    /// but it is not included in the returned gas vector.
    pub fn gas_moles_at_extent(&self, extent: f64) -> Result<Vec<f64>, ReactionExtentError> {
        if !extent.is_finite() || extent < 0.0 {
            return Err(invalid_candidate(
                "phase-forming extent must be finite and non-negative",
            ));
        }

        let mut result = Vec::with_capacity(self.gas_moles_at_absence.len());

        for (index, (&n0, &nu)) in self
            .gas_moles_at_absence
            .iter()
            .zip(&self.gas_stoichiometry)
            .enumerate()
        {
            let value = n0 + nu * extent;

            // This narrow validator deliberately solves only interior
            // gas-phase equilibria. Exact disappearance of a gas component
            // is treated as a boundary case outside its present contract.
            if !value.is_finite() || value <= 0.0 {
                return Err(ReactionExtentError::InvalidCandidate {
                    field: "phase_boundary_gas_moles",
                    message: format!("gas species {index} is non-positive at extent {extent}"),
                });
            }

            result.push(value);
        }

        Ok(result)
    }

    /// Mole number of the pure candidate phase at extent ξ.
    pub fn candidate_moles_at_extent(&self, extent: f64) -> Result<f64, ReactionExtentError> {
        if !extent.is_finite() || extent < 0.0 {
            return Err(invalid_candidate(
                "phase-forming extent must be finite and non-negative",
            ));
        }

        Ok(self.candidate_stoichiometry * extent)
    }

    /// Natural logarithm of the reaction quotient for the gas composition.
    ///
    /// For an ideal gas,
    ///
    ///     a_i = x_i P/P°.
    ///
    /// The pure candidate contributes nothing because `ln(a_s) = ln(1) = 0`.
    pub fn ln_reaction_quotient_for_gas_moles(
        &self,
        gas_moles: &[f64],
    ) -> Result<f64, ReactionExtentError> {
        if gas_moles.len() != self.gas_species.len() {
            return Err(ReactionExtentError::DimensionMismatch(format!(
                "phase-boundary quotient has {} gas mole values for {} gas species",
                gas_moles.len(),
                self.gas_species.len()
            )));
        }

        if gas_moles.iter().any(|n_i| !n_i.is_finite() || *n_i <= 0.0) {
            return Err(invalid_candidate(
                "reaction quotient requires finite, strictly positive gas moles",
            ));
        }

        let gas_total: f64 = gas_moles.iter().sum();

        if !gas_total.is_finite() || gas_total <= 0.0 {
            return Err(invalid_candidate(
                "total gas mole number must be finite and strictly positive",
            ));
        }

        let pressure_ratio = self.conditions.pressure() / self.conditions.reference_pressure();

        let mut ln_q = 0.0;

        for (&nu, &n_i) in self.gas_stoichiometry.iter().zip(gas_moles) {
            if nu == 0.0 {
                continue;
            }

            let x_i = n_i / gas_total;
            let activity = x_i * pressure_ratio;

            if !activity.is_finite() || activity <= 0.0 {
                return Err(invalid_candidate(
                    "ideal-gas activity must be finite and strictly positive",
                ));
            }

            ln_q += nu * activity.ln();
        }

        Ok(ln_q)
    }

    /// `ln(Q) - ln(K)` at reaction extent ξ.
    pub fn log_residual_at_extent(&self, extent: f64) -> Result<f64, ReactionExtentError> {
        let gas_moles = self.gas_moles_at_extent(extent)?;

        Ok(
            self.ln_reaction_quotient_for_gas_moles(&gas_moles)?
                - self.ln_equilibrium_constant()?,
        )
    }

    /// Upper positive extent allowed by gas-mole positivity.
    ///
    /// Since the candidate is created for ξ >= 0, only gas species with
    /// negative stoichiometric coefficients restrict the upper bound:
    ///
    ///     n_i(0) + ν_i ξ > 0.
    pub fn positive_extent_upper_bound(&self) -> Result<f64, ReactionExtentError> {
        let mut upper = f64::INFINITY;

        for (&n0, &nu) in self
            .gas_moles_at_absence
            .iter()
            .zip(&self.gas_stoichiometry)
        {
            if nu < 0.0 {
                upper = upper.min(n0 / (-nu));
            }
        }

        if !upper.is_finite() || upper <= 0.0 {
            return Err(ReactionExtentError::ValidationNotApplicable {
                path: "pure_phase_boundary_validator",
                message: "no finite positive phase-forming extent bound was found".to_string(),
            });
        }

        Ok(upper)
    }
}

// ============================================================================
// Boundary evidence
// ============================================================================

/// Independent thermodynamic prediction at zero candidate-phase amount.
///
/// The reaction orientation is part of the contract: positive extent creates
/// the pure candidate. Therefore:
///
/// - residual < 0  -> candidate formation lowers G;
/// - residual = 0  -> phase boundary;
/// - residual > 0  -> candidate remains absent.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PurePhaseBoundaryPrediction {
    /// The phase-forming direction is uphill at zero extent; candidate absent.
    StableInactive,
    /// The boundary residual is zero at the tolerance; an exact phase boundary.
    Boundary,
    /// Formation lowers Gibbs energy; the candidate should appear.
    ShouldAppear,
}

/// Discrete phase-state expectation implied by independent boundary evidence.
///
/// At an exact thermodynamic boundary neither active nor inactive is the only
/// correct topology: production hysteresis and accepted phase-set history own
/// that decision. The independent validator must expose this ambiguity instead
/// of smuggling `Boundary` into the inactive branch.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PurePhaseTopologyExpectation {
    /// Independent evidence requires the candidate to remain absent.
    MustBeInactive,
    /// An exact boundary; production hysteresis/history owns the decision.
    HysteresisDependent,
    /// Independent evidence requires the candidate to be present.
    MustBeActive,
}

impl PurePhaseBoundaryPrediction {
    /// Converts the continuous thermodynamic prediction into the strongest
    /// topology statement that can be made without a phase-control history.
    pub fn topology_expectation(self) -> PurePhaseTopologyExpectation {
        match self {
            Self::StableInactive => PurePhaseTopologyExpectation::MustBeInactive,
            Self::Boundary => PurePhaseTopologyExpectation::HysteresisDependent,
            Self::ShouldAppear => PurePhaseTopologyExpectation::MustBeActive,
        }
    }
}

/// Tolerance for classifying the zero-extent boundary residual.
///
/// This tolerance belongs to the independent validator. It should NOT reuse
/// `dg_create` or `dg_keep`, because those belong to the production phase
/// hysteresis policy and would destroy the independence of the validation path.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PurePhaseBoundaryTolerances {
    /// Absolute tolerance classifying the zero-extent boundary residual.
    pub max_abs_boundary_log_residual: f64,
    /// Absolute tolerance on the finite two-phase extent residual.
    pub max_abs_equilibrium_log_residual: f64,
}

impl Default for PurePhaseBoundaryTolerances {
    fn default() -> Self {
        Self {
            max_abs_boundary_log_residual: 1e-8,
            max_abs_equilibrium_log_residual: 1e-10,
        }
    }
}

impl PurePhaseBoundaryTolerances {
    /// Validates the boundary and two-phase-equilibrium acceptance tolerances.
    ///
    /// Rejects any non-finite or non-positive tolerance and returns the
    /// validated value object. This guarantees the boundary validator never
    /// compares residuals against a disabled or degenerate threshold.
    fn validate(self) -> Result<Self, ReactionExtentError> {
        for (name, value) in [
            (
                "max_abs_boundary_log_residual",
                self.max_abs_boundary_log_residual,
            ),
            (
                "max_abs_equilibrium_log_residual",
                self.max_abs_equilibrium_log_residual,
            ),
        ] {
            if !value.is_finite() || value <= 0.0 {
                return Err(invalid_problem(format!(
                    "{name} must be finite and strictly positive"
                )));
            }
        }

        Ok(self)
    }
}

/// Independent evidence for the phase boundary at ξ = 0.
#[derive(Debug, Clone, PartialEq)]
pub struct PurePhaseBoundaryReport {
    /// Name of the pure condensed candidate phase.
    pub candidate_name: String,
    /// Thermodynamic conditions of the boundary evaluation.
    pub conditions: EquilibriumConditions,

    /// ln(Q) evaluated for the gas-only state at zero candidate amount.
    pub ln_q_at_absence: f64,

    /// Thermodynamic ln(K) of the phase-forming reaction.
    pub ln_k: f64,

    /// `ln(Q) - ln(K)` at candidate absence.
    pub log_residual_at_absence: f64,

    /// Same driving force in J/mol:
    ///
    ///     Δ_r G(ξ=0) = RT [ln(Q)-ln(K)].
    pub reaction_gibbs_at_absence: f64,

    /// Classification of the zero-extent boundary state.
    pub prediction: PurePhaseBoundaryPrediction,
}

/// Evaluates the independent phase-boundary criterion without solving for a
/// finite amount of the new phase.
pub fn evaluate_pure_phase_boundary(
    problem: &PurePhaseBoundaryProblem,
    tolerances: PurePhaseBoundaryTolerances,
) -> Result<PurePhaseBoundaryReport, ReactionExtentError> {
    let tolerances = tolerances.validate()?;

    let ln_q = problem.ln_reaction_quotient_for_gas_moles(problem.gas_moles_at_absence())?;
    let ln_k = problem.ln_equilibrium_constant()?;
    let residual = ln_q - ln_k;

    let prediction = if residual.abs() <= tolerances.max_abs_boundary_log_residual {
        PurePhaseBoundaryPrediction::Boundary
    } else if residual < 0.0 {
        PurePhaseBoundaryPrediction::ShouldAppear
    } else {
        PurePhaseBoundaryPrediction::StableInactive
    };

    let reaction_gibbs = MOLAR_GAS_CONSTANT * problem.conditions().temperature() * residual;

    Ok(PurePhaseBoundaryReport {
        candidate_name: problem.candidate_name().to_string(),
        conditions: problem.conditions(),
        ln_q_at_absence: ln_q,
        ln_k,
        log_residual_at_absence: residual,
        reaction_gibbs_at_absence: reaction_gibbs,
        prediction,
    })
}

/// Settings for the independent temperature bisection used by boundary tests.
///
/// This finder operates on `ln(Q)-ln(K)` only. It is intentionally separate
/// from the scalar finite-extent solve below: one locates a phase boundary in
/// temperature, while the other locates a finite equilibrium extent at fixed
/// temperature.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PurePhaseBoundaryTemperatureSearchSettings {
    /// Maximum bisection iterations for the boundary-temperature search.
    pub max_iterations: usize,
}

impl Default for PurePhaseBoundaryTemperatureSearchSettings {
    fn default() -> Self {
        Self { max_iterations: 96 }
    }
}

impl PurePhaseBoundaryTemperatureSearchSettings {
    /// Validates the temperature bisection search settings.
    ///
    /// Requires at least one iteration so the scalar boundary-temperature root
    /// search always terminates. Returns the validated settings for chaining.
    fn validate(self) -> Result<Self, ReactionExtentError> {
        if self.max_iterations == 0 {
            return Err(invalid_problem(
                "boundary temperature search requires at least one iteration",
            ));
        }
        Ok(self)
    }
}

/// Independent root evidence for one pure-phase boundary temperature.
#[derive(Debug, Clone, PartialEq)]
pub struct PurePhaseBoundaryTemperatureRoot {
    /// Root boundary temperature in K.
    pub temperature: f64,
    /// `ln(Q)-ln(K)` residual at the root.
    pub log_residual: f64,
    /// Bisection iterations used to locate the root.
    pub iterations: usize,
}

/// Locates a bracketed pure-phase boundary using only `ln(Q)-ln(K)`.
///
/// `problem_at_temperature` must construct an independent boundary problem at
/// the requested temperature. The function deliberately accepts a factory
/// instead of mutating a reusable production problem, so a synthetic fixture
/// can supply its own analytical Gibbs functions and its own data lifecycle.
/// The endpoints must have opposite residual signs (or one endpoint must
/// already satisfy the supplied boundary tolerance).
pub fn bisect_pure_phase_boundary_temperature<F>(
    lower_temperature: f64,
    upper_temperature: f64,
    settings: PurePhaseBoundaryTemperatureSearchSettings,
    tolerances: PurePhaseBoundaryTolerances,
    problem_at_temperature: F,
) -> Result<PurePhaseBoundaryTemperatureRoot, ReactionExtentError>
where
    F: Fn(f64) -> Result<PurePhaseBoundaryProblem, ReactionExtentError>,
{
    let settings = settings.validate()?;
    let tolerances = tolerances.validate()?;
    if !lower_temperature.is_finite()
        || !upper_temperature.is_finite()
        || lower_temperature <= 0.0
        || upper_temperature <= lower_temperature
    {
        return Err(invalid_problem(
            "boundary temperature bracket must be finite, positive, and ordered",
        ));
    }

    let evaluate_at = |temperature: f64| -> Result<PurePhaseBoundaryReport, ReactionExtentError> {
        let problem = problem_at_temperature(temperature)?;
        if problem.conditions().temperature() != temperature {
            return Err(invalid_problem(format!(
                "boundary problem factory returned {} K for requested {temperature} K",
                problem.conditions().temperature()
            )));
        }
        evaluate_pure_phase_boundary(&problem, tolerances)
    };

    let lower_report = evaluate_at(lower_temperature)?;
    if lower_report.log_residual_at_absence.abs() <= tolerances.max_abs_boundary_log_residual {
        return Ok(PurePhaseBoundaryTemperatureRoot {
            temperature: lower_temperature,
            log_residual: lower_report.log_residual_at_absence,
            iterations: 0,
        });
    }
    let upper_report = evaluate_at(upper_temperature)?;
    if upper_report.log_residual_at_absence.abs() <= tolerances.max_abs_boundary_log_residual {
        return Ok(PurePhaseBoundaryTemperatureRoot {
            temperature: upper_temperature,
            log_residual: upper_report.log_residual_at_absence,
            iterations: 0,
        });
    }
    if lower_report.log_residual_at_absence.signum()
        == upper_report.log_residual_at_absence.signum()
    {
        return Err(ReactionExtentError::ValidationNotApplicable {
            path: "pure_phase_boundary_temperature",
            message: "boundary temperature bracket does not contain a sign change".into(),
        });
    }

    let mut left_temperature = lower_temperature;
    let mut left_residual = lower_report.log_residual_at_absence;
    let mut right_temperature = upper_temperature;

    for iteration in 1..=settings.max_iterations {
        let temperature = 0.5 * (left_temperature + right_temperature);
        let report = evaluate_at(temperature)?;
        let residual = report.log_residual_at_absence;
        if residual.abs() <= tolerances.max_abs_boundary_log_residual {
            return Ok(PurePhaseBoundaryTemperatureRoot {
                temperature,
                log_residual: residual,
                iterations: iteration,
            });
        }

        if left_residual.signum() != residual.signum() {
            right_temperature = temperature;
        } else {
            left_temperature = temperature;
            left_residual = residual;
        }
    }

    Err(ReactionExtentError::SolveError(
        crate::Thermodynamics::ChemEquilibrium::equilibrium_nonlinear::SolveError::MaxIterations,
    ))
}

// ============================================================================
// Optional finite two-phase equilibrium solve
// ============================================================================

/// Numerical settings for the independent scalar extent solve.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PurePhaseBoundarySolverSettings {
    /// Maximum bisection iterations for the finite-extent solve.
    pub max_iterations: usize,

    /// Number of interior samples used to find a sign-changing bracket.
    pub bracket_samples: usize,

    /// Relative shrink from the upper positivity boundary.
    pub feasibility_margin: f64,
}

impl Default for PurePhaseBoundarySolverSettings {
    fn default() -> Self {
        Self {
            max_iterations: 96,
            bracket_samples: 32,
            feasibility_margin: 1e-12,
        }
    }
}

impl PurePhaseBoundarySolverSettings {
    /// Validates the two-phase equilibrium solver settings.
    ///
    /// Requires a strictly positive iteration budget and at least two bracket
    /// samples for the extent search. Returns the validated settings so the
    /// solver can chain the check before running.
    fn validate(self) -> Result<Self, ReactionExtentError> {
        if self.max_iterations == 0 {
            return Err(invalid_problem("max_iterations must be strictly positive"));
        }

        if self.bracket_samples < 2 {
            return Err(invalid_problem("bracket_samples must be at least 2"));
        }

        if !self.feasibility_margin.is_finite()
            || self.feasibility_margin <= 0.0
            || self.feasibility_margin >= 0.5
        {
            return Err(invalid_problem(
                "feasibility_margin must lie strictly between 0 and 0.5",
            ));
        }

        Ok(self)
    }
}

/// Finite two-phase equilibrium predicted by the independent K_eq route.
#[derive(Debug, Clone, PartialEq)]
pub struct PurePhaseEquilibriumResult {
    /// Reaction extent from the gas-only boundary state.
    pub extent: f64,

    /// Gas composition at the independent two-phase equilibrium.
    pub gas_moles: Vec<f64>,

    /// Amount of the pure condensed phase.
    pub candidate_moles: f64,

    /// Final `ln(Q) - ln(K)` residual.
    pub log_residual: f64,

    /// Number of scalar bisection iterations.
    pub iterations: usize,
}

/// High-level result of the independent validator.
///
/// A stable inactive phase has no finite two-phase solve because the
/// phase-forming direction is uphill already at ξ = 0.
/// Structural and state identity of a pure-phase validation case.
///
/// The bit-level representation deliberately makes independent evidence
/// applicable only to the exact state it was evaluated for. This prevents a
/// comparison between equal-length but differently ordered species vectors or
/// between adjacent temperature points.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PurePhaseBoundaryCaseIdentity {
    /// Gas species names in declared order.
    pub gas_species: Vec<String>,
    /// Name of the pure condensed candidate phase.
    pub candidate_name: String,
    /// Exact bit pattern of each gas mole number at absence.
    pub gas_mole_bits: Vec<u64>,
    /// Exact bit pattern of each gas stoichiometric coefficient.
    pub gas_stoichiometry_bits: Vec<u64>,
    /// Exact bit pattern of the candidate stoichiometric coefficient.
    pub candidate_stoichiometry_bits: u64,
    /// Exact bit pattern of the temperature.
    pub temperature_bits: u64,
    /// Exact bit pattern of the pressure.
    pub pressure_bits: u64,
    /// Exact bit pattern of the reference pressure.
    pub reference_pressure_bits: u64,
}

/// High-level result of the independent validator.
#[derive(Debug, Clone, PartialEq)]
pub struct PurePhaseValidationResult {
    /// Identity of the problem this result was produced for.
    pub case_identity: PurePhaseBoundaryCaseIdentity,

    /// Boundary tolerances used to classify the independent result.
    pub boundary_tolerances: PurePhaseBoundaryTolerances,

    /// Zero-extent boundary evidence and classification.
    pub boundary: PurePhaseBoundaryReport,
    /// Finite two-phase result, present only when formation is favorable.
    pub equilibrium: Option<PurePhaseEquilibriumResult>,
}

/// Independent validator.
///
/// The implementation intentionally uses a very conservative scalar method.
/// Performance is irrelevant here; independence and robustness are more
/// valuable than sharing the production nonlinear machinery.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PurePhaseBoundaryValidator {
    /// Boundary and finite-equilibrium acceptance tolerances.
    pub tolerances: PurePhaseBoundaryTolerances,
    /// Numerical settings for the scalar extent solve.
    pub solver_settings: PurePhaseBoundarySolverSettings,
}

impl Default for PurePhaseBoundaryValidator {
    fn default() -> Self {
        Self {
            tolerances: PurePhaseBoundaryTolerances::default(),
            solver_settings: PurePhaseBoundarySolverSettings::default(),
        }
    }
}

impl PurePhaseBoundaryValidator {
    pub fn validate(
        &self,
        problem: &PurePhaseBoundaryProblem,
    ) -> Result<PurePhaseValidationResult, ReactionExtentError> {
        let boundary = evaluate_pure_phase_boundary(problem, self.tolerances)?;

        match boundary.prediction {
            PurePhaseBoundaryPrediction::StableInactive | PurePhaseBoundaryPrediction::Boundary => {
                Ok(PurePhaseValidationResult {
                    case_identity: problem.case_identity(),
                    boundary_tolerances: self.tolerances,
                    boundary,
                    equilibrium: None,
                })
            }

            PurePhaseBoundaryPrediction::ShouldAppear => {
                let equilibrium = self.solve_two_phase_equilibrium(problem)?;

                Ok(PurePhaseValidationResult {
                    case_identity: problem.case_identity(),
                    boundary_tolerances: self.tolerances,
                    boundary,
                    equilibrium: Some(equilibrium),
                })
            }
        }
    }

    /// Solves `ln(Q(ξ)) - ln(K) = 0` for ξ > 0.
    ///
    /// Since the boundary residual is already known to be negative, the
    /// validator searches for a positive extent at which the residual becomes
    /// positive. Once such a bracket is found, plain bisection is enough.
    ///
    /// Deliberately NOT using Newton here gives an implementation that is
    /// structurally different from both the production equilibrium solver and
    /// the existing independent K_eq safeguarded Newton solver.
    fn solve_two_phase_equilibrium(
        &self,
        problem: &PurePhaseBoundaryProblem,
    ) -> Result<PurePhaseEquilibriumResult, ReactionExtentError> {
        let settings = self.solver_settings.validate()?;
        let tolerances = self.tolerances.validate()?;

        let raw_upper = problem.positive_extent_upper_bound()?;
        let normalized_upper = 1.0 - settings.feasibility_margin;
        let upper = raw_upper * normalized_upper;

        if !upper.is_finite() || upper <= 0.0 {
            return Err(invalid_problem(
                "positive feasible extent interval collapsed",
            ));
        }

        let f0 = problem.log_residual_at_extent(0.0)?;

        if f0 >= 0.0 {
            return Err(invalid_problem(
                "finite two-phase solve requested although candidate formation is not favorable at zero extent",
            ));
        }

        // Search deterministically for the first sign-changing interval.
        // A linear scan is intentionally simple and transparent.
        // Work in eta = extent / raw_upper.  The residual itself still uses
        // physical moles, but bisection geometry must not acquire an implicit
        // one-mole scale through its termination condition.
        let mut left = 0.0;
        let mut f_left = f0;
        let mut bracket = None;

        for sample in 1..=settings.bracket_samples {
            let fraction = sample as f64 / settings.bracket_samples as f64;
            let right = normalized_upper * fraction;
            let f_right = problem.log_residual_at_extent(raw_upper * right)?;

            if f_right.abs() <= tolerances.max_abs_equilibrium_log_residual {
                return self.accept(problem, raw_upper * right, f_right, 0);
            }

            if f_left.signum() != f_right.signum() {
                bracket = Some((left, right, f_left, f_right));
                break;
            }

            left = right;
            f_left = f_right;
        }

        let (mut left, mut right, mut f_left, mut f_right) =
            bracket.ok_or_else(|| ReactionExtentError::ValidationNotApplicable {
                path: "pure_phase_boundary_validator",
                message: concat!(
                    "candidate formation is favorable at zero extent, ",
                    "but no interior two-phase equilibrium root was bracketed ",
                    "before a gas-species positivity boundary"
                )
                .to_string(),
            })?;

        for iteration in 0..settings.max_iterations {
            let mid = 0.5 * (left + right);
            let f_mid = problem.log_residual_at_extent(raw_upper * mid)?;

            if f_mid.abs() <= tolerances.max_abs_equilibrium_log_residual {
                return self.accept(problem, raw_upper * mid, f_mid, iteration + 1);
            }

            // Standard bracket update. We do not assume a particular
            // monotonic direction; only the sign change matters.
            if f_left.signum() != f_mid.signum() {
                right = mid;
                f_right = f_mid;
            } else {
                left = mid;
                f_left = f_mid;
            }

            // eta is dimensionless, so this is a scale-free floating-point
            // resolution test for the physical extent interval.
            if right <= left
                || (right - left) <= f64::EPSILON * left.abs().max(right.abs()).max(1.0)
            {
                let (extent, residual) = if f_left.abs() <= f_right.abs() {
                    (raw_upper * left, f_left)
                } else {
                    (raw_upper * right, f_right)
                };

                if residual.abs() <= tolerances.max_abs_equilibrium_log_residual {
                    return self.accept(problem, extent, residual, iteration + 1);
                }

                break;
            }
        }

        Err(ReactionExtentError::SolveError(
            crate::Thermodynamics::ChemEquilibrium::equilibrium_nonlinear::SolveError::MaxIterations,
        ))
    }

    /// Materializes the accepted two-phase equilibrium result from a solved
    /// extent and its log-residual.
    ///
    /// Reconstructs the physical gas and candidate mole numbers implied by the
    /// converged extent and packages them, together with the residual and
    /// iteration count, into an immutable [`PurePhaseEquilibriumResult`]. This
    /// is a pure projection; it performs no further iteration.
    fn accept(
        &self,
        problem: &PurePhaseBoundaryProblem,
        extent: f64,
        log_residual: f64,
        iterations: usize,
    ) -> Result<PurePhaseEquilibriumResult, ReactionExtentError> {
        let gas_moles = problem.gas_moles_at_extent(extent)?;
        let candidate_moles = problem.candidate_moles_at_extent(extent)?;

        Ok(PurePhaseEquilibriumResult {
            extent,
            gas_moles,
            candidate_moles,
            log_residual,
            iterations,
        })
    }
}

// ============================================================================
// Cross-validation against the canonical phase-control result
// ============================================================================

/// Minimal canonical evidence needed by this independent validator.
///
/// This deliberately avoids depending directly on the complete production
/// outcome type. An adapter can be written next to the production runner.
///
/// `candidate_active` is the final physical phase state, not merely presence
/// in a numerical recovery probe.
#[derive(Debug, Clone, PartialEq)]
pub struct CanonicalPurePhaseEvidence {
    /// Ordered gas component identities matching `gas_moles`.
    ///
    /// The cross-validator rejects a result with a different order even where
    /// dimensions happen to agree. Values without component identities are not
    /// valid cross-validation evidence.
    pub gas_species: Vec<String>,

    /// Identity of the pure candidate represented by the phase state and
    /// amount below.
    pub candidate_name: String,

    /// Whether the candidate is active in the accepted production solution.
    pub candidate_active: bool,

    /// Final physical gas moles, in `PurePhaseBoundaryProblem::gas_species`
    /// order. Entries may be zero after complete condensation; the separate
    /// `boundary_gas_moles` witness remains strictly positive.
    pub gas_moles: Vec<f64>,

    /// Gas-only physical moles at the state where canonical TPD was evaluated.
    ///
    /// This is usually equal to `gas_moles` for a stable inactive candidate.
    /// After phase activation it instead comes from the accepted pre-activation
    /// restart state, while `gas_moles` remains the final two-phase composition.
    /// Keeping both vectors prevents the phase-boundary driving force from
    /// being compared at one state and composition at another.
    pub boundary_gas_moles: Vec<f64>,

    /// Final physical amount of the candidate phase.
    pub candidate_moles: f64,

    /// Canonical TPD value at the gas-only boundary state, when available.
    ///
    /// This field is diagnostic only. The independent validator never uses it
    /// to construct its own prediction.
    pub boundary_minimum_tpd: Option<f64>,
}

/// Extracts the narrow canonical evidence needed by this validator from one
/// genuine pure-phase TPD report.
///
/// The `PhaseStabilityReport` must belong to the requested candidate phase,
/// represent an evaluated one-component condensed candidate, and contain a
/// numeric TPD minimum. `candidate_active` and the final mole amounts belong
/// to the later accepted equilibrium state; they are intentionally separate
/// from the gas-only boundary state at which the TPD report was evaluated.
/// This preserves the order of the physical argument instead of pretending
/// that phase appearance was known before the stability test.
pub fn canonical_pure_phase_evidence_from_stability_report(
    candidate_phase_index: usize,
    gas_species: Vec<String>,
    candidate_name: impl Into<String>,
    candidate_active: bool,
    gas_moles: Vec<f64>,
    candidate_moles: f64,
    stability: &PhaseStabilityReport,
) -> Result<CanonicalPurePhaseEvidence, ReactionExtentError> {
    let candidate_name = candidate_name.into();
    validate_component_layout(&gas_species, &candidate_name)?;

    if stability.phase.index() != candidate_phase_index {
        return Err(invalid_problem(format!(
            "TPD report belongs to phase {}, not requested candidate phase {candidate_phase_index}",
            stability.phase.index()
        )));
    }
    if stability.status != PhaseStabilityStatus::Evaluated {
        return Err(ReactionExtentError::ValidationNotApplicable {
            path: "pure_phase_boundary_canonical_adapter",
            message: format!(
                "candidate phase {candidate_phase_index} has TPD status {:?}, not Evaluated",
                stability.status
            ),
        });
    }
    if stability.layout.phase_component_indices.len() != 1 {
        return Err(ReactionExtentError::ValidationNotApplicable {
            path: "pure_phase_boundary_canonical_adapter",
            message: format!(
                "candidate phase {candidate_phase_index} has {} components; pure-phase validation requires exactly one",
                stability.layout.phase_component_indices.len()
            ),
        });
    }
    let minimum_tpd =
        stability
            .minimum_tpd
            .ok_or_else(|| ReactionExtentError::ValidationNotApplicable {
                path: "pure_phase_boundary_canonical_adapter",
                message: format!(
                    "candidate phase {candidate_phase_index} was evaluated without a TPD minimum"
                ),
            })?;
    if !minimum_tpd.is_finite()
        || gas_species.len() != gas_moles.len()
        || gas_moles
            .iter()
            .any(|moles| !moles.is_finite() || *moles <= 0.0)
        || !candidate_moles.is_finite()
        || candidate_moles < 0.0
    {
        return Err(invalid_candidate(
            "canonical pure-phase evidence contains non-finite or infeasible mole/TPD values",
        ));
    }

    Ok(CanonicalPurePhaseEvidence {
        gas_species,
        candidate_name,
        candidate_active,
        boundary_gas_moles: gas_moles.clone(),
        gas_moles,
        candidate_moles,
        boundary_minimum_tpd: Some(minimum_tpd),
    })
}

/// Comparison tolerances for the canonical-vs-K_eq phase-control check.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PurePhaseCrossValidationTolerances {
    /// Absolute part of the accepted per-species gas mole delta.
    ///
    /// This protects comparisons near zero, where a relative tolerance alone
    /// would be undefined or unhelpfully strict.
    pub max_abs_gas_mole_delta: f64,
    /// Relative part of the accepted per-species gas mole delta.
    ///
    /// A comparison passes when `|a-b| <= abs + rel * max(|a|, |b|)`.
    pub max_relative_gas_mole_delta: f64,

    /// Absolute part of the accepted candidate mole delta.
    pub max_abs_candidate_mole_delta: f64,
    /// Relative part of the accepted candidate mole delta.
    pub max_relative_candidate_mole_delta: f64,

    /// Optional consistency check between two *different units* of the same
    /// infinitesimal driving force:
    ///
    ///     Δ_rG = RT [lnQ-lnK].
    ///
    /// For the q=1 candidate this should match canonical pure-phase TPD up to
    /// reaction normalization. Leave this comparison out unless reaction
    /// normalization is controlled explicitly.
    pub max_abs_boundary_energy_delta: f64,
}

impl Default for PurePhaseCrossValidationTolerances {
    fn default() -> Self {
        Self {
            max_abs_gas_mole_delta: 1e-5,
            max_relative_gas_mole_delta: 1e-8,
            max_abs_candidate_mole_delta: 1e-5,
            max_relative_candidate_mole_delta: 1e-8,
            max_abs_boundary_energy_delta: 1e-5,
        }
    }
}

/// Result of comparing the canonical TPD/outer-loop path with the independent
/// equilibrium-constant path.
#[derive(Debug, Clone, PartialEq)]
pub struct PurePhaseCrossValidationReport {
    /// Whether the candidate is active in the canonical production solution.
    pub canonical_active: bool,
    /// Independent prediction from the zero-extent boundary residual.
    pub independent_prediction: PurePhaseBoundaryPrediction,
    /// Strongest topology statement implied by the independent prediction.
    pub topology_expectation: PurePhaseTopologyExpectation,

    /// `Some` only where independent thermodynamics implies one discrete
    /// topology. `None` means an exact boundary, whose state belongs to the
    /// production hysteresis/history test rather than this algebraic checker.
    pub topology_agreement: Option<bool>,

    /// Agreement between canonical TPD and independent `RT[ln(Q)-ln(K)]/nu`.
    /// `None` means that this comparison was not requested because no canonical
    /// TPD value was supplied.
    pub thermodynamic_agreement: Option<bool>,

    /// Agreement of finite two-phase compositions. `None` is legitimate for a
    /// stable inactive candidate and for an exact boundary, where no unique
    /// finite two-phase reference state is implied.
    pub composition_agreement: Option<bool>,

    /// Largest per-species gas mole delta between canonical and independent.
    pub max_abs_gas_mole_delta: Option<f64>,
    /// Absolute candidate mole delta between canonical and independent.
    pub abs_candidate_mole_delta: Option<f64>,

    /// Difference between canonical pure-phase TPD and the reaction-based
    /// boundary driving force, when that comparison is meaningful.
    pub abs_boundary_energy_delta: Option<f64>,

    /// Completeness-aware result of the comparison.
    pub status: PurePhaseCrossValidationStatus,

    /// Convenience summary for strict release gates. It is true only for
    /// complete evidence, never merely because no available axis disagreed.
    pub accepted: bool,
}

/// Completeness-aware outcome of a pure-phase cross-validation comparison.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PurePhaseCrossValidationStatus {
    /// Every evidence axis required by the independent prediction agreed.
    Complete,
    /// Available evidence agreed, but a required axis was absent.
    ConsistentButPartial,
    /// No meaningful axis could be compared.
    InsufficientEvidence,
    /// At least one independently comparable axis disagreed.
    Disagreed,
}

/// Runs the independent validator and compares its result with canonical
/// evidence for the same explicitly identified case.
pub fn cross_validate_pure_phase(
    problem: &PurePhaseBoundaryProblem,
    canonical: &CanonicalPurePhaseEvidence,
    validator: PurePhaseBoundaryValidator,
    tolerances: PurePhaseCrossValidationTolerances,
) -> Result<PurePhaseCrossValidationReport, ReactionExtentError> {
    let independent = validator.validate(problem)?;
    compare_pure_phase_validation(problem, canonical, &independent, tolerances)
}

/// Compares canonical evidence with an already computed independent result.
///
/// This lower-level function is useful for focused regression tests. It still
/// verifies that the supplied result belongs to `problem`, so callers cannot
/// compare equal-length vectors from another temperature, candidate, or
/// component ordering.
///
/// Important:
/// reaction normalization matters if comparing TPD [J/mol candidate] directly
/// to Δ_rG [J/mol reaction]. For arbitrary `ν_s`, the comparable quantity is
///
///     Δ_rG / ν_s.
///
/// Therefore the code normalizes the independent reaction driving force per
/// mole of candidate phase before comparing it with canonical TPD.
pub(crate) fn compare_pure_phase_validation(
    problem: &PurePhaseBoundaryProblem,
    canonical: &CanonicalPurePhaseEvidence,
    independent: &PurePhaseValidationResult,
    tolerances: PurePhaseCrossValidationTolerances,
) -> Result<PurePhaseCrossValidationReport, ReactionExtentError> {
    validate_cross_tolerances(tolerances)?;

    if independent.case_identity != problem.case_identity() {
        return Err(invalid_problem(
            "independent validation result belongs to a different pure-phase boundary case",
        ));
    }
    let recomputed_boundary =
        evaluate_pure_phase_boundary(problem, independent.boundary_tolerances)?;
    if recomputed_boundary != independent.boundary {
        return Err(invalid_problem(
            "independent validation boundary does not match the supplied pure-phase problem",
        ));
    }

    validate_component_layout(&canonical.gas_species, &canonical.candidate_name)?;

    if canonical.gas_species != problem.gas_species {
        return Err(invalid_problem(
            "canonical gas component identities do not match the phase-boundary problem order",
        ));
    }
    if canonical.candidate_name != problem.candidate_name {
        return Err(invalid_problem(
            "canonical candidate identity does not match the phase-boundary problem",
        ));
    }

    if canonical.gas_moles.len() != problem.gas_species.len() {
        return Err(ReactionExtentError::DimensionMismatch(format!(
            "canonical result has {} gas species but phase-boundary problem has {}",
            canonical.gas_moles.len(),
            problem.gas_species.len()
        )));
    }
    if canonical.boundary_gas_moles.len() != problem.gas_species.len() {
        return Err(ReactionExtentError::DimensionMismatch(format!(
            "canonical boundary state has {} gas species but phase-boundary problem has {}",
            canonical.boundary_gas_moles.len(),
            problem.gas_species.len()
        )));
    }
    if canonical
        .boundary_gas_moles
        .iter()
        .any(|moles| !moles.is_finite() || *moles <= 0.0)
    {
        return Err(invalid_problem(
            "canonical boundary gas moles must be finite and strictly positive",
        ));
    }
    if canonical
        .boundary_gas_moles
        .iter()
        .zip(problem.gas_moles_at_absence())
        .any(|(canonical_moles, independent_moles)| {
            !within_abs_relative(
                *canonical_moles,
                *independent_moles,
                tolerances.max_abs_gas_mole_delta,
                tolerances.max_relative_gas_mole_delta,
            )
        })
    {
        return Err(invalid_problem(
            "canonical TPD boundary state does not match the independent boundary inventory",
        ));
    }

    let topology_expectation = independent.boundary.prediction.topology_expectation();
    let topology_agreement = match topology_expectation {
        PurePhaseTopologyExpectation::MustBeInactive => Some(!canonical.candidate_active),
        PurePhaseTopologyExpectation::HysteresisDependent => None,
        PurePhaseTopologyExpectation::MustBeActive => Some(canonical.candidate_active),
    };

    let mut max_abs_gas_mole_delta = None;
    let mut abs_candidate_mole_delta = None;
    let mut gas_composition_agreement = None;
    let mut candidate_composition_agreement = None;

    // Composition comparison is meaningful only when the independent path
    // actually found a finite two-phase equilibrium.
    if let Some(eq) = &independent.equilibrium {
        let max_gas_delta = canonical
            .gas_moles
            .iter()
            .zip(&eq.gas_moles)
            .map(|(a, b)| (a - b).abs())
            .fold(0.0_f64, f64::max);

        max_abs_gas_mole_delta = Some(max_gas_delta);
        abs_candidate_mole_delta = Some((canonical.candidate_moles - eq.candidate_moles).abs());
        gas_composition_agreement = Some(canonical.gas_moles.iter().zip(&eq.gas_moles).all(
            |(canonical_moles, independent_moles)| {
                within_abs_relative(
                    *canonical_moles,
                    *independent_moles,
                    tolerances.max_abs_gas_mole_delta,
                    tolerances.max_relative_gas_mole_delta,
                )
            },
        ));
        candidate_composition_agreement = Some(within_abs_relative(
            canonical.candidate_moles,
            eq.candidate_moles,
            tolerances.max_abs_candidate_mole_delta,
            tolerances.max_relative_candidate_mole_delta,
        ));
    }

    let abs_boundary_energy_delta = canonical.boundary_minimum_tpd.map(|canonical_tpd| {
        // Δ_rG is per mole of reaction extent. Divide by ν_s to obtain
        // the driving force per mole of candidate phase.
        let independent_per_mole_candidate =
            independent.boundary.reaction_gibbs_at_absence / problem.candidate_stoichiometry();

        (canonical_tpd - independent_per_mole_candidate).abs()
    });

    let composition_agreement = match (gas_composition_agreement, candidate_composition_agreement) {
        (Some(gas), Some(candidate)) => Some(gas && candidate),
        // No finite two-phase equilibrium exists for a stable inactive phase,
        // and an exact boundary does not prescribe one. Both cases are not an
        // omitted or failed composition check.
        (None, None) => None,
        _ => Some(false),
    };

    let thermodynamic_agreement =
        abs_boundary_energy_delta.map(|delta| delta <= tolerances.max_abs_boundary_energy_delta);

    let axes = [
        topology_agreement,
        thermodynamic_agreement,
        composition_agreement,
    ];
    let status = if axes.iter().any(|axis| matches!(axis, Some(false))) {
        PurePhaseCrossValidationStatus::Disagreed
    } else {
        let required_axes_present = match topology_expectation {
            PurePhaseTopologyExpectation::MustBeInactive => {
                topology_agreement == Some(true) && thermodynamic_agreement == Some(true)
            }
            PurePhaseTopologyExpectation::HysteresisDependent => {
                thermodynamic_agreement == Some(true)
            }
            PurePhaseTopologyExpectation::MustBeActive => {
                topology_agreement == Some(true)
                    && thermodynamic_agreement == Some(true)
                    && composition_agreement == Some(true)
            }
        };

        if required_axes_present {
            PurePhaseCrossValidationStatus::Complete
        } else if axes.iter().any(|axis| *axis == Some(true)) {
            PurePhaseCrossValidationStatus::ConsistentButPartial
        } else {
            PurePhaseCrossValidationStatus::InsufficientEvidence
        }
    };
    let accepted = status == PurePhaseCrossValidationStatus::Complete;

    Ok(PurePhaseCrossValidationReport {
        canonical_active: canonical.candidate_active,
        independent_prediction: independent.boundary.prediction,
        topology_expectation,
        topology_agreement,
        thermodynamic_agreement,
        composition_agreement,
        max_abs_gas_mole_delta,
        abs_candidate_mole_delta,
        abs_boundary_energy_delta,
        status,
        accepted,
    })
}

// ============================================================================
// Small report formatting
// ============================================================================

impl fmt::Display for PurePhaseBoundaryPrediction {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::StableInactive => write!(f, "stable-inactive"),
            Self::Boundary => write!(f, "boundary"),
            Self::ShouldAppear => write!(f, "should-appear"),
        }
    }
}

impl fmt::Display for PurePhaseBoundaryReport {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "[boundary] candidate = {}", self.candidate_name)?;
        writeln!(f, "[boundary] T = {:.6}", self.conditions.temperature())?;
        writeln!(f, "[boundary] P = {:.6}", self.conditions.pressure())?;
        writeln!(f, "[boundary] lnQ = {:.6e}", self.ln_q_at_absence)?;
        writeln!(f, "[boundary] lnK = {:.6e}", self.ln_k)?;
        writeln!(
            f,
            "[boundary] lnQ-lnK = {:.6e}",
            self.log_residual_at_absence
        )?;
        writeln!(
            f,
            "[boundary] delta_r_G = {:.6e} J/mol-reaction",
            self.reaction_gibbs_at_absence
        )?;
        writeln!(f, "[boundary] prediction = {}", self.prediction)
    }
}

// ============================================================================
// Helpers
// ============================================================================

/// Computes an SVD rank without borrowing the canonical reaction-basis
/// implementation. The threshold scales with the largest singular value while
/// retaining an absolute floor for an all-small synthetic matrix.
pub(crate) fn independent_matrix_rank(
    matrix: &DMatrix<f64>,
    tolerances: PurePhaseBoundaryStructuralTolerances,
) -> Result<usize, ReactionExtentError> {
    if matrix.nrows() == 0 || matrix.ncols() == 0 {
        return Err(invalid_problem(
            "independent reaction-space matrix must have non-zero dimensions",
        ));
    }
    if matrix.iter().any(|value| !value.is_finite()) {
        return Err(invalid_problem(
            "independent reaction-space matrix contains a non-finite value",
        ));
    }
    let singular_values = nalgebra::linalg::SVD::new(matrix.clone(), false, false).singular_values;
    let largest = singular_values.iter().copied().fold(0.0_f64, f64::max);
    let threshold = tolerances
        .rank_absolute_tolerance
        .max(tolerances.rank_relative_tolerance * largest);
    Ok(singular_values
        .iter()
        .filter(|&&value| value > threshold)
        .count())
}

/// Returns whether two values satisfy the shared absolute-plus-relative
/// comparison contract.
///
/// The absolute term keeps zero and trace inventories meaningful; the
/// relative term keeps an otherwise identical scaled system from failing only
/// because its physical inventory is larger.
fn within_abs_relative(actual: f64, expected: f64, absolute: f64, relative: f64) -> bool {
    (actual - expected).abs() <= absolute + relative * actual.abs().max(expected.abs())
}

/// Rejects non-finite or invalid cross-validation tolerances.
///
/// Shared by the public and internal cross-validation entry points so both use
/// the same acceptance contract for mole, energy, and boundary disagreements.
fn validate_cross_tolerances(
    tolerances: PurePhaseCrossValidationTolerances,
) -> Result<(), ReactionExtentError> {
    for (name, value) in [
        ("max_abs_gas_mole_delta", tolerances.max_abs_gas_mole_delta),
        (
            "max_abs_candidate_mole_delta",
            tolerances.max_abs_candidate_mole_delta,
        ),
        (
            "max_abs_boundary_energy_delta",
            tolerances.max_abs_boundary_energy_delta,
        ),
    ] {
        if !value.is_finite() || value <= 0.0 {
            return Err(invalid_problem(format!(
                "{name} must be finite and strictly positive"
            )));
        }
    }

    for (name, value) in [
        (
            "max_relative_gas_mole_delta",
            tolerances.max_relative_gas_mole_delta,
        ),
        (
            "max_relative_candidate_mole_delta",
            tolerances.max_relative_candidate_mole_delta,
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

/// Validates the named component layout carried by canonical evidence.
///
/// Checks that the gas-species names and candidate name in the canonical
/// evidence exactly match the problem's declared identities, so cross-validation
/// never compares two compositions that merely share a vector length.
fn validate_component_layout(
    gas_species: &[String],
    candidate_name: &str,
) -> Result<(), ReactionExtentError> {
    if gas_species.is_empty() {
        return Err(invalid_problem(
            "canonical pure-phase evidence requires at least one gas component identity",
        ));
    }
    if gas_species.iter().any(|species| species.trim().is_empty()) {
        return Err(invalid_problem(
            "canonical pure-phase evidence gas component identities must be non-empty",
        ));
    }
    let mut unique_species = std::collections::BTreeSet::new();
    if gas_species
        .iter()
        .any(|species| !unique_species.insert(species.as_str()))
    {
        return Err(invalid_problem(
            "canonical pure-phase evidence gas component identities must be unique",
        ));
    }
    if candidate_name.trim().is_empty() {
        return Err(invalid_problem(
            "canonical pure-phase evidence candidate identity must be non-empty",
        ));
    }

    Ok(())
}

fn invalid_problem(message: impl Into<String>) -> ReactionExtentError {
    ReactionExtentError::InvalidProblem {
        field: "pure_phase_boundary_validator",
        message: message.into(),
    }
}

fn invalid_candidate(message: impl Into<String>) -> ReactionExtentError {
    ReactionExtentError::InvalidCandidate {
        field: "pure_phase_boundary_validator",
        message: message.into(),
    }
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use std::rc::Rc;

    fn constant_gibbs(value: f64) -> GibbsFn {
        Rc::new(move |_| value)
    }

    /// Builds a completely synthetic reaction
    ///
    ///     2 A(g) -> B(g) + S(cond)
    ///
    /// where `S` is the pure candidate phase.
    ///
    /// Nothing in the validator knows what A, B, or S physically are. A
    /// Boudouard fixture can later use the same machinery with
    ///
    ///     A = CO,
    ///     B = CO2,
    ///     S = C(graphite).
    fn synthetic_problem_for_target_k(
        target_k: f64,
        gas_moles: [f64; 2],
    ) -> PurePhaseBoundaryProblem {
        synthetic_problem_for_target_k_with_stoichiometric_scale(target_k, gas_moles, 1.0)
    }

    /// Builds the same physical reaction under a different reaction-coordinate
    /// normalization. `target_k` belongs to the unscaled reaction. The
    /// stoichiometric vector is scaled, while standard Gibbs functions remain
    /// material properties and therefore must not be scaled a second time.
    fn synthetic_problem_for_target_k_with_stoichiometric_scale(
        target_k: f64,
        gas_moles: [f64; 2],
        stoichiometric_scale: f64,
    ) -> PurePhaseBoundaryProblem {
        let temperature = 1_000.0;
        let pressure = 101_325.0;
        let reference_pressure = 101_325.0;

        // Reaction:
        //
        //     -2 A + B + S = 0
        //
        // Choose standard Gibbs functions so that the complete reaction has
        // exactly the requested K:
        //
        //     ΔG° = -RT ln K.
        //
        // For simplicity assign the entire ΔG° to the pure candidate.
        let delta_g = -MOLAR_GAS_CONSTANT * temperature * target_k.ln();

        PurePhaseBoundaryProblem::new(
            vec!["A".into(), "B".into()],
            gas_moles.to_vec(),
            vec![-2.0 * stoichiometric_scale, stoichiometric_scale],
            stoichiometric_scale,
            vec![constant_gibbs(0.0), constant_gibbs(0.0)],
            constant_gibbs(delta_g),
            EquilibriumConditions::new(temperature, pressure, reference_pressure).unwrap(),
            "S",
        )
        .unwrap()
    }

    fn strict_independent_fixture() -> PurePhaseBoundaryProblem {
        // The synthetic reaction is -2 A + B + S = 0. Give the species two
        // abstract elements such that
        //
        //     -2 [1, 1] + [1, 0] + [1, 2] = [0, 0].
        //
        // The two gas rows are independent, so gas-only composition has no
        // reaction degree of freedom. Adding S gives three species over a
        // rank-two matrix, hence exactly one phase-forming direction.
        let composition = PurePhaseBoundaryElementComposition::new(
            vec!["X".into(), "Y".into()],
            DMatrix::from_row_slice(2, 2, &[1.0, 1.0, 1.0, 0.0]),
            vec![1.0, 2.0],
        )
        .expect("synthetic composition is structurally valid");
        synthetic_problem_for_target_k(4.0, [1.0, 1.0])
            .with_element_composition(
                composition,
                PurePhaseBoundaryStructuralTolerances::default(),
            )
            .expect("synthetic phase-forming reaction conserves X and Y")
    }

    /// A controlled temperature fixture with a known boundary. At equal gas
    /// mole numbers and reference pressure, `ln(Q) = ln(2)`. Choosing
    ///
    ///     ln(K(T)) = ln(2) + slope * (T_boundary - T)
    ///
    /// makes `T_boundary` the exact zero of the independent log residual.
    /// The complete temperature dependence lives in the synthetic standard
    /// Gibbs function; no production TPD or phase-control code participates.
    fn analytic_temperature_boundary_problem(
        temperature: f64,
        boundary_temperature: f64,
        slope_per_kelvin: f64,
    ) -> PurePhaseBoundaryProblem {
        let pressure = 101_325.0;
        let ln_q_at_absence = 2.0_f64.ln();
        let candidate_gibbs: GibbsFn = Rc::new(move |evaluated_temperature| {
            let ln_k =
                ln_q_at_absence + slope_per_kelvin * (boundary_temperature - evaluated_temperature);
            -MOLAR_GAS_CONSTANT * evaluated_temperature * ln_k
        });

        let composition = PurePhaseBoundaryElementComposition::new(
            vec!["X".into(), "Y".into()],
            DMatrix::from_row_slice(2, 2, &[1.0, 1.0, 1.0, 0.0]),
            vec![1.0, 2.0],
        )
        .expect("synthetic composition is valid");
        PurePhaseBoundaryProblem::new(
            vec!["A".into(), "B".into()],
            vec![1.0, 1.0],
            vec![-2.0, 1.0],
            1.0,
            vec![constant_gibbs(0.0), constant_gibbs(0.0)],
            candidate_gibbs,
            EquilibriumConditions::new(temperature, pressure, pressure).unwrap(),
            "S",
        )
        .unwrap()
        .with_element_composition(
            composition,
            PurePhaseBoundaryStructuralTolerances::default(),
        )
        .expect("the analytic reaction conserves X and Y")
    }

    #[test]
    fn strict_independent_family_accepts_one_phase_forming_direction() {
        let reaction_space = strict_independent_fixture()
            .validate_strict_independent_family(PurePhaseBoundaryStructuralTolerances::default())
            .expect("the fixture has exactly one full phase-forming direction");

        assert_eq!(reaction_space.full_rank, 2);
        assert_eq!(reaction_space.gas_only_rank, 2);
        assert_eq!(reaction_space.full_reaction_dimension, 1);
        assert_eq!(reaction_space.gas_only_reaction_dimension, 0);
        assert_eq!(
            reaction_space.full_reaction_dimension,
            reaction_space.gas_only_reaction_dimension + 1,
            "adding one pure candidate contributes exactly one phase-forming direction"
        );
        assert!(
            reaction_space
                .element_balance_residuals
                .iter()
                .all(|residual| residual.abs() < 1e-12),
            "the independently supplied reaction must conserve every element"
        );
    }

    #[test]
    fn element_aware_problem_rejects_nonconserving_reaction_without_projection() {
        // Replacing S=[1,2] with S=[1,1] leaves a Y residual of -1. The
        // validator must reject this fixture rather than silently repairing
        // its independently supplied reaction.
        let composition = PurePhaseBoundaryElementComposition::new(
            vec!["X".into(), "Y".into()],
            DMatrix::from_row_slice(2, 2, &[1.0, 1.0, 1.0, 0.0]),
            vec![1.0, 1.0],
        )
        .unwrap();
        let result = synthetic_problem_for_target_k(4.0, [1.0, 1.0]).with_element_composition(
            composition,
            PurePhaseBoundaryStructuralTolerances::default(),
        );
        let error = match result {
            Ok(_) => panic!("unbalanced reaction must not become a validation problem"),
            Err(error) => error,
        };

        assert!(matches!(error, ReactionExtentError::InvalidProblem { .. }));
        assert!(
            error
                .to_string()
                .contains("violates elemental conservation")
        );
    }

    #[test]
    fn strict_family_reports_excess_full_reaction_dimension() {
        // All three species carry one abstract element. The supplied reaction
        // is conserved, but a rank-one full matrix over three species leaves
        // two independent reaction directions, so it cannot be used as the
        // strict one-reaction cross-validation family.
        let composition = PurePhaseBoundaryElementComposition::new(
            vec!["X".into()],
            DMatrix::from_row_slice(2, 1, &[1.0, 1.0]),
            vec![1.0],
        )
        .unwrap();
        let problem = synthetic_problem_for_target_k(4.0, [1.0, 1.0])
            .with_element_composition(
                composition,
                PurePhaseBoundaryStructuralTolerances::default(),
            )
            .unwrap();
        let reaction_space = problem
            .reaction_space(PurePhaseBoundaryStructuralTolerances::default())
            .unwrap();

        assert_eq!(reaction_space.full_reaction_dimension, 2);
        let error = problem
            .validate_strict_independent_family(PurePhaseBoundaryStructuralTolerances::default())
            .expect_err("multiple full-system directions invalidate the strict family");
        assert!(
            error
                .to_string()
                .contains("one full-system reaction dimension")
        );
    }

    #[test]
    fn strict_family_reports_remaining_gas_only_chemistry() {
        // The same rank-one gas matrix leaves one gas-only reaction direction
        // when S is absent. This is separately observable in the rank report,
        // even though the strict-family predicate also rejects the enlarged
        // full reaction space first.
        let composition = PurePhaseBoundaryElementComposition::new(
            vec!["X".into()],
            DMatrix::from_row_slice(2, 1, &[1.0, 1.0]),
            vec![1.0],
        )
        .unwrap();
        let problem = synthetic_problem_for_target_k(4.0, [1.0, 1.0])
            .with_element_composition(
                composition,
                PurePhaseBoundaryStructuralTolerances::default(),
            )
            .unwrap();

        let reaction_space = problem
            .reaction_space(PurePhaseBoundaryStructuralTolerances::default())
            .unwrap();
        assert_eq!(reaction_space.gas_only_reaction_dimension, 1);
    }

    #[test]
    fn generic_problem_requires_explicit_composition_for_strict_validation() {
        let error = synthetic_problem_for_target_k(4.0, [1.0, 1.0])
            .validate_strict_independent_family(PurePhaseBoundaryStructuralTolerances::default())
            .expect_err("generic K_eq validation has no independent element matrix");
        assert!(matches!(
            error,
            ReactionExtentError::ValidationNotApplicable { .. }
        ));
    }

    #[test]
    fn independent_temperature_bisection_recovers_analytic_phase_boundary() {
        let expected_temperature = 600.0;
        let slope_per_kelvin = 0.02;
        let tolerances = PurePhaseBoundaryTolerances {
            max_abs_boundary_log_residual: 1e-11,
            ..PurePhaseBoundaryTolerances::default()
        };

        let cold =
            analytic_temperature_boundary_problem(500.0, expected_temperature, slope_per_kelvin);
        let hot =
            analytic_temperature_boundary_problem(700.0, expected_temperature, slope_per_kelvin);
        assert_eq!(
            evaluate_pure_phase_boundary(&cold, tolerances)
                .unwrap()
                .prediction,
            PurePhaseBoundaryPrediction::ShouldAppear
        );
        assert_eq!(
            evaluate_pure_phase_boundary(&hot, tolerances)
                .unwrap()
                .prediction,
            PurePhaseBoundaryPrediction::StableInactive
        );

        let root = bisect_pure_phase_boundary_temperature(
            500.0,
            700.0,
            PurePhaseBoundaryTemperatureSearchSettings::default(),
            tolerances,
            |temperature| {
                Ok(analytic_temperature_boundary_problem(
                    temperature,
                    expected_temperature,
                    slope_per_kelvin,
                ))
            },
        )
        .expect("the known boundary is bracketed by the synthetic temperatures");

        assert!((root.temperature - expected_temperature).abs() < 1e-8);
        assert!(root.log_residual.abs() <= tolerances.max_abs_boundary_log_residual);
        assert!(root.iterations > 0);
    }

    #[test]
    fn boundary_residual_is_ln_q_minus_ln_k() {
        let problem = synthetic_problem_for_target_k(2.0, [1.0, 1.0]);

        let report =
            evaluate_pure_phase_boundary(&problem, PurePhaseBoundaryTolerances::default()).unwrap();

        let x_a: f64 = 0.5;
        let x_b: f64 = 0.5;

        let expected_ln_q = x_b.ln() - 2.0 * x_a.ln();

        assert!((report.ln_q_at_absence - expected_ln_q).abs() < 1e-12);
        assert!((report.ln_k - 2.0_f64.ln()).abs() < 1e-12);
    }

    #[test]
    fn negative_boundary_residual_predicts_phase_appearance() {
        // At 1 bar and n_A = n_B:
        //
        //     Q = x_B / x_A^2 = 0.5 / 0.25 = 2.
        //
        // Select K = 4, hence Q/K = 0.5 and
        //
        //     ln(Q/K) < 0.
        //
        // The forward reaction creating S must therefore be favorable.
        let problem = synthetic_problem_for_target_k(4.0, [1.0, 1.0]);

        let report =
            evaluate_pure_phase_boundary(&problem, PurePhaseBoundaryTolerances::default()).unwrap();

        assert_eq!(report.prediction, PurePhaseBoundaryPrediction::ShouldAppear);
        assert!(report.log_residual_at_absence < 0.0);
        assert!(report.reaction_gibbs_at_absence < 0.0);
    }

    #[test]
    fn positive_boundary_residual_predicts_stable_absence() {
        // Same gas state has Q=2, but now K=1.
        let problem = synthetic_problem_for_target_k(1.0, [1.0, 1.0]);

        let report =
            evaluate_pure_phase_boundary(&problem, PurePhaseBoundaryTolerances::default()).unwrap();

        assert_eq!(
            report.prediction,
            PurePhaseBoundaryPrediction::StableInactive
        );
        assert!(report.log_residual_at_absence > 0.0);
    }

    #[test]
    fn exactly_matched_q_and_k_marks_phase_boundary() {
        // Same gas state gives Q=2.
        let problem = synthetic_problem_for_target_k(2.0, [1.0, 1.0]);

        let report =
            evaluate_pure_phase_boundary(&problem, PurePhaseBoundaryTolerances::default()).unwrap();

        assert_eq!(report.prediction, PurePhaseBoundaryPrediction::Boundary);
        assert!(report.log_residual_at_absence.abs() < 1e-12);
    }

    #[test]
    fn favorable_phase_can_be_solved_to_finite_two_phase_equilibrium() {
        let problem = synthetic_problem_for_target_k(4.0, [1.0, 1.0]);

        let validator = PurePhaseBoundaryValidator::default();
        let result = validator.validate(&problem).unwrap();

        assert_eq!(
            result.boundary.prediction,
            PurePhaseBoundaryPrediction::ShouldAppear
        );

        let equilibrium = result
            .equilibrium
            .expect("favorable candidate should have an interior root");

        assert!(equilibrium.extent > 0.0);
        assert!(equilibrium.candidate_moles > 0.0);
        assert!(
            equilibrium.log_residual.abs() <= validator.tolerances.max_abs_equilibrium_log_residual
        );

        // Conservation along the supplied reaction direction is automatic:
        //
        // n_A = 1 - 2ξ
        // n_B = 1 + ξ
        // n_S = ξ
        assert!((equilibrium.gas_moles[0] - (1.0 - 2.0 * equilibrium.extent)).abs() < 1e-12);

        assert!((equilibrium.gas_moles[1] - (1.0 + equilibrium.extent)).abs() < 1e-12);

        assert!((equilibrium.candidate_moles - equilibrium.extent).abs() < 1e-12);
    }

    #[test]
    fn finite_equilibrium_is_invariant_under_reaction_coordinate_scaling() {
        let validator = PurePhaseBoundaryValidator::default();
        let mut reference: Option<PurePhaseEquilibriumResult> = None;

        for scale in [0.5, 1.0, 2.0] {
            let problem =
                synthetic_problem_for_target_k_with_stoichiometric_scale(4.0, [1.0, 1.0], scale);
            let result = validator
                .validate(&problem)
                .expect("a rescaled representation is still the same physical reaction");
            let equilibrium = result
                .equilibrium
                .expect("candidate formation remains favorable after normalization");

            assert!(
                equilibrium.log_residual.abs()
                    <= validator.tolerances.max_abs_equilibrium_log_residual,
                "scale={scale} must still solve the independent K_eq equation"
            );
            if let Some(reference) = &reference {
                for (scaled, unscaled) in equilibrium.gas_moles.iter().zip(&reference.gas_moles) {
                    assert!(
                        (scaled - unscaled).abs() < 1e-10,
                        "gas state changed under stoichiometric scale {scale}"
                    );
                }
                assert!(
                    (equilibrium.candidate_moles - reference.candidate_moles).abs() < 1e-10,
                    "candidate amount changed under stoichiometric scale {scale}"
                );
            } else {
                reference = Some(equilibrium);
            }
        }
    }

    #[test]
    fn stable_inactive_phase_does_not_run_finite_extent_solver() {
        let problem = synthetic_problem_for_target_k(1.0, [1.0, 1.0]);

        let result = PurePhaseBoundaryValidator::default()
            .validate(&problem)
            .unwrap();

        assert_eq!(
            result.boundary.prediction,
            PurePhaseBoundaryPrediction::StableInactive
        );
        assert!(result.equilibrium.is_none());
    }

    #[test]
    fn reaction_orientation_requires_candidate_to_be_product() {
        let result = PurePhaseBoundaryProblem::new(
            vec!["A".into(), "B".into()],
            vec![1.0, 1.0],
            vec![-2.0, 1.0],
            -1.0, // forbidden orientation
            vec![constant_gibbs(0.0), constant_gibbs(0.0)],
            constant_gibbs(0.0),
            EquilibriumConditions::new(1_000.0, 101_325.0, 101_325.0).unwrap(),
            "S",
        );

        assert!(result.is_err());
    }
}
