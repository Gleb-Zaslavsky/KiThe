//! Pure stability-analysis primitives for the phase-control workflow.
//!
//! The outer workflow owns active-set transitions, hysteresis, cycle detection,
//! and transactional publication. This module owns the immutable thermodynamic
//! state from which phase stability is evaluated: canonical chemical
//! potentials, element-potential reconstruction, and, in later P7 passes, the
//! reacting-system tangent-plane-distance minimization itself.
//!
//! Keeping this boundary pure is deliberate. A failed stability calculation
//! must not mutate an accepted composition or change a `PhaseSet`.

use crate::Thermodynamics::ChemEquilibrium::equilibrium_activity::PhaseActivityModel;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_ids::PhaseIndex;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_log_moles::{GibbsFn, Phase, R};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_nonlinear::ReactionExtentError;
use nalgebra::{
    DMatrix, DVector,
    linalg::{SVD, SymmetricEigen},
};
use std::collections::HashMap;

/// Immutable thermodynamic state used by phase-stability analysis.
///
/// Values are stored in canonical solver order. `chemical_potentials` are
/// evaluated through the same phase activity law as the residual and Jacobian
/// paths, so stability analysis cannot silently use a different standard-state
/// or pressure convention.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct CanonicalPhaseState {
    moles: Vec<f64>,
    phase_totals: Vec<f64>,
    standard_gibbs: Vec<f64>,
    chemical_potentials: Vec<f64>,
    species_phase: Vec<usize>,
    element_composition: DMatrix<f64>,
    active: Vec<bool>,
    candidates: Vec<bool>,
    temperature: f64,
    pressure: f64,
    reference_pressure: f64,
}

impl CanonicalPhaseState {
    /// Reconstructs and validates a stability snapshot from accepted log-moles.
    ///
    /// Every declared phase must have a positive total in this representation.
    /// The outer active-set workflow maintains trace coordinates for absent
    /// phases precisely so the log-mole and activity contracts remain defined.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn from_log_moles(
        log_moles: &[f64],
        gibbs: &[GibbsFn],
        phases: &[Phase],
        species_phase: &[usize],
        element_composition: &DMatrix<f64>,
        temperature: f64,
        pressure: f64,
        reference_pressure: f64,
        active: &[bool],
        candidates: &[bool],
    ) -> Result<Self, ReactionExtentError> {
        let species_count = log_moles.len();
        if gibbs.len() != species_count
            || species_phase.len() != species_count
            || element_composition.nrows() != species_count
            || active.len() != phases.len()
            || candidates.len() != phases.len()
        {
            return Err(ReactionExtentError::DimensionMismatch(format!(
                "phase stability has {species_count} log-moles, {} Gibbs functions, {} phase labels, {} element rows, {} phases, {} active flags, and {} candidate flags",
                gibbs.len(),
                species_phase.len(),
                element_composition.nrows(),
                phases.len(),
                active.len(),
                candidates.len(),
            )));
        }
        for (parameter, value) in [
            ("temperature", temperature),
            ("pressure", pressure),
            ("reference_pressure", reference_pressure),
        ] {
            if !value.is_finite() || value <= 0.0 {
                return Err(ReactionExtentError::InvalidConditions { parameter, value });
            }
        }

        let mut moles = Vec::with_capacity(species_count);
        for (index, &log_moles) in log_moles.iter().enumerate() {
            if !log_moles.is_finite() {
                return Err(ReactionExtentError::InvalidCandidate {
                    field: "phase_stability_log_moles",
                    message: format!(
                        "species {index} has a non-finite log-mole coordinate {log_moles}"
                    ),
                });
            }
            let moles_i = log_moles.exp();
            if !moles_i.is_finite() || moles_i <= 0.0 {
                return Err(ReactionExtentError::InvalidCandidate {
                    field: "phase_stability_moles",
                    message: format!("species {index} reconstructs invalid mole number {moles_i}"),
                });
            }
            moles.push(moles_i);
        }

        let mut phase_totals = vec![0.0; phases.len()];
        for (species, &phase) in species_phase.iter().enumerate() {
            let total = phase_totals.get_mut(phase).ok_or_else(|| {
                ReactionExtentError::DimensionMismatch(format!(
                    "species {species} refers to missing phase {phase}"
                ))
            })?;
            *total += moles[species];
        }
        for (phase, &phase_total) in phase_totals.iter().enumerate() {
            if !phase_total.is_finite() || phase_total <= 0.0 {
                return Err(ReactionExtentError::InvalidNPhase {
                    index: phase,
                    value: phase_total,
                });
            }
        }

        let rt = R * temperature;
        let mut standard_gibbs = Vec::with_capacity(species_count);
        let mut chemical_potentials = Vec::with_capacity(species_count);
        for species in 0..species_count {
            let phase = species_phase[species];
            let g0 = gibbs[species](temperature);
            if !g0.is_finite() {
                return Err(ReactionExtentError::InvalidDG0 {
                    species_index: species,
                    dg0: g0,
                    temperature,
                });
            }
            let log_activity = phases[phase].kind.log_activity(
                moles[species],
                phase_totals[phase],
                pressure,
                reference_pressure,
            )?;
            let chemical_potential = g0 + rt * log_activity;
            if !chemical_potential.is_finite() {
                return Err(ReactionExtentError::InvalidCandidate {
                    field: "phase_stability_chemical_potential",
                    message: format!("species {species} produced a non-finite chemical potential"),
                });
            }
            standard_gibbs.push(g0);
            chemical_potentials.push(chemical_potential);
        }

        Ok(Self {
            moles,
            phase_totals,
            standard_gibbs,
            chemical_potentials,
            species_phase: species_phase.to_vec(),
            element_composition: element_composition.clone(),
            active: active.to_vec(),
            candidates: candidates.to_vec(),
            temperature,
            pressure,
            reference_pressure,
        })
    }

    pub(crate) fn moles(&self) -> &[f64] {
        &self.moles
    }

    pub(crate) fn standard_gibbs(&self) -> &[f64] {
        &self.standard_gibbs
    }

    pub(crate) fn chemical_potentials(&self) -> &[f64] {
        &self.chemical_potentials
    }

    pub(crate) fn species_phase(&self) -> &[usize] {
        &self.species_phase
    }

    pub(crate) fn element_composition(&self) -> &DMatrix<f64> {
        &self.element_composition
    }

    pub(crate) fn active(&self) -> &[bool] {
        &self.active
    }

    pub(crate) fn candidates(&self) -> &[bool] {
        &self.candidates
    }
}

/// Auditable least-squares reconstruction of elemental chemical potentials.
///
/// A rank-deficient reference matrix is valid. In that case the SVD returns a
/// deterministic minimum-norm representative of a non-unique `lambda`; the
/// later TPD feasibility contract makes physical decisions invariant to the
/// corresponding null-space freedom.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct ElementPotentialFit {
    potentials: Vec<f64>,
    max_abs_residual: f64,
    residual_tolerance: f64,
    rank: usize,
    singular_value_tolerance: f64,
    reference_species: Vec<usize>,
}

/// Immutable linear-algebra geometry of one active reference assemblage.
///
/// Chemical potentials change at every accepted temperature point, but the
/// elemental matrix does not while the active layout is unchanged. Retaining
/// this SVD avoids rebuilding the same range/null-space description for every
/// TPD probe. The geometry is scoped to one prepared phase-control runner;
/// it never crosses unrelated system layouts.
#[derive(Debug, Clone)]
pub(crate) struct ElementalReferenceGeometry {
    reference_matrix: DMatrix<f64>,
    left_vectors: DMatrix<f64>,
    right_vectors_transpose: DMatrix<f64>,
    singular_values: Vec<f64>,
    singular_value_tolerance: f64,
    rank: usize,
    nullspace_basis: DMatrix<f64>,
}

impl ElementalReferenceGeometry {
    /// Builds the reusable SVD/range/null-space description for active rows.
    pub(crate) fn new(reference_matrix: DMatrix<f64>) -> Result<Self, ReactionExtentError> {
        if reference_matrix.nrows() == 0 || reference_matrix.ncols() == 0 {
            return Err(ReactionExtentError::DimensionMismatch(format!(
                "elemental reference geometry requires a non-empty matrix, got {}x{}",
                reference_matrix.nrows(),
                reference_matrix.ncols(),
            )));
        }
        if reference_matrix.iter().any(|value| !value.is_finite()) {
            return Err(ReactionExtentError::InvalidCandidate {
                field: "elemental_reference_geometry",
                message: "elemental reference matrix contains non-finite values".to_string(),
            });
        }

        let svd = SVD::new(reference_matrix.clone(), true, true);
        let singular_values = svd.singular_values.iter().copied().collect::<Vec<_>>();
        let singular_scale = singular_values.iter().copied().fold(0.0_f64, f64::max);
        let singular_value_tolerance = (singular_scale * 1e-12).max(1e-12);
        let rank = singular_values
            .iter()
            .filter(|&&value| value > singular_value_tolerance)
            .count();
        let left_vectors = svd.u.ok_or_else(|| ReactionExtentError::InvalidProblem {
            field: "elemental_reference_geometry",
            message: "SVD did not retain left singular vectors".to_string(),
        })?;
        let right_vectors_transpose =
            svd.v_t.ok_or_else(|| ReactionExtentError::InvalidProblem {
                field: "elemental_reference_geometry",
                message: "SVD did not retain right singular vectors".to_string(),
            })?;
        let element_count = reference_matrix.ncols();
        // `nalgebra` keeps a thin V^T for rectangular SVDs. The full element
        // null space is therefore recovered from A^T A only once here, not
        // guessed from absent SVD rows on every TPD evaluation.
        let gram = reference_matrix.transpose() * &reference_matrix;
        let eigen = SymmetricEigen::new(gram);
        let eigen_tolerance = singular_value_tolerance * singular_value_tolerance;
        let null_columns = (0..element_count)
            .filter(|&index| eigen.eigenvalues[index].abs() <= eigen_tolerance)
            .collect::<Vec<_>>();
        let nullspace_basis = DMatrix::from_fn(element_count, null_columns.len(), |row, column| {
            eigen.eigenvectors[(row, null_columns[column])]
        });

        Ok(Self {
            reference_matrix,
            left_vectors,
            right_vectors_transpose,
            singular_values,
            singular_value_tolerance,
            rank,
            nullspace_basis,
        })
    }

    pub(crate) fn fit_potentials(
        &self,
        reference_mu: &DVector<f64>,
        reference_species: &[usize],
    ) -> Result<ElementPotentialFit, ReactionExtentError> {
        if reference_mu.len() != self.reference_matrix.nrows()
            || reference_species.len() != self.reference_matrix.nrows()
            || reference_mu.iter().any(|value| !value.is_finite())
        {
            return Err(ReactionExtentError::DimensionMismatch(format!(
                "element-potential geometry has {} reference rows, {} chemical potentials, and {} reference species",
                self.reference_matrix.nrows(),
                reference_mu.len(),
                reference_species.len(),
            )));
        }
        let mut lambda = DVector::zeros(self.reference_matrix.ncols());
        for index in 0..self.singular_values.len() {
            let singular = self.singular_values[index];
            if singular <= self.singular_value_tolerance {
                continue;
            }
            let coefficient = self.left_vectors.column(index).dot(reference_mu) / singular;
            for element in 0..lambda.len() {
                lambda[element] += self.right_vectors_transpose[(index, element)] * coefficient;
            }
        }
        let residual = &self.reference_matrix * &lambda - reference_mu;
        let max_abs_residual = residual
            .iter()
            .fold(0.0_f64, |maximum, value| maximum.max(value.abs()));
        let chemical_potential_scale = reference_mu
            .iter()
            .copied()
            .fold(1.0_f64, |scale, value| scale.max(value.abs()));
        let residual_tolerance = (chemical_potential_scale * 1e-8).max(1e-6);
        if max_abs_residual > residual_tolerance {
            return Err(ReactionExtentError::InvalidCandidate {
                field: "element_potentials",
                message: format!(
                    "active chemical potentials cannot be represented by elemental potentials: residual {max_abs_residual} exceeds {residual_tolerance}"
                ),
            });
        }
        Ok(ElementPotentialFit {
            potentials: lambda.iter().copied().collect(),
            max_abs_residual,
            residual_tolerance,
            rank: self.rank,
            singular_value_tolerance: self.singular_value_tolerance,
            reference_species: reference_species.to_vec(),
        })
    }

    fn candidate_constraints(
        &self,
        candidate_element_composition: &DMatrix<f64>,
    ) -> Result<DMatrix<f64>, ReactionExtentError> {
        if candidate_element_composition.ncols() != self.reference_matrix.ncols()
            || candidate_element_composition
                .iter()
                .any(|value| !value.is_finite())
        {
            return Err(ReactionExtentError::DimensionMismatch(format!(
                "elemental reference geometry has {} element columns but candidate matrix is {}x{}",
                self.reference_matrix.ncols(),
                candidate_element_composition.nrows(),
                candidate_element_composition.ncols(),
            )));
        }
        let raw_constraints = candidate_element_composition * &self.nullspace_basis;
        let rows = raw_constraints
            .column_iter()
            .map(|column| column.iter().copied().collect::<Vec<_>>())
            .collect::<Vec<_>>();
        independent_row_basis(
            &rows,
            candidate_element_composition.nrows(),
            self.singular_value_tolerance,
        )
    }

    fn assess_feasibility(
        &self,
        candidate_element_composition: &DMatrix<f64>,
        composition: &[f64],
    ) -> Result<ElementalDirectionFeasibility, ReactionExtentError> {
        if candidate_element_composition.ncols() != self.reference_matrix.ncols()
            || candidate_element_composition.nrows() != composition.len()
        {
            return Err(ReactionExtentError::DimensionMismatch(format!(
                "elemental-direction geometry has active {}x{}, candidate {}x{}, and {} composition values",
                self.reference_matrix.nrows(),
                self.reference_matrix.ncols(),
                candidate_element_composition.nrows(),
                candidate_element_composition.ncols(),
                composition.len(),
            )));
        }
        validate_simplex(composition)?;
        let composition = DVector::from_column_slice(composition);
        let candidate_element_totals = candidate_element_composition.transpose() * composition;
        let mut coefficients = DVector::zeros(self.reference_matrix.nrows());
        for index in 0..self.singular_values.len() {
            let singular = self.singular_values[index];
            if singular <= self.singular_value_tolerance {
                continue;
            }
            let coefficient = (0..candidate_element_totals.len())
                .map(|element| {
                    self.right_vectors_transpose[(index, element)]
                        * candidate_element_totals[element]
                })
                .sum::<f64>()
                / singular;
            for species in 0..coefficients.len() {
                coefficients[species] += self.left_vectors[(species, index)] * coefficient;
            }
        }
        let residual = self.reference_matrix.transpose() * coefficients - &candidate_element_totals;
        let max_abs_residual = residual
            .iter()
            .fold(0.0_f64, |maximum, value| maximum.max(value.abs()));
        let elemental_scale = candidate_element_totals
            .iter()
            .copied()
            .fold(1.0_f64, |scale, value| scale.max(value.abs()));
        let residual_tolerance = (elemental_scale * 1e-10).max(1e-12);
        Ok(ElementalDirectionFeasibility {
            candidate_element_totals: candidate_element_totals.iter().copied().collect(),
            max_abs_residual,
            residual_tolerance,
        })
    }
}

/// Deterministic cache statistics for the pure TPD geometry layer.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub(crate) struct PhaseStabilityGeometryCacheStats {
    pub(crate) entries: usize,
    pub(crate) builds: usize,
    pub(crate) reuses: usize,
}

/// Runner-scoped cache of active elemental reference geometries.
#[derive(Debug, Default)]
pub(crate) struct PhaseStabilityGeometryCache {
    entries: HashMap<Vec<usize>, ElementalReferenceGeometry>,
    builds: usize,
    reuses: usize,
}

impl PhaseStabilityGeometryCache {
    pub(crate) fn geometry_for(
        &mut self,
        state: &CanonicalPhaseState,
        reference_species: &[usize],
    ) -> Result<&ElementalReferenceGeometry, ReactionExtentError> {
        let species_count = state.chemical_potentials().len();
        if reference_species.is_empty() {
            return Err(ReactionExtentError::ValidationNotApplicable {
                path: "element_potentials",
                message: "element-potential reconstruction requires active reference species"
                    .to_string(),
            });
        }
        let mut seen = vec![false; species_count];
        for &species in reference_species {
            if species >= species_count || seen[species] {
                return Err(ReactionExtentError::InvalidProblem {
                    field: "element_potential_reference_species",
                    message: format!(
                        "reference species index {species} is out of bounds or duplicated"
                    ),
                });
            }
            seen[species] = true;
        }
        let key = reference_species.to_vec();
        if self.entries.contains_key(&key) {
            self.reuses += 1;
        } else {
            let element_count = state.element_composition().ncols();
            let matrix = DMatrix::from_fn(reference_species.len(), element_count, |row, col| {
                state.element_composition()[(reference_species[row], col)]
            });
            self.entries
                .insert(key.clone(), ElementalReferenceGeometry::new(matrix)?);
            self.builds += 1;
        }
        // The entry was inserted or verified above; no fallback geometry is
        // constructed on a cache miss after an error.
        Ok(self
            .entries
            .get(&key)
            .expect("elemental geometry inserted above"))
    }

    pub(crate) fn statistics(&self) -> PhaseStabilityGeometryCacheStats {
        PhaseStabilityGeometryCacheStats {
            entries: self.entries.len(),
            builds: self.builds,
            reuses: self.reuses,
        }
    }
}

#[allow(dead_code)] // Typed reports consume the remaining diagnostics during P7.5.
impl ElementPotentialFit {
    pub(crate) fn potentials(&self) -> &[f64] {
        &self.potentials
    }

    pub(crate) fn max_abs_residual(&self) -> f64 {
        self.max_abs_residual
    }

    pub(crate) fn residual_tolerance(&self) -> f64 {
        self.residual_tolerance
    }

    pub(crate) fn rank(&self) -> usize {
        self.rank
    }

    pub(crate) fn reference_species(&self) -> &[usize] {
        &self.reference_species
    }
}

/// Feasibility evidence for one candidate phase composition.
///
/// The candidate elemental direction is admissible only when it lies in the
/// range of the active elemental map. This remains meaningful when the active
/// map is rank deficient: the fit residual, rather than rank alone, decides
/// feasibility.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct ElementalDirectionFeasibility {
    candidate_element_totals: Vec<f64>,
    max_abs_residual: f64,
    residual_tolerance: f64,
}

#[allow(dead_code)] // Candidate composition is retained for the typed report migration.
impl ElementalDirectionFeasibility {
    pub(crate) fn candidate_element_totals(&self) -> &[f64] {
        &self.candidate_element_totals
    }

    pub(crate) fn max_abs_residual(&self) -> f64 {
        self.max_abs_residual
    }

    pub(crate) fn residual_tolerance(&self) -> f64 {
        self.residual_tolerance
    }

    pub(crate) fn is_feasible(&self) -> bool {
        self.max_abs_residual <= self.residual_tolerance
    }
}

/// Evaluates whether a simplex composition of a candidate phase preserves an
/// elemental direction accessible to the active assemblage.
///
/// `active_element_composition` and `candidate_element_composition` are both
/// species-by-element matrices. The returned residual is from
/// `A_active^T * y ~= A_candidate^T * x`.
#[allow(dead_code)] // Retained as an independent one-shot feasibility oracle for tests and adapters.
pub(crate) fn assess_elemental_direction_feasibility(
    active_element_composition: &DMatrix<f64>,
    candidate_element_composition: &DMatrix<f64>,
    composition: &[f64],
) -> Result<ElementalDirectionFeasibility, ReactionExtentError> {
    if active_element_composition.nrows() == 0
        || active_element_composition.ncols() == 0
        || active_element_composition.ncols() != candidate_element_composition.ncols()
        || candidate_element_composition.nrows() != composition.len()
    {
        return Err(ReactionExtentError::DimensionMismatch(format!(
            "elemental-direction feasibility has active {}x{}, candidate {}x{}, and {} composition values",
            active_element_composition.nrows(),
            active_element_composition.ncols(),
            candidate_element_composition.nrows(),
            candidate_element_composition.ncols(),
            composition.len(),
        )));
    }
    if active_element_composition
        .iter()
        .any(|value| !value.is_finite())
        || candidate_element_composition
            .iter()
            .any(|value| !value.is_finite())
    {
        return Err(ReactionExtentError::InvalidCandidate {
            field: "elemental_direction",
            message: "elemental composition contains non-finite values".to_string(),
        });
    }
    validate_simplex(composition)?;

    let geometry = ElementalReferenceGeometry::new(active_element_composition.clone())?;
    geometry.assess_feasibility(candidate_element_composition, composition)
}

/// Minimum of an ideal candidate-phase tangent-plane-distance objective.
///
/// `component_reduced_potentials` are the quantities
/// `g_i^0 + R*T*phase_offset - a_i*lambda`. Keeping them in the result makes
/// the analytical minimum auditable and gives future constrained minimizers a
/// common input representation.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct IdealTpdMinimum {
    minimum_tpd: f64,
    incipient_composition: Vec<f64>,
    component_reduced_potentials: Vec<f64>,
}

/// Evidence emitted by the constrained ideal TPD minimizer.
///
/// The candidate composition is constrained to remain elementally reachable
/// from the active assemblage. The stored residuals are independent checks of
/// that contract and of the first-order optimum; neither is inferred merely
/// from a nonlinear iteration status.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct ConstrainedIdealTpdMinimum {
    minimum: IdealTpdMinimum,
    feasibility: ElementalDirectionFeasibility,
    independent_constraint_count: usize,
    active_component_count: usize,
    iterations: usize,
    max_abs_constraint_residual: f64,
    constraint_tolerance: f64,
    max_abs_kkt_residual: f64,
}

#[allow(dead_code)] // Detailed minimizer evidence is intentionally crate-private for now.
impl ConstrainedIdealTpdMinimum {
    pub(crate) fn minimum(&self) -> &IdealTpdMinimum {
        &self.minimum
    }

    pub(crate) fn feasibility(&self) -> &ElementalDirectionFeasibility {
        &self.feasibility
    }

    pub(crate) fn independent_constraint_count(&self) -> usize {
        self.independent_constraint_count
    }

    pub(crate) fn active_component_count(&self) -> usize {
        self.active_component_count
    }

    pub(crate) fn iterations(&self) -> usize {
        self.iterations
    }

    pub(crate) fn max_abs_constraint_residual(&self) -> f64 {
        self.max_abs_constraint_residual
    }

    pub(crate) fn constraint_tolerance(&self) -> f64 {
        self.constraint_tolerance
    }

    pub(crate) fn max_abs_kkt_residual(&self) -> f64 {
        self.max_abs_kkt_residual
    }
}

#[allow(dead_code)] // Full result accessors are covered by unit tests pending P7.5 reports.
impl IdealTpdMinimum {
    pub(crate) fn minimum_tpd(&self) -> f64 {
        self.minimum_tpd
    }

    pub(crate) fn incipient_composition(&self) -> &[f64] {
        &self.incipient_composition
    }

    pub(crate) fn component_reduced_potentials(&self) -> &[f64] {
        &self.component_reduced_potentials
    }
}

/// Fully specified ideal TPD objective for one declared candidate phase.
///
/// The same objective serves both the full-rank closed form and the
/// rank-deficient feasible-simplex minimizer. Its numerical implementation is
/// deliberately independent from phase-control mutation and accepts no trace
/// mole amount as part of the thermodynamic definition.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct IdealTpdProblem {
    activity_model: PhaseActivityModel,
    standard_gibbs: Vec<f64>,
    candidate_element_composition: DMatrix<f64>,
    element_potentials: Vec<f64>,
    temperature: f64,
    pressure: f64,
    reference_pressure: f64,
}

#[allow(dead_code)] // The closed-form helper remains an audited independent regression path.
impl IdealTpdProblem {
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn new(
        activity_model: PhaseActivityModel,
        standard_gibbs: Vec<f64>,
        candidate_element_composition: DMatrix<f64>,
        element_potentials: Vec<f64>,
        temperature: f64,
        pressure: f64,
        reference_pressure: f64,
    ) -> Result<Self, ReactionExtentError> {
        if standard_gibbs.is_empty()
            || candidate_element_composition.nrows() != standard_gibbs.len()
            || candidate_element_composition.ncols() != element_potentials.len()
        {
            return Err(ReactionExtentError::DimensionMismatch(format!(
                "ideal TPD has {} standard Gibbs values, a {}x{} element matrix, and {} elemental potentials",
                standard_gibbs.len(),
                candidate_element_composition.nrows(),
                candidate_element_composition.ncols(),
                element_potentials.len(),
            )));
        }
        if standard_gibbs.iter().any(|value| !value.is_finite())
            || candidate_element_composition
                .iter()
                .any(|value| !value.is_finite())
            || element_potentials.iter().any(|value| !value.is_finite())
        {
            return Err(ReactionExtentError::InvalidCandidate {
                field: "ideal_tpd",
                message: "ideal TPD inputs must be finite".to_string(),
            });
        }
        for (parameter, value) in [
            ("temperature", temperature),
            ("pressure", pressure),
            ("reference_pressure", reference_pressure),
        ] {
            if !value.is_finite() || value <= 0.0 {
                return Err(ReactionExtentError::InvalidConditions { parameter, value });
            }
        }
        activity_model.log_phase_offset(pressure, reference_pressure)?;

        Ok(Self {
            activity_model,
            standard_gibbs,
            candidate_element_composition,
            element_potentials,
            temperature,
            pressure,
            reference_pressure,
        })
    }

    /// Solves the unconstrained ideal simplex problem by stable log-sum-exp.
    ///
    /// For a one-component phase this is exactly the historical pure-phase
    /// expression. No trace amount and no log-mole floor enters this result.
    pub(crate) fn unconstrained_minimum(&self) -> Result<IdealTpdMinimum, ReactionExtentError> {
        let reduced = self.reduced_potentials()?;
        let composition = self.softmax_composition(
            &reduced,
            &DVector::zeros(0),
            &DMatrix::zeros(0, reduced.len()),
        )?;
        self.minimum_for_composition(reduced, composition)
    }

    /// Minimizes the ideal TPD over the elementally feasible candidate simplex.
    ///
    /// For a rank-deficient active elemental map, the left null space defines
    /// the exact linear constraints `B*x = 0`. The ideal TPD has a convex
    /// entropy form, so its dual variables can be solved with a deterministic
    /// Newton method whose Jacobian is the analytical softmax covariance. No
    /// penalty term or post-hoc normalization changes the constrained problem.
    ///
    /// A deterministic phase-I linear program identifies the relative interior
    /// support before the softmax solve. Boundary components remain exact
    /// mathematical zeroes in the result; a floor appears only when an
    /// activated phase is later written into log-mole coordinates.
    pub(crate) fn constrained_minimum(
        &self,
        active_element_composition: &DMatrix<f64>,
    ) -> Result<ConstrainedIdealTpdMinimum, ReactionExtentError> {
        if active_element_composition.nrows() == 0
            || active_element_composition.ncols() != self.candidate_element_composition.ncols()
        {
            return Err(ReactionExtentError::DimensionMismatch(format!(
                "constrained ideal TPD has active {}x{} and candidate {}x{} element matrices",
                active_element_composition.nrows(),
                active_element_composition.ncols(),
                self.candidate_element_composition.nrows(),
                self.candidate_element_composition.ncols(),
            )));
        }
        if active_element_composition
            .iter()
            .any(|value| !value.is_finite())
        {
            return Err(ReactionExtentError::InvalidCandidate {
                field: "constrained_ideal_tpd",
                message: "active element composition contains non-finite values".to_string(),
            });
        }

        let geometry = ElementalReferenceGeometry::new(active_element_composition.clone())?;
        self.constrained_minimum_with_geometry(&geometry)
    }

    /// Evaluates the constrained ideal minimum using a precomputed active
    /// elemental geometry. The caller may retain that geometry across outer
    /// iterations and temperature points, but the objective itself remains a
    /// pure function of the current thermodynamic state.
    pub(crate) fn constrained_minimum_with_geometry(
        &self,
        geometry: &ElementalReferenceGeometry,
    ) -> Result<ConstrainedIdealTpdMinimum, ReactionExtentError> {
        if geometry.reference_matrix.ncols() != self.candidate_element_composition.ncols() {
            return Err(ReactionExtentError::DimensionMismatch(format!(
                "constrained ideal TPD geometry has {} elements but candidate matrix has {}",
                geometry.reference_matrix.ncols(),
                self.candidate_element_composition.ncols(),
            )));
        }
        let reduced = self.reduced_potentials()?;
        let constraints = geometry.candidate_constraints(&self.candidate_element_composition)?;
        let constraint_scale = constraints
            .iter()
            .copied()
            .fold(1.0_f64, |scale, value| scale.max(value.abs()));
        let constraint_tolerance = (constraint_scale * 1e-10).max(1e-12);

        let support = feasible_component_support(&constraints)?;
        let reduced_support = support
            .iter()
            .map(|&component| reduced[component])
            .collect::<Vec<_>>();
        let support_constraints =
            DMatrix::from_fn(constraints.nrows(), support.len(), |row, column| {
                constraints[(row, support[column])]
            });
        let (support_composition, dual, iterations, _) = solve_ideal_tpd_dual(
            &reduced_support,
            R * self.temperature,
            &support_constraints,
            constraint_tolerance,
        )?;
        let mut composition = vec![0.0; reduced.len()];
        for (&component, &value) in support.iter().zip(&support_composition) {
            composition[component] = value;
        }
        let max_abs_constraint_residual = (&constraints * DVector::from_column_slice(&composition))
            .iter()
            .fold(0.0_f64, |maximum, value| maximum.max(value.abs()));
        validate_simplex(&composition)?;
        if max_abs_constraint_residual > constraint_tolerance {
            return Err(ReactionExtentError::InvalidCandidate {
                field: "constrained_ideal_tpd",
                message: format!(
                    "candidate elemental constraints retain residual {max_abs_constraint_residual:e}, exceeding {constraint_tolerance:e}"
                ),
            });
        }

        let feasibility =
            geometry.assess_feasibility(&self.candidate_element_composition, &composition)?;
        if !feasibility.is_feasible() {
            return Err(ReactionExtentError::InvalidCandidate {
                field: "phase_stability_elemental_direction",
                message: format!(
                    "constrained ideal TPD returned an elementally infeasible composition: residual {} exceeds {}",
                    feasibility.max_abs_residual(),
                    feasibility.residual_tolerance(),
                ),
            });
        }

        let max_abs_kkt_residual = ideal_tpd_kkt_residual(
            &reduced_support,
            R * self.temperature,
            &support_composition,
            &support_constraints,
            &dual,
        )?;
        let kkt_tolerance = ((R * self.temperature).abs() * 1e-9).max(1e-7);
        if max_abs_kkt_residual > kkt_tolerance {
            return Err(ReactionExtentError::InvalidCandidate {
                field: "constrained_ideal_tpd",
                message: format!(
                    "ideal TPD KKT residual {max_abs_kkt_residual:e} exceeds {kkt_tolerance:e}"
                ),
            });
        }

        Ok(ConstrainedIdealTpdMinimum {
            minimum: self.minimum_for_composition(reduced, composition)?,
            feasibility,
            independent_constraint_count: constraints.nrows(),
            active_component_count: support.len(),
            iterations,
            max_abs_constraint_residual,
            constraint_tolerance,
            max_abs_kkt_residual,
        })
    }

    fn reduced_potentials(&self) -> Result<Vec<f64>, ReactionExtentError> {
        let rt = R * self.temperature;
        let phase_offset = self
            .activity_model
            .log_phase_offset(self.pressure, self.reference_pressure)?;
        let reduced = self
            .standard_gibbs
            .iter()
            .enumerate()
            .map(|(component, &g0)| {
                let elemental_reference = (0..self.element_potentials.len())
                    .map(|element| {
                        self.candidate_element_composition[(component, element)]
                            * self.element_potentials[element]
                    })
                    .sum::<f64>();
                g0 + rt * phase_offset - elemental_reference
            })
            .collect::<Vec<_>>();
        if reduced.iter().any(|value| !value.is_finite()) {
            return Err(ReactionExtentError::InvalidCandidate {
                field: "ideal_tpd",
                message: "ideal TPD reduced potentials are non-finite".to_string(),
            });
        }
        Ok(reduced)
    }

    fn softmax_composition(
        &self,
        reduced: &[f64],
        dual: &DVector<f64>,
        constraints: &DMatrix<f64>,
    ) -> Result<Vec<f64>, ReactionExtentError> {
        softmax_composition(reduced, R * self.temperature, dual, constraints)
    }

    fn minimum_for_composition(
        &self,
        reduced: Vec<f64>,
        composition: Vec<f64>,
    ) -> Result<IdealTpdMinimum, ReactionExtentError> {
        validate_simplex(&composition)?;
        let rt = R * self.temperature;
        let minimum_tpd = composition
            .iter()
            .zip(&reduced)
            .map(|(&x, &reduced_potential)| {
                if x == 0.0 {
                    0.0
                } else {
                    x * (reduced_potential + rt * x.ln())
                }
            })
            .sum::<f64>();
        if !minimum_tpd.is_finite() {
            return Err(ReactionExtentError::InvalidCandidate {
                field: "ideal_tpd",
                message: "ideal TPD minimum is non-finite".to_string(),
            });
        }
        Ok(IdealTpdMinimum {
            minimum_tpd,
            incipient_composition: composition,
            component_reduced_potentials: reduced,
        })
    }
}

/// Retains the first linearly independent constraint rows in their declared
/// order. The right-hand side is zero for every elemental null-space row, so
/// dropping a dependent row cannot alter the feasible set.
fn independent_row_basis(
    raw_rows: &[Vec<f64>],
    column_count: usize,
    tolerance: f64,
) -> Result<DMatrix<f64>, ReactionExtentError> {
    let mut selected = Vec::<Vec<f64>>::new();
    let mut current_rank = 0;
    for row in raw_rows {
        if row.len() != column_count || row.iter().any(|value| !value.is_finite()) {
            return Err(ReactionExtentError::InvalidCandidate {
                field: "constrained_ideal_tpd",
                message: "elemental null-space constraint is malformed".to_string(),
            });
        }
        let mut candidate = selected.clone();
        candidate.push(row.clone());
        let matrix = DMatrix::from_fn(candidate.len(), column_count, |r, c| candidate[r][c]);
        let singular_scale = SVD::new(matrix, false, false)
            .singular_values
            .iter()
            .copied()
            .fold(0.0_f64, f64::max);
        let rank_tolerance = (singular_scale * 1e-12).max(tolerance);
        let rank = SVD::new(
            DMatrix::from_fn(candidate.len(), column_count, |r, c| candidate[r][c]),
            false,
            false,
        )
        .singular_values
        .iter()
        .filter(|&&value| value > rank_tolerance)
        .count();
        if rank > current_rank {
            selected.push(row.clone());
            current_rank = rank;
        }
    }
    Ok(DMatrix::from_fn(selected.len(), column_count, |r, c| {
        selected[r][c]
    }))
}

/// Finds the relative-interior support of `B*x = 0, x >= 0, sum(x) = 1`.
///
/// Each component is retained only if a deterministic linear program can make
/// it strictly positive. The union of those components is the relative
/// interior support of the feasible polytope, so the entropy TPD minimum on
/// that face has finite softmax coordinates without inventing a composition
/// floor.
fn feasible_component_support(
    constraints: &DMatrix<f64>,
) -> Result<Vec<usize>, ReactionExtentError> {
    let component_count = constraints.ncols();
    if component_count == 0 {
        return Err(ReactionExtentError::DimensionMismatch(
            "TPD feasibility requires at least one candidate component".to_string(),
        ));
    }
    if constraints.nrows() == 0 {
        return Ok((0..component_count).collect());
    }

    let mut coefficients = Vec::with_capacity(2 * constraints.nrows() + 2);
    let mut bounds = Vec::with_capacity(2 * constraints.nrows() + 2);
    for row in 0..constraints.nrows() {
        let positive = (0..component_count)
            .map(|column| constraints[(row, column)])
            .collect::<Vec<_>>();
        let negative = positive.iter().map(|value| -value).collect::<Vec<_>>();
        coefficients.push(positive);
        bounds.push(0.0);
        coefficients.push(negative);
        bounds.push(0.0);
    }
    coefficients.push(vec![1.0; component_count]);
    bounds.push(1.0);
    coefficients.push(vec![-1.0; component_count]);
    bounds.push(-1.0);

    let zero_objective = vec![0.0; component_count];
    let feasibility =
        DeterministicSimplex::new(&coefficients, &bounds, &zero_objective)?.solve()?;
    if matches!(feasibility, LinearProgramOutcome::Infeasible) {
        return Err(ReactionExtentError::InvalidCandidate {
            field: "phase_stability_elemental_direction",
            message: "candidate phase has no elementally feasible simplex composition".to_string(),
        });
    }
    if matches!(feasibility, LinearProgramOutcome::Unbounded) {
        return Err(ReactionExtentError::InvalidProblem {
            field: "constrained_ideal_tpd",
            message: "bounded simplex feasibility problem was unexpectedly unbounded".to_string(),
        });
    }

    let mut support = Vec::new();
    for component in 0..component_count {
        let mut objective = vec![0.0; component_count];
        objective[component] = 1.0;
        match DeterministicSimplex::new(&coefficients, &bounds, &objective)?.solve()? {
            LinearProgramOutcome::Optimal { value, .. } if value > 1e-12 => support.push(component),
            LinearProgramOutcome::Optimal { .. } => {}
            LinearProgramOutcome::Infeasible => {
                return Err(ReactionExtentError::InvalidProblem {
                    field: "constrained_ideal_tpd",
                    message: "a feasible simplex became infeasible while inspecting support"
                        .to_string(),
                });
            }
            LinearProgramOutcome::Unbounded => {
                return Err(ReactionExtentError::InvalidProblem {
                    field: "constrained_ideal_tpd",
                    message: "bounded simplex support problem was unexpectedly unbounded"
                        .to_string(),
                });
            }
        }
    }
    if support.is_empty() {
        return Err(ReactionExtentError::InvalidCandidate {
            field: "phase_stability_elemental_direction",
            message: "candidate phase has no positive component on its feasible simplex"
                .to_string(),
        });
    }
    Ok(support)
}

#[derive(Debug)]
enum LinearProgramOutcome {
    Optimal { value: f64 },
    Infeasible,
    Unbounded,
}

/// Small deterministic two-phase simplex solver used only for TPD support
/// discovery. It has no thermodynamic policy: callers still independently
/// validate the resulting composition and TPD KKT conditions.
struct DeterministicSimplex {
    constraints: usize,
    variables: usize,
    basis: Vec<isize>,
    nonbasis: Vec<isize>,
    tableau: Vec<Vec<f64>>,
}

impl DeterministicSimplex {
    fn new(
        coefficients: &[Vec<f64>],
        bounds: &[f64],
        objective: &[f64],
    ) -> Result<Self, ReactionExtentError> {
        let constraints = coefficients.len();
        let variables = objective.len();
        if constraints == 0
            || variables == 0
            || bounds.len() != constraints
            || coefficients
                .iter()
                .any(|row| row.len() != variables || row.iter().any(|value| !value.is_finite()))
            || bounds.iter().any(|value| !value.is_finite())
            || objective.iter().any(|value| !value.is_finite())
        {
            return Err(ReactionExtentError::InvalidProblem {
                field: "constrained_ideal_tpd",
                message: "linear feasibility program has invalid dimensions or non-finite values"
                    .to_string(),
            });
        }
        let mut tableau = vec![vec![0.0; variables + 2]; constraints + 2];
        for row in 0..constraints {
            for column in 0..variables {
                tableau[row][column] = coefficients[row][column];
            }
            tableau[row][variables] = -1.0;
            tableau[row][variables + 1] = bounds[row];
        }
        for column in 0..variables {
            tableau[constraints][column] = -objective[column];
        }
        tableau[constraints + 1][variables] = 1.0;
        Ok(Self {
            constraints,
            variables,
            basis: (0..constraints)
                .map(|row| (variables + row) as isize)
                .collect(),
            nonbasis: (0..=variables)
                .map(|column| {
                    if column == variables {
                        -1
                    } else {
                        column as isize
                    }
                })
                .collect(),
            tableau,
        })
    }

    fn solve(mut self) -> Result<LinearProgramOutcome, ReactionExtentError> {
        const EPS: f64 = 1e-10;
        let rhs = self.variables + 1;
        let artificial = self.variables;
        let mut pivot_row = 0;
        for row in 1..self.constraints {
            if self.tableau[row][rhs] < self.tableau[pivot_row][rhs] {
                pivot_row = row;
            }
        }
        if self.tableau[pivot_row][rhs] < -EPS {
            self.pivot(pivot_row, artificial)?;
            if !self.simplex(true)? || self.tableau[self.constraints + 1][rhs] < -EPS {
                return Ok(LinearProgramOutcome::Infeasible);
            }
            if self.tableau[self.constraints + 1][rhs].abs() > EPS {
                return Ok(LinearProgramOutcome::Infeasible);
            }
            if let Some(row) = self.basis.iter().position(|&value| value == -1) {
                let mut column = 0;
                for candidate in 1..=self.variables {
                    if self.tableau[row][candidate] < self.tableau[row][column] - EPS
                        || ((self.tableau[row][candidate] - self.tableau[row][column]).abs() <= EPS
                            && self.nonbasis[candidate] < self.nonbasis[column])
                    {
                        column = candidate;
                    }
                }
                self.pivot(row, column)?;
            }
        }
        if !self.simplex(false)? {
            return Ok(LinearProgramOutcome::Unbounded);
        }
        Ok(LinearProgramOutcome::Optimal {
            value: self.tableau[self.constraints][rhs],
        })
    }

    fn simplex(&mut self, phase_one: bool) -> Result<bool, ReactionExtentError> {
        const EPS: f64 = 1e-10;
        let objective_row = if phase_one {
            self.constraints + 1
        } else {
            self.constraints
        };
        let rhs = self.variables + 1;
        loop {
            let mut entering = None;
            for column in 0..=self.variables {
                if !phase_one && self.nonbasis[column] == -1 {
                    continue;
                }
                match entering {
                    None => entering = Some(column),
                    Some(best)
                        if self.tableau[objective_row][column]
                            < self.tableau[objective_row][best] - EPS
                            || ((self.tableau[objective_row][column]
                                - self.tableau[objective_row][best])
                                .abs()
                                <= EPS
                                && self.nonbasis[column] < self.nonbasis[best]) =>
                    {
                        entering = Some(column)
                    }
                    _ => {}
                }
            }
            let entering = entering.ok_or_else(|| ReactionExtentError::InvalidProblem {
                field: "constrained_ideal_tpd",
                message: "linear feasibility simplex has no entering column".to_string(),
            })?;
            if self.tableau[objective_row][entering] >= -EPS {
                return Ok(true);
            }
            let mut leaving = None;
            for row in 0..self.constraints {
                if self.tableau[row][entering] <= EPS {
                    continue;
                }
                match leaving {
                    None => leaving = Some(row),
                    Some(best) => {
                        let candidate_ratio = self.tableau[row][rhs] / self.tableau[row][entering];
                        let best_ratio = self.tableau[best][rhs] / self.tableau[best][entering];
                        if candidate_ratio < best_ratio - EPS
                            || ((candidate_ratio - best_ratio).abs() <= EPS
                                && self.basis[row] < self.basis[best])
                        {
                            leaving = Some(row);
                        }
                    }
                }
            }
            let Some(leaving) = leaving else {
                return Ok(false);
            };
            self.pivot(leaving, entering)?;
        }
    }

    fn pivot(&mut self, row: usize, column: usize) -> Result<(), ReactionExtentError> {
        let pivot = self.tableau[row][column];
        if !pivot.is_finite() || pivot.abs() <= 1e-14 {
            return Err(ReactionExtentError::InvalidProblem {
                field: "constrained_ideal_tpd",
                message: "linear feasibility simplex selected an invalid pivot".to_string(),
            });
        }
        let rows = self.constraints + 2;
        let columns = self.variables + 2;
        for other_row in 0..rows {
            if other_row == row {
                continue;
            }
            for other_column in 0..columns {
                if other_column == column {
                    continue;
                }
                self.tableau[other_row][other_column] -=
                    self.tableau[row][other_column] * self.tableau[other_row][column] / pivot;
            }
        }
        for other_column in 0..columns {
            if other_column != column {
                self.tableau[row][other_column] /= pivot;
            }
        }
        for other_row in 0..rows {
            if other_row != row {
                self.tableau[other_row][column] /= -pivot;
            }
        }
        self.tableau[row][column] = 1.0 / pivot;
        std::mem::swap(&mut self.basis[row], &mut self.nonbasis[column]);
        Ok(())
    }
}

/// Solves the ideal TPD dual `B*softmax(-(d + B^T*nu)/(R*T)) = 0`.
fn solve_ideal_tpd_dual(
    reduced: &[f64],
    rt: f64,
    constraints: &DMatrix<f64>,
    tolerance: f64,
) -> Result<(Vec<f64>, DVector<f64>, usize, f64), ReactionExtentError> {
    if constraints.nrows() == 0 {
        let composition = softmax_composition(reduced, rt, &DVector::zeros(0), constraints)?;
        return Ok((composition, DVector::zeros(0), 0, 0.0));
    }
    let mut dual = DVector::zeros(constraints.nrows());
    let maximum_iterations = 64;
    for iteration in 0..maximum_iterations {
        let composition = softmax_composition(reduced, rt, &dual, constraints)?;
        let composition_vector = DVector::from_column_slice(&composition);
        let residual = constraints * &composition_vector;
        let max_abs_residual = residual
            .iter()
            .fold(0.0_f64, |maximum, value| maximum.max(value.abs()));
        if max_abs_residual <= tolerance {
            return Ok((composition, dual, iteration, max_abs_residual));
        }

        let covariance = DMatrix::from_fn(composition.len(), composition.len(), |row, col| {
            if row == col {
                composition[row] * (1.0 - composition[col])
            } else {
                -composition[row] * composition[col]
            }
        });
        let jacobian = -(constraints * covariance * constraints.transpose()) / rt;
        let step = SVD::new(jacobian, true, true)
            .solve(&(-residual), tolerance)
            .map_err(|message| ReactionExtentError::InvalidProblem {
                field: "constrained_ideal_tpd",
                message: format!("could not solve ideal TPD dual Newton step: {message}"),
            })?;
        if step.iter().any(|value| !value.is_finite()) {
            return Err(ReactionExtentError::InvalidCandidate {
                field: "constrained_ideal_tpd",
                message: "ideal TPD dual Newton step is non-finite".to_string(),
            });
        }

        let mut accepted = false;
        let mut step_scale = 1.0;
        for _ in 0..24 {
            let candidate_dual = &dual + step_scale * &step;
            let candidate_composition =
                softmax_composition(reduced, rt, &candidate_dual, constraints)?;
            let candidate_residual =
                constraints * DVector::from_column_slice(&candidate_composition);
            let candidate_max = candidate_residual
                .iter()
                .fold(0.0_f64, |maximum, value| maximum.max(value.abs()));
            if candidate_max.is_finite() && candidate_max < max_abs_residual {
                dual = candidate_dual;
                accepted = true;
                break;
            }
            step_scale *= 0.5;
        }
        if !accepted {
            return Err(ReactionExtentError::InvalidCandidate {
                field: "constrained_ideal_tpd",
                message: format!(
                    "could not reduce elemental constraint residual from {max_abs_residual:e}; the feasible simplex may be empty or boundary-only"
                ),
            });
        }
    }
    Err(ReactionExtentError::InvalidCandidate {
        field: "constrained_ideal_tpd",
        message: format!(
            "dual Newton exceeded {maximum_iterations} iterations; the feasible simplex may be boundary-only"
        ),
    })
}

fn softmax_composition(
    reduced: &[f64],
    rt: f64,
    dual: &DVector<f64>,
    constraints: &DMatrix<f64>,
) -> Result<Vec<f64>, ReactionExtentError> {
    if !rt.is_finite()
        || rt <= 0.0
        || constraints.ncols() != reduced.len()
        || dual.len() != constraints.nrows()
    {
        return Err(ReactionExtentError::DimensionMismatch(format!(
            "ideal TPD softmax has {} reduced potentials, {}x{} constraints, and {} dual values",
            reduced.len(),
            constraints.nrows(),
            constraints.ncols(),
            dual.len(),
        )));
    }
    let shifted = constraints.transpose() * dual;
    let log_weights = reduced
        .iter()
        .enumerate()
        .map(|(index, &value)| -(value + shifted[index]) / rt)
        .collect::<Vec<_>>();
    let maximum = log_weights
        .iter()
        .copied()
        .fold(f64::NEG_INFINITY, f64::max);
    let sum = log_weights
        .iter()
        .map(|value| (value - maximum).exp())
        .sum::<f64>();
    if !sum.is_finite() || sum <= 0.0 {
        return Err(ReactionExtentError::InvalidCandidate {
            field: "constrained_ideal_tpd",
            message: "ideal TPD softmax normalization is invalid".to_string(),
        });
    }
    let log_normalizer = maximum + sum.ln();
    let composition = log_weights
        .iter()
        .map(|value| (value - log_normalizer).exp())
        .collect::<Vec<_>>();
    validate_simplex(&composition)?;
    Ok(composition)
}

fn ideal_tpd_kkt_residual(
    reduced: &[f64],
    rt: f64,
    composition: &[f64],
    constraints: &DMatrix<f64>,
    dual: &DVector<f64>,
) -> Result<f64, ReactionExtentError> {
    let dual_term = constraints.transpose() * dual;
    let stationarity = reduced
        .iter()
        .enumerate()
        .map(|(index, &value)| value + rt * (composition[index].ln() + 1.0) + dual_term[index])
        .collect::<Vec<_>>();
    let mean = stationarity.iter().sum::<f64>() / stationarity.len() as f64;
    let residual = stationarity
        .iter()
        .map(|value| (value - mean).abs())
        .fold(0.0_f64, f64::max);
    if !residual.is_finite() {
        return Err(ReactionExtentError::InvalidCandidate {
            field: "constrained_ideal_tpd",
            message: "ideal TPD KKT residual is non-finite".to_string(),
        });
    }
    Ok(residual)
}

fn validate_simplex(composition: &[f64]) -> Result<(), ReactionExtentError> {
    if composition.is_empty()
        || composition
            .iter()
            .any(|value| !value.is_finite() || *value < 0.0)
    {
        return Err(ReactionExtentError::InvalidCandidate {
            field: "candidate_composition",
            message: "candidate composition must be finite, non-negative, and non-empty"
                .to_string(),
        });
    }
    let sum = composition.iter().sum::<f64>();
    if (sum - 1.0).abs() > 1e-10 {
        return Err(ReactionExtentError::InvalidCandidate {
            field: "candidate_composition",
            message: format!("candidate composition must sum to one, got {sum}"),
        });
    }
    Ok(())
}

/// Fits elemental potentials for selected active species in a canonical state.
#[allow(dead_code)] // Retained as a pure one-shot reference path beside the runner cache.
pub(crate) fn reconstruct_element_potentials(
    state: &CanonicalPhaseState,
    reference_species: &[usize],
) -> Result<ElementPotentialFit, ReactionExtentError> {
    let species_count = state.chemical_potentials().len();
    if reference_species.is_empty() {
        return Err(ReactionExtentError::ValidationNotApplicable {
            path: "element_potentials",
            message: "element-potential reconstruction requires active reference species"
                .to_string(),
        });
    }
    let mut seen = vec![false; species_count];
    for &species in reference_species {
        if species >= species_count || seen[species] {
            return Err(ReactionExtentError::InvalidProblem {
                field: "element_potential_reference_species",
                message: format!(
                    "reference species index {species} is out of bounds or duplicated"
                ),
            });
        }
        seen[species] = true;
    }

    let element_count = state.element_composition().ncols();
    if element_count == 0 {
        return Err(ReactionExtentError::DimensionMismatch(
            "element-potential reconstruction requires at least one element column".to_string(),
        ));
    }
    let reference_matrix = DMatrix::from_fn(reference_species.len(), element_count, |row, col| {
        state.element_composition()[(reference_species[row], col)]
    });
    let reference_mu = DVector::from_iterator(
        reference_species.len(),
        reference_species
            .iter()
            .map(|&species| state.chemical_potentials()[species]),
    );
    ElementalReferenceGeometry::new(reference_matrix)?
        .fit_potentials(&reference_mu, reference_species)
}

/// Fits `A * lambda ~= mu` and verifies the residual independently of rank.
#[allow(dead_code)] // Retained as a pure one-shot reference path beside the runner cache.
pub(crate) fn fit_element_potentials(
    reference_matrix: &DMatrix<f64>,
    reference_mu: &DVector<f64>,
    reference_species: &[usize],
) -> Result<ElementPotentialFit, ReactionExtentError> {
    if reference_matrix.nrows() == 0
        || reference_matrix.ncols() == 0
        || reference_matrix.nrows() != reference_mu.len()
        || reference_matrix.nrows() != reference_species.len()
    {
        return Err(ReactionExtentError::DimensionMismatch(format!(
            "element-potential fit has {} rows, {} columns, {} chemical potentials, and {} reference species",
            reference_matrix.nrows(),
            reference_matrix.ncols(),
            reference_mu.len(),
            reference_species.len(),
        )));
    }
    if reference_matrix.iter().any(|value| !value.is_finite())
        || reference_mu.iter().any(|value| !value.is_finite())
    {
        return Err(ReactionExtentError::InvalidCandidate {
            field: "element_potentials",
            message: "element-potential equations contain non-finite values".to_string(),
        });
    }

    ElementalReferenceGeometry::new(reference_matrix.clone())?
        .fit_potentials(reference_mu, reference_species)
}

/// Presentation-safe evidence for the elemental-potential reconstruction.
#[derive(Debug, Clone, PartialEq)]
pub struct ElementPotentialReport {
    pub potentials: Vec<f64>,
    pub max_abs_residual: f64,
    pub rank: usize,
}

impl From<ElementPotentialFit> for ElementPotentialReport {
    fn from(fit: ElementPotentialFit) -> Self {
        Self {
            potentials: fit.potentials,
            max_abs_residual: fit.max_abs_residual,
            rank: fit.rank,
        }
    }
}

/// Thermodynamic conditions under which one TPD report was evaluated.
///
/// Values are copied into the evidence rather than inferred later from a
/// mutable solver object. This matters for `P,H`, where the accepted
/// temperature can differ from the initial temperature seed.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PhaseStabilityConditions {
    pub temperature: f64,
    pub pressure: f64,
    pub reference_pressure: f64,
}

/// Canonical component ordering used by one phase-stability report.
///
/// `phase_component_indices` are indices in the full resolved species layout,
/// so an incipient composition is never ambiguous even if the same molecular
/// formula appears in multiple declared phases.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PhaseStabilityLayout {
    pub system_species_count: usize,
    pub system_phase_count: usize,
    pub element_count: usize,
    pub phase_component_indices: Vec<usize>,
}

/// Auditable elemental-direction feasibility evidence for a TPD minimizer.
#[derive(Debug, Clone, PartialEq)]
pub struct ElementalFeasibilityReport {
    pub candidate_element_totals: Vec<f64>,
    pub max_abs_residual: f64,
    pub residual_tolerance: f64,
}

impl From<&ElementalDirectionFeasibility> for ElementalFeasibilityReport {
    fn from(feasibility: &ElementalDirectionFeasibility) -> Self {
        Self {
            candidate_element_totals: feasibility.candidate_element_totals.clone(),
            max_abs_residual: feasibility.max_abs_residual,
            residual_tolerance: feasibility.residual_tolerance,
        }
    }
}

/// Numerical evidence emitted by the constrained ideal TPD minimizer.
#[derive(Debug, Clone, PartialEq)]
pub struct TpdMinimizerReport {
    pub independent_constraint_count: usize,
    pub active_component_count: usize,
    pub iterations: usize,
    pub max_abs_constraint_residual: f64,
    pub constraint_tolerance: f64,
    pub max_abs_kkt_residual: f64,
}

impl From<&ConstrainedIdealTpdMinimum> for TpdMinimizerReport {
    fn from(minimum: &ConstrainedIdealTpdMinimum) -> Self {
        Self {
            independent_constraint_count: minimum.independent_constraint_count,
            active_component_count: minimum.active_component_count,
            iterations: minimum.iterations,
            max_abs_constraint_residual: minimum.max_abs_constraint_residual,
            constraint_tolerance: minimum.constraint_tolerance,
            max_abs_kkt_residual: minimum.max_abs_kkt_residual,
        }
    }
}

/// Why an individual declared phase was or was not evaluated by the TPD path.
///
/// This status deliberately describes workflow policy rather than inventing a
/// second physical phase-model hierarchy. The actual activity law remains the
/// [`PhaseActivityModel`] attached to the resolved phase.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PhaseStabilityStatus {
    /// The phase has a valid TPD minimum and associated evidence.
    Evaluated,
    /// The caller explicitly excluded this phase from outer-loop transitions.
    ExcludedByPolicy,
    /// The current domain represents all ideal-gas species as one fixed
    /// assemblage, so this separately declared gas phase is not a candidate.
    FixedGasAssemblage,
    /// The phase is already active but no independent reference assemblage is
    /// available for a meaningful TPD comparison.
    ActiveWithoutReferenceAssemblage,
}

/// Typed stability evidence for one declared phase.
///
/// `minimum_tpd` is expressed in J/mol and is populated only after a genuine
/// constrained TPD minimization. `None` therefore means "not evaluated"; it
/// never masquerades as a neutral numerical stability result.
#[derive(Debug, Clone, PartialEq)]
pub struct PhaseStabilityReport {
    pub phase: PhaseIndex,
    pub active: bool,
    pub status: PhaseStabilityStatus,
    pub conditions: PhaseStabilityConditions,
    pub layout: PhaseStabilityLayout,
    pub minimum_tpd: Option<f64>,
    /// TPD minimizer in declared phase-component order, normalized to one.
    /// It is absent exactly when [`status`](Self::status) is not
    /// [`PhaseStabilityStatus::Evaluated`].
    pub incipient_composition: Option<Vec<f64>>,
    pub element_potentials: Option<ElementPotentialReport>,
    pub elemental_feasibility: Option<ElementalFeasibilityReport>,
    pub minimizer: Option<TpdMinimizerReport>,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_activity::PhaseActivityModel;
    use std::rc::Rc;

    #[test]
    fn canonical_snapshot_uses_the_shared_ideal_activity_contract() {
        let phases = vec![Phase {
            kind: PhaseActivityModel::IdealGas,
            species: vec![0, 1],
        }];
        let gibbs: Vec<GibbsFn> = vec![Rc::new(|_| 10.0), Rc::new(|_| 20.0)];
        let state = CanonicalPhaseState::from_log_moles(
            &[0.25_f64.ln(), 0.75_f64.ln()],
            &gibbs,
            &phases,
            &[0, 0],
            &DMatrix::from_row_slice(2, 1, &[1.0, 1.0]),
            1000.0,
            101_325.0,
            101_325.0,
            &[true],
            &[true],
        )
        .unwrap();

        assert_eq!(state.phase_totals, vec![1.0]);
        assert!(
            (state.chemical_potentials()[0] - (10.0 + R * 1000.0 * 0.25_f64.ln())).abs() < 1e-9
        );
        assert!(
            (state.chemical_potentials()[1] - (20.0 + R * 1000.0 * 0.75_f64.ln())).abs() < 1e-9
        );
    }

    #[test]
    fn rank_deficient_element_potential_fit_is_accepted_when_residual_is_small() {
        let matrix = DMatrix::from_row_slice(2, 2, &[1.0, 0.0, 2.0, 0.0]);
        let fit =
            fit_element_potentials(&matrix, &DVector::from_vec(vec![5.0, 10.0]), &[3, 7]).unwrap();

        assert_eq!(fit.rank(), 1);
        assert_eq!(fit.potentials().len(), 2);
        assert!(fit.max_abs_residual() < 1e-12);
        assert!(fit.max_abs_residual() <= fit.residual_tolerance());
        assert!(fit.singular_value_tolerance > 0.0);
        assert_eq!(fit.reference_species(), &[3, 7]);
    }

    #[test]
    fn inconsistent_element_potential_fit_is_rejected_with_typed_error() {
        let matrix = DMatrix::from_row_slice(2, 1, &[1.0, 1.0]);
        let result = fit_element_potentials(&matrix, &DVector::from_vec(vec![5.0, 8.0]), &[0, 1]);

        assert!(matches!(
            result,
            Err(ReactionExtentError::InvalidCandidate {
                field: "element_potentials",
                ..
            })
        ));
    }

    #[test]
    fn geometry_cache_reuses_active_element_factorization_without_changing_fit() {
        let phases = vec![Phase {
            kind: PhaseActivityModel::IdealGas,
            species: vec![0],
        }];
        let gibbs: Vec<GibbsFn> = vec![Rc::new(|_| 42.0)];
        let state = CanonicalPhaseState::from_log_moles(
            &[0.0],
            &gibbs,
            &phases,
            &[0],
            &DMatrix::from_row_slice(1, 1, &[1.0]),
            1000.0,
            101_325.0,
            101_325.0,
            &[true],
            &[true],
        )
        .unwrap();
        let reference_mu = DVector::from_vec(vec![42.0]);
        let mut cache = PhaseStabilityGeometryCache::default();

        let first = cache
            .geometry_for(&state, &[0])
            .unwrap()
            .fit_potentials(&reference_mu, &[0])
            .unwrap();
        let second = cache
            .geometry_for(&state, &[0])
            .unwrap()
            .fit_potentials(&reference_mu, &[0])
            .unwrap();
        let statistics = cache.statistics();

        assert_eq!(first.potentials(), second.potentials());
        assert_eq!(statistics.entries, 1);
        assert_eq!(statistics.builds, 1);
        assert_eq!(statistics.reuses, 1);
    }

    #[test]
    fn ideal_tpd_pure_phase_matches_the_historical_chemical_potential_formula() {
        let problem = IdealTpdProblem::new(
            PhaseActivityModel::IdealSolution,
            vec![12.5],
            DMatrix::from_row_slice(1, 2, &[2.0, 1.0]),
            vec![3.0, -4.0],
            900.0,
            101_325.0,
            101_325.0,
        )
        .unwrap();

        let minimum = problem.unconstrained_minimum().unwrap();
        assert!((minimum.minimum_tpd() - 10.5).abs() < 1e-12);
        assert_eq!(minimum.incipient_composition(), &[1.0]);
        assert_eq!(minimum.component_reduced_potentials(), &[10.5]);
    }

    #[test]
    fn ideal_tpd_multicomponent_minimum_matches_closed_form_composition() {
        let temperature = 1000.0;
        let rt = R * temperature;
        let problem = IdealTpdProblem::new(
            PhaseActivityModel::IdealSolution,
            vec![0.0, rt * 2.0_f64.ln()],
            DMatrix::from_row_slice(2, 1, &[1.0, 1.0]),
            vec![0.0],
            temperature,
            101_325.0,
            101_325.0,
        )
        .unwrap();

        let minimum = problem.unconstrained_minimum().unwrap();
        assert!((minimum.minimum_tpd() + rt * 1.5_f64.ln()).abs() < 1e-9);
        assert!((minimum.incipient_composition()[0] - 2.0 / 3.0).abs() < 1e-12);
        assert!((minimum.incipient_composition()[1] - 1.0 / 3.0).abs() < 1e-12);
    }

    #[test]
    fn ideal_tpd_is_invariant_under_candidate_component_permutation() {
        let temperature = 1000.0;
        let rt = R * temperature;
        let original = IdealTpdProblem::new(
            PhaseActivityModel::IdealSolution,
            vec![0.0, rt * 2.0_f64.ln()],
            DMatrix::from_row_slice(2, 1, &[1.0, 1.0]),
            vec![0.0],
            temperature,
            101_325.0,
            101_325.0,
        )
        .unwrap()
        .unconstrained_minimum()
        .unwrap();
        let permuted = IdealTpdProblem::new(
            PhaseActivityModel::IdealSolution,
            vec![rt * 2.0_f64.ln(), 0.0],
            DMatrix::from_row_slice(2, 1, &[1.0, 1.0]),
            vec![0.0],
            temperature,
            101_325.0,
            101_325.0,
        )
        .unwrap()
        .unconstrained_minimum()
        .unwrap();

        assert!((original.minimum_tpd() - permuted.minimum_tpd()).abs() < 1e-12);
        assert!(
            (original.incipient_composition()[0] - permuted.incipient_composition()[1]).abs()
                < 1e-12
        );
        assert!(
            (original.incipient_composition()[1] - permuted.incipient_composition()[0]).abs()
                < 1e-12
        );
    }

    #[test]
    fn rank_deficient_active_assemblage_distinguishes_feasible_and_infeasible_directions() {
        let active = DMatrix::from_row_slice(1, 2, &[2.0, 0.0]);
        let candidate = DMatrix::from_row_slice(2, 2, &[2.0, 0.0, 0.0, 2.0]);

        let feasible =
            assess_elemental_direction_feasibility(&active, &candidate, &[1.0, 0.0]).unwrap();
        assert!(feasible.is_feasible());
        assert!(feasible.max_abs_residual() <= feasible.residual_tolerance());
        assert_eq!(feasible.candidate_element_totals(), &[2.0, 0.0]);

        let infeasible =
            assess_elemental_direction_feasibility(&active, &candidate, &[0.5, 0.5]).unwrap();
        assert!(!infeasible.is_feasible());
        assert!(infeasible.max_abs_residual() > infeasible.residual_tolerance());
    }

    #[test]
    fn constrained_ideal_tpd_solves_a_rank_deficient_interior_simplex() {
        let temperature = 900.0;
        let rt = R * temperature;
        // The active assemblage can exchange H and O only in a 1:1 ratio.
        // The candidate must therefore use its H-only and O-only components
        // equally; the admissible optimum is the interior point x=(1/2,1/2).
        let active = DMatrix::from_row_slice(1, 2, &[1.0, 1.0]);
        let problem = IdealTpdProblem::new(
            PhaseActivityModel::IdealSolution,
            vec![0.0, 0.0],
            DMatrix::from_row_slice(2, 2, &[2.0, 0.0, 0.0, 2.0]),
            vec![0.0, 0.0],
            temperature,
            101_325.0,
            101_325.0,
        )
        .unwrap();

        let result = problem.constrained_minimum(&active).unwrap();
        assert_eq!(result.independent_constraint_count(), 1);
        assert!(result.iterations() <= 64);
        assert!((result.minimum().minimum_tpd() + rt * 2.0_f64.ln()).abs() < 1e-8);
        assert!((result.minimum().incipient_composition()[0] - 0.5).abs() < 1e-10);
        assert!((result.minimum().incipient_composition()[1] - 0.5).abs() < 1e-10);
        assert!(result.feasibility().is_feasible());
        assert!(
            result.max_abs_constraint_residual() <= result.constraint_tolerance(),
            "constraint residual {} must not exceed {}",
            result.max_abs_constraint_residual(),
            result.constraint_tolerance(),
        );
        assert!(result.max_abs_kkt_residual() < rt * 1e-9);
    }

    #[test]
    fn constrained_ideal_tpd_is_invariant_to_rank_deficient_lambda_representatives() {
        let temperature = 750.0;
        let active = DMatrix::from_row_slice(1, 2, &[1.0, 1.0]);
        let candidate_elements = DMatrix::from_row_slice(2, 2, &[2.0, 0.0, 0.0, 2.0]);
        let minimum = |lambda| {
            IdealTpdProblem::new(
                PhaseActivityModel::IdealSolution,
                vec![0.0, 0.0],
                candidate_elements.clone(),
                lambda,
                temperature,
                101_325.0,
                101_325.0,
            )
            .unwrap()
            .constrained_minimum(&active)
            .unwrap()
        };

        let reference = minimum(vec![0.0, 0.0]);
        // [3, -3] is a null-space perturbation for the active H/O=1:1
        // assemblage. It changes individual reduced potentials but not a
        // feasible candidate direction or the constrained minimum.
        let shifted = minimum(vec![3.0, -3.0]);
        assert!((reference.minimum().minimum_tpd() - shifted.minimum().minimum_tpd()).abs() < 1e-8);
        assert_eq!(
            reference.minimum().incipient_composition(),
            shifted.minimum().incipient_composition()
        );
    }

    #[test]
    fn constrained_ideal_tpd_solves_a_boundary_face_without_putting_a_floor_in_the_objective() {
        let active = DMatrix::from_row_slice(1, 2, &[2.0, 0.0]);
        let problem = IdealTpdProblem::new(
            PhaseActivityModel::IdealSolution,
            vec![0.0, 0.0],
            DMatrix::from_row_slice(2, 2, &[2.0, 0.0, 0.0, 2.0]),
            vec![0.0, 0.0],
            1000.0,
            101_325.0,
            101_325.0,
        )
        .unwrap();

        let result = problem.constrained_minimum(&active).unwrap();
        assert_eq!(result.active_component_count(), 1);
        assert_eq!(result.minimum().incipient_composition(), &[1.0, 0.0]);
        assert!(result.feasibility().is_feasible());
        assert!(result.max_abs_constraint_residual() <= result.constraint_tolerance());
    }
}
