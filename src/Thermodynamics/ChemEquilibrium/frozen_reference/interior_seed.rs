//! Test-only construction of element-conserving interior composition seeds.
//!
//! This module is deliberately weaker than an equilibrium calculation.  It
//! knows only the active species, their element matrix, and the conserved
//! inventory.  In particular it does not inspect Gibbs energies, equilibrium
//! constants, or a previously accepted solution.  Its purpose is to tell
//! whether a fresh nonlinear seed can be moved away from a boundary while
//! remaining on the exact elemental-balance manifold.

use nalgebra::{DMatrix, DVector};

use crate::Thermodynamics::ChemEquilibrium::equilibrium_nonlinear::ReactionExtentError;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_problem::LogMolesInitialGuess;

/// Result category for a requested element-conserving interior seed.
///
/// The categories are diagnostic evidence, not production solve outcomes.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum InteriorSeedAvailability {
    /// Every active species satisfies the requested positive interior floor.
    FullInteriorAvailable,
    /// A positive seed exists only after lowering the requested floor.
    PartialInteriorOnly,
}

/// Structural composition target used before the affine projection.
///
/// Both choices derive only from the current feed and active species. Neither
/// introduces thermochemistry or a hidden equilibrium calculation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum InteriorSeedTarget {
    /// Treat all active species symmetrically.
    Uniform,
    /// Stay as close as possible to the caller's feed composition.
    InputAnchored,
}

/// Input-only settings for [`build_element_feasible_interior_seed`].
#[derive(Debug, Clone, Copy)]
pub(crate) struct InteriorSeedSettings {
    /// Answer-independent target to project onto the elemental manifold.
    pub(crate) target: InteriorSeedTarget,
    /// Requested active-species lower bound as a fraction of active inventory.
    pub(crate) requested_minimum_fraction: f64,
    /// Smallest fraction tried before declaring the affine problem unavailable.
    pub(crate) minimum_fraction: f64,
    /// Relative linear-algebra tolerance in normalized mole coordinates.
    pub(crate) tolerance: f64,
    /// Log-mole value assigned to physically inactive species.
    pub(crate) trace_floor: f64,
}

impl Default for InteriorSeedSettings {
    fn default() -> Self {
        Self {
            target: InteriorSeedTarget::Uniform,
            requested_minimum_fraction: 1.0e-6,
            minimum_fraction: 1.0e-12,
            tolerance: 1.0e-10,
            trace_floor: 1.0e-30,
        }
    }
}

/// A full-layout diagnostic seed with exact physical elemental inventory.
#[derive(Debug, Clone)]
pub(crate) struct ElementFeasibleInteriorSeed {
    /// Physical mole vector. Inactive phases are exactly zero here.
    pub(crate) physical_moles: Vec<f64>,
    /// Positive full-layout numerical seed; inactive phases use `trace_floor`.
    pub(crate) log_moles: LogMolesInitialGuess,
    /// Sum of caller-supplied physical moles in the active phase set.
    pub(crate) total_inventory_scale: f64,
    /// Fraction originally requested by the caller.
    pub(crate) requested_minimum_fraction: f64,
    /// Fraction actually achieved by the bounded affine projection.
    pub(crate) achieved_minimum_fraction: f64,
    /// Maximum absolute physical element-balance error after reconstruction.
    pub(crate) max_element_balance_error: f64,
    /// Number of species allowed to receive physical positive seed amounts.
    pub(crate) active_species_count: usize,
    /// Whether the requested floor was feasible without relaxation.
    pub(crate) availability: InteriorSeedAvailability,
}

/// Builds an answer-independent positive seed over the current active species.
///
/// It minimizes distance to a uniform, dimensionless composition under
/// `A_active^T n = b` and lower bounds `n_i >= epsilon`. The small active-set
/// loop is a bounded affine projection, not an LP implementation and not a
/// Gibbs minimizer. If the requested lower bound is infeasible, it is reduced
/// decade by decade down to `minimum_fraction`.
#[allow(clippy::too_many_arguments)]
pub(crate) fn build_element_feasible_interior_seed(
    element_composition: &DMatrix<f64>,
    element_totals: &[f64],
    species_phase: &[usize],
    active_phase_mask: &[bool],
    input_moles: &[f64],
    settings: InteriorSeedSettings,
) -> Result<ElementFeasibleInteriorSeed, ReactionExtentError> {
    validate_inputs(
        element_composition,
        element_totals,
        species_phase,
        active_phase_mask,
        input_moles,
        settings,
    )?;
    let active = species_phase
        .iter()
        .enumerate()
        .filter_map(|(index, &phase)| active_phase_mask[phase].then_some(index))
        .collect::<Vec<_>>();
    let scale = active.iter().map(|&index| input_moles[index]).sum::<f64>();
    if !scale.is_finite() || scale <= 0.0 {
        return Err(ReactionExtentError::InvalidProblem {
            field: "interior_seed_input_moles",
            message: "active physical inventory must have a finite positive total".to_owned(),
        });
    }

    let target = DVector::from_iterator(
        element_totals.len(),
        element_totals.iter().map(|&total| total / scale),
    );
    let active_matrix = DMatrix::from_fn(active.len(), element_totals.len(), |row, column| {
        element_composition[(active[row], column)]
    });
    let mut fraction = settings.requested_minimum_fraction;
    loop {
        let reference =
            reference_composition(settings.target, &active, input_moles, scale, fraction);
        if let Some(normalized) = bounded_affine_projection(
            &active_matrix,
            &target,
            &reference,
            fraction,
            settings.tolerance,
        )? {
            let mut physical_moles = vec![0.0; input_moles.len()];
            for (local, &global) in active.iter().enumerate() {
                physical_moles[global] = normalized[local] * scale;
            }
            let log_moles =
                LogMolesInitialGuess::from_moles(&physical_moles, settings.trace_floor)?;
            let max_element_balance_error =
                max_element_balance_error(element_composition, &physical_moles, element_totals);
            let achieved_minimum_fraction =
                normalized.iter().copied().fold(f64::INFINITY, f64::min);
            let availability = if (fraction - settings.requested_minimum_fraction).abs()
                <= settings.tolerance * settings.requested_minimum_fraction.max(1.0)
            {
                InteriorSeedAvailability::FullInteriorAvailable
            } else {
                InteriorSeedAvailability::PartialInteriorOnly
            };
            return Ok(ElementFeasibleInteriorSeed {
                physical_moles,
                log_moles,
                total_inventory_scale: scale,
                requested_minimum_fraction: settings.requested_minimum_fraction,
                achieved_minimum_fraction,
                max_element_balance_error,
                active_species_count: active.len(),
                availability,
            });
        }
        if fraction <= settings.minimum_fraction {
            break;
        }
        fraction = (fraction * 0.1).max(settings.minimum_fraction);
    }
    Err(ReactionExtentError::InvalidProblem {
        field: "element_feasible_interior_seed",
        message: format!(
            "no positive interior seed is available for the active species down to fraction {:e}",
            settings.minimum_fraction
        ),
    })
}

/// Forms a strictly positive, normalized target before the affine projection.
///
/// A raw feed can lie on a simplex face: most candidate species commonly have
/// zero input.  Passing that boundary point to a bounded projection makes the
/// set of saturated lower bounds depend on roundoff.  Injecting the requested
/// *fractional* floor first keeps the reference interior and preserves exact
/// extensive scaling because all operations occur after normalization.
fn reference_composition(
    target: InteriorSeedTarget,
    active: &[usize],
    input_moles: &[f64],
    scale: f64,
    fraction: f64,
) -> DVector<f64> {
    match target {
        InteriorSeedTarget::Uniform => {
            DVector::from_element(active.len(), 1.0 / active.len() as f64)
        }
        InteriorSeedTarget::InputAnchored => {
            let values = active
                .iter()
                .map(|&index| (input_moles[index] / scale).max(fraction))
                .collect::<Vec<_>>();
            let total = values.iter().sum::<f64>();
            DVector::from_iterator(values.len(), values.into_iter().map(|value| value / total))
        }
    }
}

fn validate_inputs(
    element_composition: &DMatrix<f64>,
    element_totals: &[f64],
    species_phase: &[usize],
    active_phase_mask: &[bool],
    input_moles: &[f64],
    settings: InteriorSeedSettings,
) -> Result<(), ReactionExtentError> {
    if element_composition.nrows() != input_moles.len()
        || species_phase.len() != input_moles.len()
        || element_composition.ncols() != element_totals.len()
    {
        return Err(ReactionExtentError::DimensionMismatch(format!(
            "interior seed has {} species, {} phase assignments, {} element rows, and {} totals",
            input_moles.len(),
            species_phase.len(),
            element_composition.ncols(),
            element_totals.len(),
        )));
    }
    if active_phase_mask.is_empty()
        || species_phase
            .iter()
            .any(|&phase| phase >= active_phase_mask.len())
    {
        return Err(ReactionExtentError::DimensionMismatch(
            "interior seed phase mask does not cover every species".to_owned(),
        ));
    }
    if !settings.requested_minimum_fraction.is_finite()
        || !settings.minimum_fraction.is_finite()
        || !settings.tolerance.is_finite()
        || !settings.trace_floor.is_finite()
        || settings.requested_minimum_fraction <= 0.0
        || settings.minimum_fraction <= 0.0
        || settings.minimum_fraction > settings.requested_minimum_fraction
        || settings.tolerance <= 0.0
        || settings.trace_floor <= 0.0
    {
        return Err(ReactionExtentError::InvalidProblem {
            field: "interior_seed_settings",
            message:
                "interior seed fractions, tolerance, and trace floor must be finite and positive"
                    .to_owned(),
        });
    }
    if input_moles
        .iter()
        .any(|moles| !moles.is_finite() || *moles < 0.0)
        || element_totals.iter().any(|total| !total.is_finite())
        || element_composition
            .iter()
            .any(|coefficient| !coefficient.is_finite())
    {
        return Err(ReactionExtentError::InvalidProblem {
            field: "interior_seed_input",
            message: "interior seed inputs must be finite, with non-negative moles".to_owned(),
        });
    }
    Ok(())
}

fn bounded_affine_projection(
    active_matrix: &DMatrix<f64>,
    element_target: &DVector<f64>,
    reference: &DVector<f64>,
    floor: f64,
    tolerance: f64,
) -> Result<Option<DVector<f64>>, ReactionExtentError> {
    let species_count = active_matrix.nrows();
    let element_count = active_matrix.ncols();
    let mut fixed = vec![false; species_count];

    for _ in 0..=species_count {
        let free = fixed
            .iter()
            .enumerate()
            .filter_map(|(index, &is_fixed)| (!is_fixed).then_some(index))
            .collect::<Vec<_>>();
        if free.is_empty() {
            return Ok(None);
        }
        let free_matrix = DMatrix::from_fn(element_count, free.len(), |row, column| {
            active_matrix[(free[column], row)]
        });
        let mut rhs = element_target.clone();
        for (index, &is_fixed) in fixed.iter().enumerate() {
            if is_fixed {
                for element in 0..element_count {
                    rhs[element] -= active_matrix[(index, element)] * floor;
                }
            }
        }
        let free_reference =
            DVector::from_iterator(free.len(), free.iter().map(|&index| reference[index]));
        let correction_rhs = &rhs - &free_matrix * &free_reference;
        let gram = &free_matrix * free_matrix.transpose();
        let multipliers = gram
            .svd(true, true)
            .solve(&correction_rhs, tolerance)
            .map_err(|message| ReactionExtentError::InvalidProblem {
                field: "element_feasible_interior_seed",
                message: format!("bounded affine projection SVD failed: {message}"),
            })?;
        let free_values = free_reference + free_matrix.transpose() * multipliers;
        let mut normalized = DVector::from_element(species_count, floor);
        for (local, &global) in free.iter().enumerate() {
            normalized[global] = free_values[local];
        }
        let residual = active_matrix.transpose() * &normalized - element_target;
        let max_residual = residual
            .iter()
            .fold(0.0_f64, |max, value| max.max(value.abs()));
        if max_residual > tolerance * 100.0 {
            return Ok(None);
        }
        let Some((worst, &minimum)) = normalized
            .iter()
            .enumerate()
            .filter(|(index, _)| !fixed[*index])
            .min_by(|(_, left), (_, right)| left.partial_cmp(right).unwrap())
        else {
            return Ok(None);
        };
        if minimum >= floor - tolerance {
            return Ok(Some(normalized));
        }
        fixed[worst] = true;
    }
    Ok(None)
}

fn max_element_balance_error(
    element_composition: &DMatrix<f64>,
    moles: &[f64],
    totals: &[f64],
) -> f64 {
    (element_composition.transpose() * DVector::from_column_slice(moles)
        - DVector::from_column_slice(totals))
    .iter()
    .fold(0.0_f64, |max, value| max.max(value.abs()))
}
