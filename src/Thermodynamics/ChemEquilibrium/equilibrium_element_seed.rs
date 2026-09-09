//! Element-conserving numerical seeds for elemental-input equilibrium problems.
//!
//! This module knows only the selected real-component matrix `A` and the
//! physical inventory `b`. It does not inspect Gibbs energies, choose phases,
//! or approximate equilibrium. Its sole job is to construct a positive
//! log-mole starting coordinate without turning a formal elemental carrier
//! into a fictional chemical species.

use nalgebra::{DMatrix, DVector};

use crate::Thermodynamics::ChemEquilibrium::equilibrium_nonlinear::{
    EquilibriumPreparationError, ReactionExtentError,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_problem::LogMolesInitialGuess;

/// Settings for the answer-independent feasible-seed projection.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ElementFeasibleSeedSettings {
    /// Requested minimum fraction for every selected component.
    pub requested_minimum_fraction: f64,
    /// Smallest fraction tried before the feasible interior is declared absent.
    pub minimum_fraction: f64,
    /// Relative linear-algebra tolerance in normalized coordinates.
    pub tolerance: f64,
    /// Positive coordinate used only for components that must stay physically zero.
    pub trace_floor: f64,
}

impl Default for ElementFeasibleSeedSettings {
    fn default() -> Self {
        Self {
            requested_minimum_fraction: 1.0e-6,
            minimum_fraction: 1.0e-12,
            tolerance: 1.0e-10,
            trace_floor: 1.0e-30,
        }
    }
}

/// A numerical seed whose physical moles conserve the supplied elemental
/// inventory before any log-coordinate trace floor is applied.
#[derive(Debug, Clone)]
pub struct ElementFeasibleSeed {
    physical_moles: Vec<f64>,
    log_moles: LogMolesInitialGuess,
    inventory_scale: f64,
    achieved_minimum_fraction: f64,
    max_element_balance_error: f64,
}

impl ElementFeasibleSeed {
    pub fn physical_moles(&self) -> &[f64] {
        &self.physical_moles
    }

    pub fn log_moles(&self) -> &LogMolesInitialGuess {
        &self.log_moles
    }

    pub fn inventory_scale(&self) -> f64 {
        self.inventory_scale
    }

    pub fn achieved_minimum_fraction(&self) -> f64 {
        self.achieved_minimum_fraction
    }

    pub fn max_element_balance_error(&self) -> f64 {
        self.max_element_balance_error
    }
}

/// Constructs a deterministic, answer-independent seed satisfying `A^T n=b`.
///
/// The projection starts from a uniform composition, applies a lower bound to
/// every real selected component, and relaxes that bound by decades only when
/// the requested interior is geometrically unavailable. It is a bounded affine
/// projection, not a Gibbs minimizer or an alternative equilibrium solver.
pub fn build_element_feasible_seed(
    element_composition: &DMatrix<f64>,
    element_totals: &[f64],
    settings: ElementFeasibleSeedSettings,
) -> Result<ElementFeasibleSeed, ReactionExtentError> {
    validate_inputs(element_composition, element_totals, settings)?;
    let scale = element_totals.iter().copied().fold(0.0_f64, f64::max);
    let target = DVector::from_iterator(
        element_totals.len(),
        element_totals.iter().map(|total| total / scale),
    );
    let reference = DVector::from_element(
        element_composition.nrows(),
        1.0 / element_composition.nrows() as f64,
    );
    // Any real species containing an absent element must be physically zero.
    // Leaving it merely "small" would alter the closed inventory before the
    // log-coordinate trace policy is even applied.
    let forced_zero = (0..element_composition.nrows())
        .map(|species| {
            element_totals.iter().enumerate().any(|(element, total)| {
                *total == 0.0 && element_composition[(species, element)] > 0.0
            })
        })
        .collect::<Vec<_>>();
    let mut floor = settings.requested_minimum_fraction;
    loop {
        if let Some(normalized) = bounded_affine_projection(
            element_composition,
            &target,
            &reference,
            &forced_zero,
            floor,
            settings.tolerance,
        )? {
            let physical_moles = normalized
                .iter()
                .map(|value| value * scale)
                .collect::<Vec<_>>();
            let log_moles =
                LogMolesInitialGuess::from_moles(&physical_moles, settings.trace_floor)?;
            let max_element_balance_error =
                max_balance_error(element_composition, &physical_moles, element_totals);
            return Ok(ElementFeasibleSeed {
                physical_moles,
                log_moles,
                inventory_scale: scale,
                achieved_minimum_fraction: normalized.iter().copied().fold(f64::INFINITY, f64::min),
                max_element_balance_error,
            });
        }
        if floor <= settings.minimum_fraction {
            break;
        }
        floor = (floor * 0.1).max(settings.minimum_fraction);
    }
    // A closed inventory can be feasible only on the boundary of the selected
    // species cone. That is a physical fact, not a failure of the element
    // input. Preserve those zero moles and let the existing log-coordinate
    // policy supply numerical trace values afterwards.
    if let Some(normalized) = bounded_affine_projection(
        element_composition,
        &target,
        &reference,
        &forced_zero,
        0.0,
        settings.tolerance,
    )? {
        let physical_moles = normalized
            .iter()
            .map(|value| value * scale)
            .collect::<Vec<_>>();
        let log_moles = LogMolesInitialGuess::from_moles(&physical_moles, settings.trace_floor)?;
        let max_element_balance_error =
            max_balance_error(element_composition, &physical_moles, element_totals);
        return Ok(ElementFeasibleSeed {
            physical_moles,
            log_moles,
            inventory_scale: scale,
            achieved_minimum_fraction: 0.0,
            max_element_balance_error,
        });
    }
    Err(ReactionExtentError::Preparation(
        EquilibriumPreparationError::NonRepresentableInventory {
            message: format!(
                "selected real species cannot provide an element-conserving positive seed down to fraction {:e}",
                settings.minimum_fraction
            ),
        },
    ))
}

fn validate_inputs(
    matrix: &DMatrix<f64>,
    totals: &[f64],
    settings: ElementFeasibleSeedSettings,
) -> Result<(), ReactionExtentError> {
    if matrix.nrows() == 0 {
        return Err(ReactionExtentError::Preparation(
            EquilibriumPreparationError::EmptySpeciesUniverse,
        ));
    }
    if matrix.ncols() == 0 || matrix.ncols() != totals.len() {
        return Err(ReactionExtentError::DimensionMismatch(format!(
            "element seed has {} species, {} element columns, and {} totals",
            matrix.nrows(),
            matrix.ncols(),
            totals.len()
        )));
    }
    if totals
        .iter()
        .any(|total| !total.is_finite() || *total < 0.0)
        || !totals.iter().any(|total| *total > 0.0)
        || matrix
            .iter()
            .any(|value| !value.is_finite() || *value < 0.0)
    {
        return Err(ReactionExtentError::InvalidProblem {
            field: "element_feasible_seed",
            message: "element totals must be non-negative with one positive entry, and A must be finite/non-negative".into(),
        });
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
            field: "element_feasible_seed_settings",
            message: "seed fractions, tolerance, and trace floor must be finite and positive"
                .into(),
        });
    }
    Ok(())
}

fn bounded_affine_projection(
    matrix: &DMatrix<f64>,
    target: &DVector<f64>,
    reference: &DVector<f64>,
    forced_zero: &[bool],
    floor: f64,
    tolerance: f64,
) -> Result<Option<DVector<f64>>, ReactionExtentError> {
    let mut fixed = forced_zero.to_vec();
    for _ in 0..=matrix.nrows() {
        let free = fixed
            .iter()
            .enumerate()
            .filter_map(|(index, fixed)| (!*fixed).then_some(index))
            .collect::<Vec<_>>();
        if free.is_empty() {
            return Ok(None);
        }
        let free_matrix = DMatrix::from_fn(matrix.ncols(), free.len(), |row, column| {
            matrix[(free[column], row)]
        });
        let mut rhs = target.clone();
        for (species, is_fixed) in fixed.iter().enumerate() {
            if *is_fixed {
                for element in 0..matrix.ncols() {
                    let fixed_value = if forced_zero[species] { 0.0 } else { floor };
                    rhs[element] -= matrix[(species, element)] * fixed_value;
                }
            }
        }
        let free_reference =
            DVector::from_iterator(free.len(), free.iter().map(|index| reference[*index]));
        let gram = &free_matrix * free_matrix.transpose();
        let multipliers = gram
            .svd(true, true)
            .solve(&(rhs - &free_matrix * &free_reference), tolerance)
            .map_err(|message| {
                ReactionExtentError::Preparation(
                    EquilibriumPreparationError::FeasibleSeedConstruction {
                        message: format!("bounded affine projection SVD failed: {message}"),
                    },
                )
            })?;
        let values = free_reference + free_matrix.transpose() * multipliers;
        let mut normalized = DVector::zeros(matrix.nrows());
        for (species, is_fixed) in fixed.iter().enumerate() {
            if *is_fixed && !forced_zero[species] {
                normalized[species] = floor;
            }
        }
        for (local, global) in free.iter().enumerate() {
            normalized[*global] = values[local];
        }
        let residual = matrix.transpose() * &normalized - target;
        if residual.iter().any(|value| !value.is_finite())
            || residual
                .iter()
                .fold(0.0_f64, |max, value| max.max(value.abs()))
                > tolerance * 100.0
        {
            return Ok(None);
        }
        let Some((worst, minimum)) = normalized
            .iter()
            .enumerate()
            .filter(|(index, _)| !fixed[*index])
            .min_by(|(_, left), (_, right)| left.partial_cmp(right).unwrap())
        else {
            return Ok(None);
        };
        if *minimum >= floor - tolerance {
            return Ok(Some(normalized));
        }
        fixed[worst] = true;
    }
    Ok(None)
}

fn max_balance_error(matrix: &DMatrix<f64>, moles: &[f64], totals: &[f64]) -> f64 {
    (matrix.transpose() * DVector::from_column_slice(moles) - DVector::from_column_slice(totals))
        .iter()
        .fold(0.0_f64, |max, value| max.max(value.abs()))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn feasible_seed_uses_only_real_species_and_preserves_inventory() {
        let matrix = DMatrix::from_row_slice(3, 2, &[2.0, 0.0, 0.0, 2.0, 2.0, 1.0]);
        let seed = build_element_feasible_seed(
            &matrix,
            &[2.0, 1.0],
            ElementFeasibleSeedSettings::default(),
        )
        .unwrap();
        assert!(seed.physical_moles().iter().all(|moles| *moles > 0.0));
        assert!(seed.max_element_balance_error() < 1.0e-8);
        assert_eq!(seed.log_moles().as_slice().len(), 3);
    }

    #[test]
    fn nonrepresentable_inventory_fails_before_a_solver_can_run() {
        let matrix = DMatrix::from_row_slice(1, 2, &[2.0, 0.0]);
        assert!(matches!(
            build_element_feasible_seed(
                &matrix,
                &[1.0, 1.0],
                ElementFeasibleSeedSettings::default()
            ),
            Err(ReactionExtentError::Preparation(
                EquilibriumPreparationError::NonRepresentableInventory { .. }
            ))
        ));
    }

    #[test]
    fn empty_species_universe_is_a_typed_preparation_failure() {
        let matrix = DMatrix::zeros(0, 2);
        assert!(matches!(
            build_element_feasible_seed(
                &matrix,
                &[1.0, 1.0],
                ElementFeasibleSeedSettings::default()
            ),
            Err(ReactionExtentError::Preparation(
                EquilibriumPreparationError::EmptySpeciesUniverse
            ))
        ));
    }

    #[test]
    fn boundary_only_inventory_preserves_physical_zeroes_before_log_trace_seeding() {
        let matrix = DMatrix::from_row_slice(2, 2, &[2.0, 0.0, 0.0, 2.0]);
        let seed = build_element_feasible_seed(
            &matrix,
            &[2.0, 0.0],
            ElementFeasibleSeedSettings::default(),
        )
        .unwrap();

        assert!((seed.physical_moles()[0] - 1.0).abs() < 1.0e-10);
        assert_eq!(seed.physical_moles()[1], 0.0);
        assert!(seed.log_moles().as_slice()[1].is_finite());
        assert_eq!(seed.achieved_minimum_fraction(), 0.0);
        assert!(seed.max_element_balance_error() < 1.0e-10);
    }

    #[test]
    fn feasible_seed_is_row_permutation_invariant_and_handles_extreme_scales() {
        let matrix = DMatrix::from_row_slice(3, 2, &[2.0, 0.0, 0.0, 2.0, 2.0, 1.0]);
        let permuted_matrix = DMatrix::from_row_slice(3, 2, &[2.0, 1.0, 2.0, 0.0, 0.0, 2.0]);
        let settings = ElementFeasibleSeedSettings::default();
        let reference = build_element_feasible_seed(&matrix, &[2.0, 1.0], settings).unwrap();
        let permuted =
            build_element_feasible_seed(&permuted_matrix, &[2.0, 1.0], settings).unwrap();

        for (actual, expected) in permuted.physical_moles().iter().zip([
            reference.physical_moles()[2],
            reference.physical_moles()[0],
            reference.physical_moles()[1],
        ]) {
            assert!((actual - expected).abs() < 1.0e-10);
        }

        for scale in [1.0e-8, 1.0e8] {
            let seed = build_element_feasible_seed(&matrix, &[2.0 * scale, 1.0 * scale], settings)
                .unwrap();
            assert_eq!(seed.inventory_scale(), 2.0 * scale);
            for (actual, expected) in seed.physical_moles().iter().zip(reference.physical_moles()) {
                assert!((actual / scale - expected).abs() < 1.0e-10);
            }
            assert!(seed.max_element_balance_error() / scale < 1.0e-8);
            assert!(
                seed.log_moles()
                    .as_slice()
                    .iter()
                    .all(|value| value.is_finite())
            );
        }
    }
}
