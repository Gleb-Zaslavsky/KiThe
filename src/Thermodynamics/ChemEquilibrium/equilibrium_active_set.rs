//! Immutable projection of one fixed phase set into nonlinear coordinates.
//!
//! Phase transitions are decided by the outer loop. This module freezes that
//! decision for one backend cascade and owns the only global/local index map
//! used to project seeds and scatter accepted solutions.

use crate::Thermodynamics::ChemEquilibrium::equilibrium_ids::{PhaseIndex, SpeciesId};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_log_moles::{Phase, species_to_phase_map};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_nonlinear::{
    ReactionBasis, ReactionExtentError, compute_reaction_basis,
};
use nalgebra::{DMatrix, DVector, linalg::SVD};

/// Prepared local coordinate system for one immutable active-phase set.
///
/// Phase transitions are decided by the outer loop (see
/// [`equilibrium_workflows`](super::equilibrium_workflows)). This module freezes
/// that decision for one backend cascade and owns the only global/local index map
/// used to project seeds and scatter accepted solutions.
#[derive(Debug, Clone)]
pub(crate) struct ActiveSetProjection {
    /// Active phases in the global phase ordering.
    pub active_phases: Vec<PhaseIndex>,
    /// Active species in the global species ordering.
    pub active_species: Vec<SpeciesId>,
    /// Phase descriptors for the active phases only, with local species indices.
    pub phases: Vec<Phase>,
    /// Maps each active species to its active phase index (local coordinates).
    pub species_phase: Vec<usize>,
    /// Element composition matrix restricted to active species only.
    /// Shape: `(active_species × elements)`.
    pub element_composition: DMatrix<f64>,
    /// SVD-derived reaction basis for the reduced active-species system.
    pub reaction_basis: ReactionBasis,
    /// Total number of species in the full (unprojected) system.
    full_species_count: usize,
}

impl ActiveSetProjection {
    pub fn build(
        phases: &[Phase],
        species_phase: &[usize],
        element_composition: &DMatrix<f64>,
        active: &[bool],
        rank_tolerance: f64,
    ) -> Result<Self, ReactionExtentError> {
        let species_count = species_phase.len();
        if active.len() != phases.len() || element_composition.nrows() != species_count {
            return Err(ReactionExtentError::DimensionMismatch(format!(
                "active-set projection has {} phase flags, {} phases, {} species labels, and {} element rows",
                active.len(),
                phases.len(),
                species_count,
                element_composition.nrows()
            )));
        }

        let active_phases = active
            .iter()
            .enumerate()
            .filter_map(|(phase, &is_active)| is_active.then_some(phase))
            .map(|phase| PhaseIndex::new(phase, phases.len()))
            .collect::<Result<Vec<_>, _>>()?;
        if active_phases.is_empty() {
            return Err(ReactionExtentError::InvalidProblem {
                field: "phase_active_set",
                message: "at least one phase must remain active".to_string(),
            });
        }

        let active_species = species_phase
            .iter()
            .enumerate()
            .filter_map(|(species, &phase)| active[phase].then_some(species))
            .map(|species| SpeciesId::new(species, species_count))
            .collect::<Result<Vec<_>, _>>()?;
        if active_species.is_empty() {
            return Err(ReactionExtentError::InvalidProblem {
                field: "phase_active_set",
                message: "active phases contain no species".to_string(),
            });
        }

        let mut full_to_local = vec![None; species_count];
        for (local, species) in active_species.iter().enumerate() {
            full_to_local[species.index()] = Some(local);
        }
        let mut local_phases = Vec::with_capacity(active_phases.len());
        for phase_id in &active_phases {
            let phase = &phases[phase_id.index()];
            let species = phase
                .species
                .iter()
                .map(|&full_species| {
                    full_to_local[full_species].ok_or_else(|| {
                        ReactionExtentError::InvalidProblem {
                            field: "phase_active_set",
                            message: format!(
                                "active phase {} contains species {full_species} outside its projection",
                                phase_id.index()
                            ),
                        }
                    })
                })
                .collect::<Result<Vec<_>, _>>()?;
            local_phases.push(Phase {
                kind: phase.kind,
                species,
            });
        }

        let local_elements = DMatrix::from_fn(
            active_species.len(),
            element_composition.ncols(),
            |row, col| element_composition[(active_species[row].index(), col)],
        );
        let reaction_basis = compute_reaction_basis(&local_elements, rank_tolerance)?;
        if reaction_basis.num_reactions + local_elements.ncols() != active_species.len() {
            return Err(ReactionExtentError::InvalidProblem {
                field: "phase_active_set",
                message: format!(
                    "active species span rank {} for {} conserved element columns; this phase set cannot represent the full closed-system inventory in the square log-moles formulation",
                    reaction_basis.rank,
                    local_elements.ncols()
                ),
            });
        }
        let local_species_phase = species_to_phase_map(&local_phases, active_species.len())?;

        Ok(Self {
            active_phases,
            active_species,
            phases: local_phases,
            species_phase: local_species_phase,
            element_composition: local_elements,
            reaction_basis,
            full_species_count: species_count,
        })
    }

    pub fn project_log_moles(&self, full: &[f64]) -> Result<Vec<f64>, ReactionExtentError> {
        if full.len() != self.full_species_count {
            return Err(ReactionExtentError::DimensionMismatch(format!(
                "cannot project {} log-moles into a {}-species active set",
                full.len(),
                self.full_species_count
            )));
        }
        Ok(self
            .active_species
            .iter()
            .map(|species| full[species.index()])
            .collect())
    }

    /// Verifies that the reduced active components can represent the closed
    /// system's original element inventory.
    ///
    /// The full problem stores species in rows and elements in columns, so a
    /// feasible active inventory requires `totals` to lie in the range of
    /// `A_active^T`. A trace seed cannot repair a missing elemental direction;
    /// reject that phase set before a nonlinear backend receives it.
    pub fn validate_element_totals_representable(
        &self,
        totals: &[f64],
        tolerance: f64,
    ) -> Result<(), ReactionExtentError> {
        if totals.len() != self.element_composition.ncols() {
            return Err(ReactionExtentError::DimensionMismatch(format!(
                "active-set feasibility has {} element totals for {} element columns",
                totals.len(),
                self.element_composition.ncols()
            )));
        }
        if !tolerance.is_finite() || tolerance <= 0.0 {
            return Err(ReactionExtentError::InvalidProblem {
                field: "phase_active_set",
                message: "feasibility tolerance must be finite and strictly positive".to_string(),
            });
        }
        if totals.iter().any(|value| !value.is_finite()) {
            return Err(ReactionExtentError::InvalidProblem {
                field: "phase_active_set",
                message: "original element totals must be finite".to_string(),
            });
        }

        let active_element_map = self.element_composition.transpose();
        let target = DVector::from_column_slice(totals);
        let svd = SVD::new(active_element_map.clone(), true, true);
        let coefficients = svd.solve(&target, tolerance).map_err(|message| {
            ReactionExtentError::InvalidProblem {
                field: "phase_active_set",
                message: format!("cannot solve active-set element feasibility system: {message}"),
            }
        })?;
        let residual = active_element_map * coefficients - target;
        let residual_norm = residual.norm();
        let allowed = tolerance * totals.iter().map(|value| value.abs()).sum::<f64>().max(1.0);
        if !residual_norm.is_finite() || residual_norm > allowed {
            return Err(ReactionExtentError::InvalidProblem {
                field: "phase_active_set",
                message: format!(
                    "active phases cannot represent the original element inventory: residual {residual_norm:e} exceeds {allowed:e}"
                ),
            });
        }
        Ok(())
    }

    pub fn scatter_log_moles(
        &self,
        local: &[f64],
        inactive_floor: f64,
    ) -> Result<Vec<f64>, ReactionExtentError> {
        if local.len() != self.active_species.len() {
            return Err(ReactionExtentError::DimensionMismatch(format!(
                "cannot scatter {} local values from {} active species",
                local.len(),
                self.active_species.len()
            )));
        }
        if !inactive_floor.is_finite() {
            return Err(ReactionExtentError::InvalidProblem {
                field: "inactive_log_mole_floor",
                message: "inactive log-mole floor must be finite".to_string(),
            });
        }
        let mut full = vec![inactive_floor; self.full_species_count];
        for (local_index, species) in self.active_species.iter().enumerate() {
            full[species.index()] = local[local_index];
        }
        Ok(full)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_log_moles::PhaseKind;

    #[test]
    fn projection_preserves_declared_order_and_scatter_uses_floor() {
        let phases = vec![
            Phase {
                kind: PhaseKind::IdealGas,
                species: vec![0, 1],
            },
            Phase {
                kind: PhaseKind::IdealSolution,
                species: vec![2],
            },
        ];
        let elements = DMatrix::from_row_slice(3, 2, &[1.0, 0.0, 0.0, 1.0, 1.0, 0.0]);
        let projection =
            ActiveSetProjection::build(&phases, &[0, 0, 1], &elements, &[true, false], 1e-10)
                .unwrap();

        assert_eq!(
            projection.active_species,
            vec![SpeciesId::new(0, 3).unwrap(), SpeciesId::new(1, 3).unwrap()]
        );
        assert_eq!(
            projection.project_log_moles(&[1.0, 2.0, 3.0]).unwrap(),
            vec![1.0, 2.0]
        );
        assert_eq!(
            projection.scatter_log_moles(&[4.0, 5.0], -700.0).unwrap(),
            vec![4.0, 5.0, -700.0]
        );
    }

    #[test]
    fn projection_rejects_an_empty_active_set() {
        let phases = vec![Phase {
            kind: PhaseKind::IdealGas,
            species: vec![0],
        }];
        let error = ActiveSetProjection::build(
            &phases,
            &[0],
            &DMatrix::from_row_slice(1, 1, &[1.0]),
            &[false],
            1e-10,
        )
        .unwrap_err();

        assert!(matches!(
            error,
            ReactionExtentError::InvalidProblem {
                field: "phase_active_set",
                ..
            }
        ));
    }

    #[test]
    fn projection_rejects_an_active_set_missing_a_conserved_element() {
        let phases = vec![
            Phase {
                kind: PhaseKind::IdealGas,
                species: vec![0],
            },
            Phase {
                kind: PhaseKind::IdealSolution,
                species: vec![1],
            },
        ];
        // Species zero contains only H; species one contains only C.
        let elements = DMatrix::from_row_slice(2, 2, &[1.0, 0.0, 0.0, 1.0]);
        let error = ActiveSetProjection::build(&phases, &[0, 1], &elements, &[true, false], 1e-10)
            .unwrap_err();
        assert!(matches!(
            error,
            ReactionExtentError::InvalidProblem {
                field: "phase_active_set",
                ..
            }
        ));
    }

    #[test]
    fn projection_keeps_global_order_and_round_trips_through_scatter() {
        let phases = vec![
            Phase {
                kind: PhaseKind::IdealGas,
                species: vec![0],
            },
            Phase {
                kind: PhaseKind::IdealSolution,
                species: vec![1, 2],
            },
            Phase {
                kind: PhaseKind::IdealSolution,
                species: vec![3],
            },
        ];
        let elements = DMatrix::from_row_slice(
            4,
            2,
            &[
                1.0, 0.0, //
                0.0, 1.0, //
                1.0, 1.0, //
                0.0, 1.0,
            ],
        );
        let projection = ActiveSetProjection::build(
            &phases,
            &[0, 1, 1, 2],
            &elements,
            &[true, false, true],
            1e-10,
        )
        .unwrap();

        assert_eq!(
            projection.active_phases,
            vec![
                PhaseIndex::new(0, 3).unwrap(),
                PhaseIndex::new(2, 3).unwrap(),
            ]
        );
        assert_eq!(
            projection.active_species,
            vec![SpeciesId::new(0, 4).unwrap(), SpeciesId::new(3, 4).unwrap()]
        );
        assert_eq!(projection.phases[0].species, vec![0]);
        assert_eq!(projection.phases[1].species, vec![1]);

        let local = projection
            .project_log_moles(&[10.0, 11.0, 12.0, 13.0])
            .unwrap();
        assert_eq!(local, vec![10.0, 13.0]);
        let scattered = projection.scatter_log_moles(&local, -700.0).unwrap();
        assert_eq!(scattered, vec![10.0, -700.0, -700.0, 13.0]);
    }

    #[test]
    fn projection_rejects_dependent_element_columns_in_the_active_set() {
        let phases = vec![
            Phase {
                kind: PhaseKind::IdealGas,
                species: vec![0],
            },
            Phase {
                kind: PhaseKind::IdealSolution,
                species: vec![1],
            },
        ];
        // Two identical element columns make the reduced basis rank-deficient.
        let elements = DMatrix::from_row_slice(2, 2, &[1.0, 1.0, 0.0, 0.0]);
        let error = ActiveSetProjection::build(&phases, &[0, 1], &elements, &[true, true], 1e-10)
            .unwrap_err();

        assert!(matches!(
            error,
            ReactionExtentError::InvalidProblem {
                field: "phase_active_set",
                ..
            }
        ));
    }

    #[test]
    fn projection_with_multiple_active_phases_keeps_local_phase_order() {
        let phases = vec![
            Phase {
                kind: PhaseKind::IdealGas,
                species: vec![0, 1],
            },
            Phase {
                kind: PhaseKind::IdealSolution,
                species: vec![2],
            },
            Phase {
                kind: PhaseKind::IdealSolution,
                species: vec![3, 4],
            },
        ];
        let elements = DMatrix::from_row_slice(
            5,
            3,
            &[
                1.0, 0.0, 0.0, //
                0.0, 1.0, 0.0, //
                0.0, 0.0, 1.0, //
                1.0, 1.0, 0.0, //
                0.0, 1.0, 1.0,
            ],
        );
        let projection = ActiveSetProjection::build(
            &phases,
            &[0, 0, 1, 2, 2],
            &elements,
            &[true, true, false],
            1e-10,
        )
        .unwrap();

        assert_eq!(projection.phases.len(), 2);
        assert_eq!(projection.phases[0].species, vec![0, 1]);
        assert_eq!(projection.phases[1].species, vec![2]);
        assert_eq!(
            projection.active_phases,
            vec![
                PhaseIndex::new(0, 3).unwrap(),
                PhaseIndex::new(1, 3).unwrap(),
            ]
        );
    }

    #[test]
    fn projection_handles_sparse_global_indices_round_trip_and_flooring() {
        let phases = vec![
            Phase {
                kind: PhaseKind::IdealGas,
                species: vec![0, 2],
            },
            Phase {
                kind: PhaseKind::IdealSolution,
                species: vec![1],
            },
            Phase {
                kind: PhaseKind::IdealSolution,
                species: vec![3, 4],
            },
        ];
        let elements = DMatrix::from_row_slice(
            5,
            3,
            &[
                1.0, 0.0, 0.0, //
                0.0, 1.0, 0.0, //
                1.0, 1.0, 0.0, //
                0.0, 0.0, 1.0, //
                0.0, 1.0, 1.0,
            ],
        );
        let projection = ActiveSetProjection::build(
            &phases,
            &[0, 1, 0, 2, 2],
            &elements,
            &[true, false, true],
            1e-10,
        )
        .unwrap();

        assert_eq!(
            projection.active_phases,
            vec![
                PhaseIndex::new(0, 3).unwrap(),
                PhaseIndex::new(2, 3).unwrap(),
            ]
        );
        assert_eq!(
            projection.active_species,
            vec![
                SpeciesId::new(0, 5).unwrap(),
                SpeciesId::new(2, 5).unwrap(),
                SpeciesId::new(3, 5).unwrap(),
                SpeciesId::new(4, 5).unwrap()
            ]
        );
        assert_eq!(
            projection
                .project_log_moles(&[10.0, 11.0, 12.0, 13.0, 14.0])
                .unwrap(),
            vec![10.0, 12.0, 13.0, 14.0]
        );
        assert_eq!(
            projection
                .scatter_log_moles(&[20.0, 21.0, 22.0, 23.0], -777.0)
                .unwrap(),
            vec![20.0, -777.0, 21.0, 22.0, 23.0]
        );
    }

    #[test]
    fn projection_round_trips_projected_local_values_without_reordering() {
        let phases = vec![
            Phase {
                kind: PhaseKind::IdealGas,
                species: vec![0, 2],
            },
            Phase {
                kind: PhaseKind::IdealSolution,
                species: vec![1, 3],
            },
        ];
        let elements = DMatrix::from_row_slice(
            4,
            2,
            &[
                1.0, 0.0, //
                0.0, 1.0, //
                1.0, 1.0, //
                0.0, 0.0,
            ],
        );
        let projection =
            ActiveSetProjection::build(&phases, &[0, 1, 0, 1], &elements, &[true, true], 1e-10)
                .unwrap();

        let local = vec![0.5, 1.5, 2.5, 3.5];
        let full = projection.scatter_log_moles(&local, -999.0).unwrap();
        assert_eq!(full, vec![0.5, 1.5, 2.5, 3.5]);
        assert_eq!(projection.project_log_moles(&full).unwrap(), local);
    }

    #[test]
    fn projection_reports_the_expected_reduced_basis_dimensions() {
        let phases = vec![Phase {
            kind: PhaseKind::IdealGas,
            species: vec![0, 1, 2],
        }];
        let elements = DMatrix::from_row_slice(
            3,
            2,
            &[
                1.0, 0.0, //
                0.0, 1.0, //
                1.0, 1.0,
            ],
        );
        let projection =
            ActiveSetProjection::build(&phases, &[0, 0, 0], &elements, &[true], 1e-10).unwrap();

        assert_eq!(projection.reaction_basis.rank, 2);
        assert_eq!(projection.reaction_basis.num_reactions, 1);
        assert_eq!(projection.active_species.len(), 3);
        assert_eq!(projection.phases[0].species, vec![0, 1, 2]);
    }

    #[test]
    fn projection_rejects_an_impossible_sparse_active_set_before_backend_use() {
        let phases = vec![
            Phase {
                kind: PhaseKind::IdealGas,
                species: vec![0, 3],
            },
            Phase {
                kind: PhaseKind::IdealSolution,
                species: vec![1],
            },
            Phase {
                kind: PhaseKind::IdealSolution,
                species: vec![2],
            },
        ];
        let elements = DMatrix::from_row_slice(
            4,
            2,
            &[
                1.0, 0.0, //
                0.0, 1.0, //
                0.0, 1.0, //
                1.0, 0.0,
            ],
        );
        let error = ActiveSetProjection::build(
            &phases,
            &[0, 1, 2, 0],
            &elements,
            &[true, false, false],
            1e-10,
        )
        .unwrap_err();

        assert!(matches!(
            error,
            ReactionExtentError::InvalidProblem {
                field: "phase_active_set",
                ..
            }
        ));
    }
}
