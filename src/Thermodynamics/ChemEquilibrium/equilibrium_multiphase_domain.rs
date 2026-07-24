//! Typed physical boundary for a multiphase equilibrium problem.
//!
//! The legacy solver still accepts indexed `Phase` values.  This module is the
//! boundary before that projection: phase identities are semantic, components
//! are phase-qualified, and physical initial moles exist before a positive
//! log-space trace seed is selected.
//!
//! ```text
//! PhaseSpec[] --> MultiphaseEquilibriumLayout --> solver index projection
//!                        |
//!                        +--> MultiphaseInitialComposition
//!                                  |
//!                                  +--> physical element totals
//! ```

use crate::Thermodynamics::ChemEquilibrium::equilibrium_nonlinear::ReactionExtentError;
use crate::Thermodynamics::User_PhaseOrSolution::{PhaseModel, PhaseSpec};
use crate::Thermodynamics::phase_layout::{PhaseComponentId, SystemLayout};
use crate::Thermodynamics::physical_state::PhysicalState;
use nalgebra::DMatrix;
use std::collections::HashSet;

/// Validated phase/component order for a closed fixed-pressure, fixed-temperature
/// equilibrium problem.
///
/// Currently implemented activity models are one ideal gas phase and any
/// number of one-component pure condensed phases.  A liquid or solid solution
/// needs its own standard-state and activity-coefficient contract and is
/// rejected here rather than silently treated as ideal.
/// Validated phase/component order for a closed fixed-pressure, fixed-temperature
/// equilibrium problem.
///
/// Currently implemented activity models are one ideal gas phase and any
/// number of one-component pure condensed phases. A liquid or solid solution
/// needs its own standard-state and activity-coefficient contract and is
/// rejected here rather than silently treated as ideal.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MultiphaseEquilibriumLayout {
    /// Phase specifications in canonical semantic solver order.
    phase_specs: Vec<PhaseSpec>,
    /// Resolved system layout with component ordering and phase ranges.
    system_layout: SystemLayout,
    /// Stable process-independent fingerprint of the declared layout.
    fingerprint: u64,
}

impl MultiphaseEquilibriumLayout {
    /// Validates physical model support and creates a deterministic component order.
    pub fn new(mut phase_specs: Vec<PhaseSpec>) -> Result<Self, ReactionExtentError> {
        if phase_specs.is_empty() {
            return Err(ReactionExtentError::InvalidProblem {
                field: "phase_specs",
                message: "a multiphase equilibrium problem requires at least one phase".into(),
            });
        }

        // Resolution already uses the semantic phase order. Apply the same
        // order here so direct construction and resolved-system construction
        // cannot assign different numeric component indices to the same input.
        phase_specs.sort_by(|left, right| left.id().cmp(right.id()));

        let mut phase_ids = HashSet::with_capacity(phase_specs.len());
        let mut gas_phase_count = 0usize;

        for phase in &phase_specs {
            let phase_id = phase.id().clone();
            if !phase_ids.insert(phase_id.clone()) {
                return Err(ReactionExtentError::InvalidProblem {
                    field: "phase_specs",
                    message: format!("duplicate phase identity {:?}", phase_id.as_option()),
                });
            }

            match (phase.physical_state(), phase.model()) {
                (PhysicalState::Gas, PhaseModel::IdealGas) => gas_phase_count += 1,
                (_, PhaseModel::PureCondensed) if phase.components().len() == 1 => {}
                (_, PhaseModel::PureCondensed) => {
                    return Err(ReactionExtentError::ValidationNotApplicable {
                        path: "multiphase equilibrium layout",
                        message: format!(
                            "pure condensed phase {:?} must contain exactly one component",
                            phase_id.as_option()
                        ),
                    });
                }
                (state, model) => {
                    return Err(ReactionExtentError::ValidationNotApplicable {
                        path: "multiphase equilibrium layout",
                        message: format!(
                            "activity model {:?} is not supported for physical state {:?}",
                            model, state
                        ),
                    });
                }
            }
        }

        if gas_phase_count > 1 {
            return Err(ReactionExtentError::ValidationNotApplicable {
                path: "multiphase equilibrium layout",
                message: "at most one ideal-gas phase is currently implemented".into(),
            });
        }

        let system_layout = layout_from_ordered_specs(&phase_specs);
        let fingerprint = fingerprint_specs(&phase_specs);
        Ok(Self {
            phase_specs,
            system_layout,
            fingerprint,
        })
    }

    /// Phase specifications in canonical semantic solver order.
    pub fn phase_specs(&self) -> &[PhaseSpec] {
        &self.phase_specs
    }

    /// Phase-qualified solver-component ordering.
    pub fn system_layout(&self) -> &SystemLayout {
        &self.system_layout
    }

    /// Stable process-independent fingerprint of the declared layout.
    pub fn fingerprint(&self) -> u64 {
        self.fingerprint
    }

    /// Number of phase-qualified components.
    pub fn component_count(&self) -> usize {
        self.system_layout.component_count()
    }
}

/// Physical mole numbers aligned to a [`MultiphaseEquilibriumLayout`].
///
/// A zero is meaningful here: the component is absent from the physical
/// initial composition but remains a candidate for equilibrium.  Conversion
/// to a positive log-mole seed belongs to the numerical boundary and must not
/// change this value object.
#[derive(Debug, Clone, PartialEq)]
/// Physical mole numbers aligned to a [`MultiphaseEquilibriumLayout`].
///
/// A zero is meaningful here: the component is absent from the physical
/// initial composition but remains a candidate for equilibrium. Conversion
/// to a positive log-mole seed belongs to the numerical boundary and must not
/// change this value object.
pub struct MultiphaseInitialComposition {
    /// Fingerprint of the layout this composition was created for.
    /// Used to detect mismatches between composition and layout at runtime.
    layout_fingerprint: u64,
    /// Physical mole numbers in canonical layout component order.
    /// Zero values indicate absent components.
    moles: Vec<f64>,
}

impl MultiphaseInitialComposition {
    /// Creates a composition from values in canonical layout order.
    pub fn from_dense(
        layout: &MultiphaseEquilibriumLayout,
        moles: Vec<f64>,
    ) -> Result<Self, ReactionExtentError> {
        if moles.len() != layout.component_count() {
            return Err(ReactionExtentError::DimensionMismatch(format!(
                "initial composition has {} entries but layout has {} components",
                moles.len(),
                layout.component_count()
            )));
        }
        validate_moles(&moles)?;
        Ok(Self {
            layout_fingerprint: layout.fingerprint(),
            moles,
        })
    }

    /// Creates a composition from explicit phase-qualified component amounts.
    ///
    /// Missing layout components are physical zeroes, never an implicit phase
    /// exclusion. Duplicate sparse entries are rejected to avoid ambiguous
    /// accidental accumulation.
    pub fn from_sparse(
        layout: &MultiphaseEquilibriumLayout,
        entries: Vec<(PhaseComponentId, f64)>,
    ) -> Result<Self, ReactionExtentError> {
        let mut moles = vec![0.0; layout.component_count()];
        let mut assigned = HashSet::with_capacity(entries.len());
        for (component, amount) in entries {
            let Some(index) = layout.system_layout().component_index(&component) else {
                return Err(ReactionExtentError::InvalidProblem {
                    field: "initial_composition",
                    message: format!(
                        "component '{}' is not in the phase layout",
                        component.label()
                    ),
                });
            };
            if !assigned.insert(component.clone()) {
                return Err(ReactionExtentError::InvalidProblem {
                    field: "initial_composition",
                    message: format!(
                        "component '{}' is specified more than once",
                        component.label()
                    ),
                });
            }
            moles[index] = amount;
        }
        Self::from_dense(layout, moles)
    }

    /// Physical mole numbers in canonical component order.
    pub fn moles(&self) -> &[f64] {
        &self.moles
    }

    /// Confirms that this composition belongs to the supplied layout.
    pub fn validate_for(
        &self,
        layout: &MultiphaseEquilibriumLayout,
    ) -> Result<(), ReactionExtentError> {
        if self.layout_fingerprint != layout.fingerprint() {
            return Err(ReactionExtentError::InvalidProblem {
                field: "initial_composition",
                message: "composition belongs to a different multiphase layout".into(),
            });
        }
        Ok(())
    }

    /// Computes conserved physical element totals before any trace-mole seed
    /// is introduced for the logarithmic nonlinear coordinate.
    ///
    /// `element_composition` follows the canonical `(component, element)`
    /// orientation used by `EquilibriumProblem`.
    pub fn element_totals(
        &self,
        layout: &MultiphaseEquilibriumLayout,
        element_composition: &DMatrix<f64>,
    ) -> Result<Vec<f64>, ReactionExtentError> {
        self.validate_for(layout)?;
        if element_composition.nrows() != self.moles.len() {
            return Err(ReactionExtentError::DimensionMismatch(format!(
                "element composition has {} component rows but composition has {} entries",
                element_composition.nrows(),
                self.moles.len()
            )));
        }
        if element_composition.iter().any(|value| !value.is_finite()) {
            return Err(ReactionExtentError::InvalidProblem {
                field: "element_composition",
                message: "element composition must contain only finite values".into(),
            });
        }
        let mut totals = vec![0.0; element_composition.ncols()];
        for component in 0..element_composition.nrows() {
            for element in 0..element_composition.ncols() {
                totals[element] +=
                    self.moles[component] * element_composition[(component, element)];
            }
        }
        Ok(totals)
    }
}

fn validate_moles(moles: &[f64]) -> Result<(), ReactionExtentError> {
    let mut has_positive_amount = false;
    for (index, amount) in moles.iter().copied().enumerate() {
        if !amount.is_finite() || amount < 0.0 {
            return Err(ReactionExtentError::InvalidProblem {
                field: "initial_composition",
                message: format!("component {index} must have a finite non-negative amount"),
            });
        }
        has_positive_amount |= amount > 0.0;
    }
    if !has_positive_amount {
        return Err(ReactionExtentError::InvalidProblem {
            field: "initial_composition",
            message: "at least one component must have a positive physical mole amount".into(),
        });
    }
    Ok(())
}

fn layout_from_ordered_specs(specs: &[PhaseSpec]) -> SystemLayout {
    SystemLayout::from_ordered_phase_components(
        specs
            .iter()
            .map(|spec| (spec.id().clone(), spec.components().to_vec()))
            .collect(),
    )
}

fn fingerprint_specs(specs: &[PhaseSpec]) -> u64 {
    const FNV_OFFSET: u64 = 0xcbf29ce484222325;
    const FNV_PRIME: u64 = 0x00000100000001b3;
    let mut hash = FNV_OFFSET;
    let mut write = |bytes: &[u8]| {
        for byte in bytes {
            hash ^= u64::from(*byte);
            hash = hash.wrapping_mul(FNV_PRIME);
        }
    };
    for spec in specs {
        match spec.id().as_option() {
            Some(id) => write(id.as_bytes()),
            None => write(b"<anonymous>"),
        }
        write(&[0]);
        write(&[match spec.physical_state() {
            PhysicalState::Gas => 1,
            PhysicalState::Liquid => 2,
            PhysicalState::Solid => 3,
            PhysicalState::Condensed => 4,
        }]);
        write(&[match spec.model() {
            PhaseModel::IdealGas => 1,
            PhaseModel::PureCondensed => 2,
        }]);
        for component in spec.components() {
            write(component.as_bytes());
            write(&[0]);
        }
        write(&[0xff]);
    }
    hash
}
