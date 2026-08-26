//! Canonical phase-qualified identity for equilibrium solver components.
//!
//! A component is not merely a substance name: `gas::H2O` and `liquid::H2O`
//! are distinct unknowns with potentially different standard-state records.
//! This module is deliberately independent from bridge construction and from
//! nonlinear solvers so both layers share one identity contract.

use crate::Thermodynamics::ChemEquilibrium::equilibrium_activity::PhaseActivityModel;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_ids::PhaseIndex;
use crate::Thermodynamics::User_PhaseOrSolution::PhaseModel;
use crate::Thermodynamics::phase_layout::{PhaseComponentId, PhaseId};
use crate::Thermodynamics::physical_state::PhysicalState;
use std::ops::Range;

/// Solver-facing identity and thermodynamic interpretation of one component.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct EquilibriumComponentDescriptor {
    id: PhaseComponentId,
    physical_state: PhysicalState,
    phase_model: PhaseModel,
    activity_model: PhaseActivityModel,
}

impl EquilibriumComponentDescriptor {
    /// Creates one fully qualified equilibrium component.
    pub fn new(
        id: PhaseComponentId,
        physical_state: PhysicalState,
        phase_model: PhaseModel,
        activity_model: PhaseActivityModel,
    ) -> Self {
        Self {
            id,
            physical_state,
            phase_model,
            activity_model,
        }
    }

    /// Phase-qualified component identity in canonical solver order.
    pub fn id(&self) -> &PhaseComponentId {
        &self.id
    }

    /// Bare substance name used for thermochemical record lookup.
    pub fn substance(&self) -> &str {
        &self.id.substance
    }

    /// Stable user-facing label that remains unique across phases.
    pub fn label(&self) -> String {
        self.id.label()
    }

    /// Physical state selected during thermochemical lookup.
    pub fn physical_state(&self) -> PhysicalState {
        self.physical_state
    }

    /// Domain-level phase model declared by the data subsystem.
    pub fn phase_model(&self) -> PhaseModel {
        self.phase_model
    }

    /// Canonical activity law consumed by equilibrium residuals.
    pub fn activity_model(&self) -> PhaseActivityModel {
        self.activity_model
    }
}

/// Compatibility conversion for existing one-phase gas callers.
///
/// New multiphase code must construct [`EquilibriumComponentDescriptor`]
/// explicitly. A bare string is projected into the anonymous single gas phase
/// only so old one-phase tests and examples keep a well-defined identity while
/// they migrate.
impl From<String> for EquilibriumComponentDescriptor {
    fn from(substance: String) -> Self {
        Self::new(
            PhaseComponentId::new(PhaseId::new(None), substance),
            PhysicalState::Gas,
            PhaseModel::IdealGas,
            PhaseActivityModel::IdealGas,
        )
    }
}

/// Ordered semantic description of one equilibrium phase.
///
/// Its component range is the canonical phase-to-component mapping. Numeric
/// `Phase { kind, species }` values exist only as a derived projection for
/// residual and Jacobian implementations that have not yet migrated.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct EquilibriumPhaseDescriptor {
    id: PhaseId,
    index: PhaseIndex,
    physical_state: PhysicalState,
    phase_model: PhaseModel,
    activity_model: PhaseActivityModel,
    component_range: Range<usize>,
}

impl EquilibriumPhaseDescriptor {
    /// Creates a fully qualified phase descriptor in canonical phase order.
    pub fn new(
        id: PhaseId,
        index: PhaseIndex,
        physical_state: PhysicalState,
        phase_model: PhaseModel,
        activity_model: PhaseActivityModel,
        component_range: Range<usize>,
    ) -> Self {
        Self {
            id,
            index,
            physical_state,
            phase_model,
            activity_model,
            component_range,
        }
    }

    /// Semantic phase identity retained through lookup, solve, and reporting.
    pub fn id(&self) -> &PhaseId {
        &self.id
    }

    /// Dense numeric phase coordinate derived from canonical order.
    pub fn index(&self) -> PhaseIndex {
        self.index
    }

    /// Physical state selected for thermochemical record lookup.
    pub fn physical_state(&self) -> PhysicalState {
        self.physical_state
    }

    /// Domain-level phase model.
    pub fn phase_model(&self) -> PhaseModel {
        self.phase_model
    }

    /// Canonical activity law used by the numerical formulation.
    pub fn activity_model(&self) -> PhaseActivityModel {
        self.activity_model
    }

    /// Contiguous component range in canonical component order.
    pub fn component_range(&self) -> Range<usize> {
        self.component_range.clone()
    }
}
