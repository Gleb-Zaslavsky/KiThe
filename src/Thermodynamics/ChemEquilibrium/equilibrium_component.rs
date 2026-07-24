//! Canonical phase-qualified identity for equilibrium solver components.
//!
//! A component is not merely a substance name: `gas::H2O` and `liquid::H2O`
//! are distinct unknowns with potentially different standard-state records.
//! This module is deliberately independent from bridge construction and from
//! nonlinear solvers so both layers share one identity contract.

use crate::Thermodynamics::ChemEquilibrium::equilibrium_activity::PhaseActivityModel;
use crate::Thermodynamics::User_PhaseOrSolution::PhaseModel;
use crate::Thermodynamics::phase_layout::{PhaseComponentId, PhaseId};
use crate::Thermodynamics::physical_state::PhysicalState;

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
