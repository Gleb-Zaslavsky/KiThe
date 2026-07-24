//! Typed declarations and resolved records for phase thermodynamics.
//!
//! This module deliberately contains no numeric property evaluation and no
//! cache state. It defines what a phase is before lookup, then validates the
//! layout-aligned payload produced by resolution.
//!
//! ```text
//! HIERARCHY
//!
//!             PhaseSpec
//!       /        |            \               
//! PhaseId PhasePhysicalState PhaseModel
//! ```

use std::collections::{HashMap, HashSet};
use std::fmt;

use crate::Thermodynamics::User_substances::{Phases, SubsData};
use crate::Thermodynamics::User_substances2::SearchSummaryReport;
use crate::Thermodynamics::phase_layout::{PhaseId, SystemLayout};
use crate::Thermodynamics::physical_state::PhysicalState;

use super::SubstanceSystemFactoryError;

/// Physical state of a phase. This is deliberately separate from the
/// thermodynamic model: a liquid is not automatically an ideal solution.
/// Backwards-compatible phase-level name for the canonical lookup state.
pub type PhasePhysicalState = PhysicalState;

/// Activity/mixing model actually implemented by the phase-property layer.
/// More elaborate solution models must be added with their standard-state and
/// activity-coefficient contracts; they must not be implied by a phase name.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PhaseModel {
    IdealGas,
    PureCondensed,
}

/// Declarative definition of one phase before library data is resolved.
/// Components are ordered because that order is part of the numerical solver
/// contract and is preserved by `SystemLayout` after resolution.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PhaseSpec {
    id: PhaseId,
    components: Vec<String>,
    physical_state: PhasePhysicalState,
    model: PhaseModel,
}

impl PhaseSpec {
    pub fn new(
        id: PhaseId,
        components: Vec<String>,
        physical_state: PhasePhysicalState,
        model: PhaseModel,
    ) -> Result<Self, SubstanceSystemFactoryError> {
        let spec = Self {
            id,
            components,
            physical_state,
            model,
        };
        spec.validate()?;
        Ok(spec)
    }

    pub fn ideal_gas(
        id: PhaseId,
        components: Vec<String>,
    ) -> Result<Self, SubstanceSystemFactoryError> {
        Self::new(
            id,
            components,
            PhasePhysicalState::Gas,
            PhaseModel::IdealGas,
        )
    }

    pub fn pure_condensed(
        id: PhaseId,
        components: Vec<String>,
        physical_state: PhasePhysicalState,
    ) -> Result<Self, SubstanceSystemFactoryError> {
        Self::new(id, components, physical_state, PhaseModel::PureCondensed)
    }

    pub fn id(&self) -> &PhaseId {
        &self.id
    }

    pub fn components(&self) -> &[String] {
        &self.components
    }

    pub fn physical_state(&self) -> PhasePhysicalState {
        self.physical_state
    }

    pub fn model(&self) -> PhaseModel {
        self.model
    }

    /// Derives the per-substance lookup-state projection required by `SubsData`.
    ///
    /// `PhaseSpec` remains the canonical owner of physical state. This map is
    /// a compatibility-shaped projection for lower-level record lookup, not
    /// an independent declaration that callers should maintain by hand.
    pub(crate) fn legacy_component_phase_map(&self) -> HashMap<String, Option<Phases>> {
        let state = Some(self.physical_state);
        self.components
            .iter()
            .cloned()
            .map(|component| (component, state))
            .collect()
    }

    pub(crate) fn validate(&self) -> Result<(), SubstanceSystemFactoryError> {
        if let Some(name) = self.id.as_option() {
            if name.trim().is_empty() {
                return Err(SubstanceSystemFactoryError::InvalidSpecification {
                    field: "phase id".to_string(),
                    reason: "named phases must have a non-empty identifier".to_string(),
                });
            }
        }
        if self.components.is_empty() {
            return Err(SubstanceSystemFactoryError::InvalidSpecification {
                field: "phase components".to_string(),
                reason: "each phase must contain at least one substance".to_string(),
            });
        }
        let mut seen = HashSet::with_capacity(self.components.len());
        for component in &self.components {
            if component.trim().is_empty() || !seen.insert(component) {
                return Err(SubstanceSystemFactoryError::InvalidSpecification {
                    field: "phase components".to_string(),
                    reason: format!(
                        "component '{}' is empty or duplicated within phase",
                        component
                    ),
                });
            }
        }
        match (self.physical_state, self.model) {
            (PhasePhysicalState::Gas, PhaseModel::IdealGas)
            | (PhasePhysicalState::Liquid, PhaseModel::PureCondensed)
            | (PhasePhysicalState::Solid, PhaseModel::PureCondensed)
            | (PhasePhysicalState::Condensed, PhaseModel::PureCondensed) => Ok(()),
            (state, model) => Err(SubstanceSystemFactoryError::InvalidSpecification {
                field: "phase model".to_string(),
                reason: format!(
                    "model {:?} is not implemented for physical state {:?}",
                    model, state
                ),
            }),
        }
    }
}

/// Fully resolved phase-system input. It connects validated declarations to
/// the exact `SubsData` records selected from the thermodynamic libraries and
/// to the solver component order derived from those records.
///
/// This is intentionally distinct from `PhaseSystem`, which owns mutable
/// computed caches. Resolution is a fallible boundary; evaluation is not
/// allowed to silently change what was resolved.
///
/// The accompanying report preserves the lookup provenance selected for each
/// phase, including whether NIST fallback was enabled for this resolution.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PhaseResolutionSummary {
    phase: PhaseId,
    search: SearchSummaryReport,
}

impl PhaseResolutionSummary {
    pub fn phase(&self) -> &PhaseId {
        &self.phase
    }

    pub fn search(&self) -> &SearchSummaryReport {
        &self.search
    }
}

/// Immutable lookup provenance paired with a resolved phase-system payload.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ResolvedPhaseSystemReport {
    nist_fallback_enabled: bool,
    phases: Vec<PhaseResolutionSummary>,
}

impl ResolvedPhaseSystemReport {
    /// Whether unresolved local records were allowed to query NIST.
    pub fn nist_fallback_enabled(&self) -> bool {
        self.nist_fallback_enabled
    }

    /// Per-phase lookup summaries in canonical solver phase order.
    pub fn phases(&self) -> &[PhaseResolutionSummary] {
        &self.phases
    }

    /// Finds lookup provenance for one phase identity.
    pub fn phase(&self, phase: &PhaseId) -> Option<&PhaseResolutionSummary> {
        self.phases.iter().find(|summary| summary.phase == *phase)
    }
}

#[derive(Clone)]
pub struct ResolvedPhaseSystem {
    phase_specs: Vec<PhaseSpec>,
    phase_data: HashMap<Option<String>, SubsData>,
    layout: SystemLayout,
    report: ResolvedPhaseSystemReport,
}

impl ResolvedPhaseSystem {
    pub fn new(
        phase_specs: Vec<PhaseSpec>,
        phase_data: HashMap<Option<String>, SubsData>,
    ) -> Result<Self, SubstanceSystemFactoryError> {
        Self::new_with_nist_fallback_policy(phase_specs, phase_data, false)
    }

    /// Creates a resolved system while retaining the lookup policy that
    /// produced the phase payloads.
    pub(crate) fn new_with_nist_fallback_policy(
        mut phase_specs: Vec<PhaseSpec>,
        phase_data: HashMap<Option<String>, SubsData>,
        nist_fallback_enabled: bool,
    ) -> Result<Self, SubstanceSystemFactoryError> {
        phase_specs.sort_by(|left, right| left.id.cmp(&right.id));
        let expected_keys = phase_specs
            .iter()
            .map(|spec| spec.id.as_option().clone())
            .collect::<HashSet<_>>();
        let resolved_keys = phase_data.keys().cloned().collect::<HashSet<_>>();
        if expected_keys != resolved_keys {
            return Err(SubstanceSystemFactoryError::InvalidSpecification {
                field: "resolved phases".to_string(),
                reason: "resolved phase records do not match the requested phase ids".to_string(),
            });
        }
        for phase in &phase_specs {
            let key = phase.id.as_option();
            let data = phase_data.get(key).ok_or_else(|| {
                SubstanceSystemFactoryError::InvalidSpecification {
                    field: "resolved phases".to_string(),
                    reason: format!("missing resolved record for phase {:?}", phase.id),
                }
            })?;
            if data.substances != phase.components {
                return Err(SubstanceSystemFactoryError::InvalidSpecification {
                    field: "resolved phase components".to_string(),
                    reason: format!(
                        "resolved component order for phase {:?} differs from its specification",
                        phase.id
                    ),
                });
            }
        }
        let layout = SystemLayout::from_phase_map(&phase_data);
        let phase_reports = phase_specs
            .iter()
            .map(|phase| {
                let key = phase.id.as_option();
                let data = phase_data.get(key).ok_or_else(|| {
                    SubstanceSystemFactoryError::InvalidSpecification {
                        field: "resolved phases".to_string(),
                        reason: format!(
                            "missing resolved record while reporting phase {:?}",
                            phase.id
                        ),
                    }
                })?;
                Ok(PhaseResolutionSummary {
                    phase: phase.id.clone(),
                    search: data.search_summary_report(),
                })
            })
            .collect::<Result<Vec<_>, SubstanceSystemFactoryError>>()?;
        let report = ResolvedPhaseSystemReport {
            nist_fallback_enabled,
            phases: phase_reports,
        };
        Ok(Self {
            phase_specs,
            phase_data,
            layout,
            report,
        })
    }

    pub fn phase_specs(&self) -> &[PhaseSpec] {
        &self.phase_specs
    }

    pub fn layout(&self) -> &SystemLayout {
        &self.layout
    }

    pub fn phase_data(&self) -> &HashMap<Option<String>, SubsData> {
        &self.phase_data
    }

    /// Typed lookup provenance for the resolved payloads.
    pub fn report(&self) -> &ResolvedPhaseSystemReport {
        &self.report
    }

    pub(crate) fn into_parts(
        self,
    ) -> (
        Vec<PhaseSpec>,
        HashMap<Option<String>, SubsData>,
        SystemLayout,
        ResolvedPhaseSystemReport,
    ) {
        (self.phase_specs, self.phase_data, self.layout, self.report)
    }
}

impl fmt::Debug for ResolvedPhaseSystem {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("ResolvedPhaseSystem")
            .field("phase_specs", &self.phase_specs)
            .field("layout", &self.layout)
            .finish_non_exhaustive()
    }
}
