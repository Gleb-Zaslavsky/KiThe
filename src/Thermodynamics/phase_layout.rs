//! Canonical ordering for components of a resolved phase system.
//!
//! `SystemLayout` establishes the phase sequence, the substance sequence
//! within each phase, and the phase-qualified identity of every solver
//! component.
//!
//! ```text
//!                 SystemLayout
//!            /      |         \
//!          phases   components  phase_ranges
//!                      |
//!              PhaseComponentId
//! ```
//!
//! A thermodynamic component is identified by **both** its phase and
//! substance. Consequently, the same substance in gas and liquid is two
//! distinct solver components. `SystemLayout` is the single ordering contract
//! used by numeric vectors, symbolic variables, element matrices, and cached
//! thermodynamic snapshots.
//!
//! `HashMap`-backed phase payloads are always projected through the helpers in
//! this module before they reach a solver-facing path. Named phases are sorted
//! deterministically; the anonymous `None` phase represents a single-phase
//! system. Its deterministic order keeps numeric vectors, symbolic variables,
//! and phase-keyed cache maps aligned without depending on `HashMap`
//! iteration order.

use crate::Thermodynamics::User_substances::SubsData;
use std::collections::HashMap;
use std::ops::Range;

#[derive(Debug, Clone, PartialEq, Eq, Hash, PartialOrd, Ord)]
/// Stable identity of one thermodynamic phase.
///
/// `None` is reserved for a one-phase system; named phases use `Some(name)`.
pub struct PhaseId(pub Option<String>);

impl PhaseId {
    /// Creates a phase identity from the legacy optional phase name.
    pub fn new(phase: Option<String>) -> Self {
        Self(phase)
    }

    /// Borrows the legacy optional name without allocating.
    pub fn as_option(&self) -> &Option<String> {
        &self.0
    }

    /// Consumes the identity and returns the legacy optional name.
    pub fn into_option(self) -> Option<String> {
        self.0
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
/// Solver identity of a component: one substance in one phase.
pub struct PhaseComponentId {
    /// Phase that owns the component.
    pub phase: PhaseId,
    /// Canonical substance name within the phase.
    pub substance: String,
}

impl PhaseComponentId {
    /// Creates a component identity from its phase and substance name.
    pub fn new(phase: PhaseId, substance: impl Into<String>) -> Self {
        Self {
            phase,
            substance: substance.into(),
        }
    }

    /// Human-readable, collision-free component label.
    ///
    /// Named phases use `phase::substance`; single-phase labels are simply the
    /// substance name for compatibility with existing output.
    pub fn label(&self) -> String {
        match self.phase.as_option() {
            Some(phase) => format!("{}::{}", phase, self.substance),
            None => self.substance.clone(),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
/// Ordered, immutable component/phase projection used by solver-facing code.
///
/// Public vectors are ordered by phase and then by the component order declared
/// for that phase. `phase_ranges` maps each phase to its contiguous component
/// slice; callers should prefer the accessor methods rather than reconstructing
/// ranges from maps.
pub struct SystemLayout {
    /// Canonical phase order.
    pub phases: Vec<PhaseId>,
    /// Canonical component order aligned with vectors and symbolic variables.
    pub components: Vec<PhaseComponentId>,
    /// Contiguous component ranges for each phase.
    phase_ranges: HashMap<PhaseId, Range<usize>>,
}

impl SystemLayout {
    /// Builds a layout from already validated, explicitly ordered phase
    /// component groups.
    ///
    /// Unlike [`Self::from_phase_map`], this constructor never sorts phases:
    /// declaration order is preserved because callers already own a typed
    /// ordering contract (for example a multiphase equilibrium problem).
    pub fn from_ordered_phase_components(phase_components: Vec<(PhaseId, Vec<String>)>) -> Self {
        let mut phases = Vec::with_capacity(phase_components.len());
        let mut components = Vec::new();
        let mut phase_ranges = HashMap::with_capacity(phase_components.len());

        for (phase, substances) in phase_components {
            let start = components.len();
            components.extend(
                substances
                    .into_iter()
                    .map(|substance| PhaseComponentId::new(phase.clone(), substance)),
            );
            phase_ranges.insert(phase.clone(), start..components.len());
            phases.push(phase);
        }

        Self {
            phases,
            components,
            phase_ranges,
        }
    }

    /// Builds the canonical layout for a one-phase payload.
    pub fn from_single_phase(substances: &[String]) -> Self {
        let phase = PhaseId::new(None);
        let components = substances
            .iter()
            .cloned()
            .map(|substance| PhaseComponentId::new(phase.clone(), substance))
            .collect::<Vec<_>>();
        let mut phase_ranges = HashMap::new();
        phase_ranges.insert(phase.clone(), 0..components.len());
        Self {
            phases: vec![phase],
            components,
            phase_ranges,
        }
    }

    /// Projects unordered phase payloads into a deterministic layout.
    ///
    /// Phase names are sorted while the declared substance order inside each
    /// `SubsData` payload is preserved.
    ///
    /// HashMap is unordered, so this function ensures that the resulting SystemLayout is deterministic and consistent across runs.
    pub fn from_phase_map(phases: &HashMap<Option<String>, SubsData>) -> Self {
        let mut ordered_phases: Vec<PhaseId> = phases.keys().cloned().map(PhaseId::new).collect();
        ordered_phases.sort();

        let mut components = Vec::new();
        let mut phase_ranges = HashMap::with_capacity(ordered_phases.len());

        for phase in &ordered_phases {
            let start = components.len();
            if let Some(subsdata) = phases.get(&phase.0) {
                for substance in subsdata.substances.iter().cloned() {
                    components.push(PhaseComponentId::new(phase.clone(), substance));
                }
            }
            phase_ranges.insert(phase.clone(), start..components.len());
        }

        Self {
            phases: ordered_phases,
            components,
            phase_ranges,
        }
    }

    /// Number of solver components across all phases.
    pub fn component_count(&self) -> usize {
        self.components.len()
    }

    /// Contiguous component range belonging to `phase`.
    pub fn phase_component_range(&self, phase: &PhaseId) -> Option<&Range<usize>> {
        self.phase_ranges.get(phase)
    }

    /// Position of `phase` in the canonical phase order.
    pub fn phase_index(&self, phase: &PhaseId) -> Option<usize> {
        self.phases.iter().position(|candidate| candidate == phase)
    }

    /// Position of a fully-qualified component identity.
    pub fn component_index(&self, component: &PhaseComponentId) -> Option<usize> {
        self.components
            .iter()
            .position(|candidate| candidate == component)
    }

    /// Position of a component identified by phase and substance.
    pub fn component_index_by_parts(&self, phase: &PhaseId, substance: &str) -> Option<usize> {
        self.components
            .iter()
            .position(|candidate| &candidate.phase == phase && candidate.substance == substance)
    }

    /// Canonical ordered phases.
    pub fn phases(&self) -> &[PhaseId] {
        &self.phases
    }

    /// Canonical ordered components.
    pub fn components(&self) -> &[PhaseComponentId] {
        &self.components
    }

    /// Components owned by one phase, preserving their declared order.
    pub fn components_for_phase(&self, phase: &PhaseId) -> Option<&[PhaseComponentId]> {
        self.phase_component_range(phase)
            .map(|range| &self.components[range.clone()])
    }

    /// Human-readable label for the component at `index`.
    pub fn component_label(&self, index: usize) -> Option<String> {
        self.components
            .get(index)
            .map(|component| component.label())
    }

    /// Human-readable labels for all components in canonical order.
    pub fn component_labels(&self) -> Vec<String> {
        (0..self.components.len())
            .filter_map(|index| self.component_label(index))
            .collect()
    }
}

/// Returns phase payloads in the deterministic order used by `SystemLayout`.
pub fn ordered_phase_entries<'a>(
    phases: &'a HashMap<Option<String>, SubsData>,
) -> Vec<(PhaseId, &'a SubsData)> {
    let mut ordered: Vec<(PhaseId, &'a SubsData)> = phases
        .iter()
        .map(|(phase, subsdata)| (PhaseId(phase.clone()), subsdata))
        .collect();
    ordered.sort_by(|left, right| left.0.cmp(&right.0));
    ordered
}

/// Returns only the ordered phase keys for iteration without borrowing payloads.
pub fn ordered_phase_keys(phases: &HashMap<Option<String>, SubsData>) -> Vec<Option<String>> {
    ordered_phase_entries(phases)
        .into_iter()
        .map(|(phase, _)| phase.0)
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_phase_id_and_component_id_helpers() {
        let phase = PhaseId::new(Some("gas".to_string()));
        assert_eq!(phase.as_option(), &Some("gas".to_string()));
        let component = PhaseComponentId::new(phase.clone(), "O2");
        assert_eq!(component.label(), "gas::O2".to_string());
        assert_eq!(component.phase, phase);
        assert_eq!(component.substance, "O2".to_string());
    }

    #[test]
    fn test_system_layout_uses_typed_order_and_phase_ranges() {
        let mut phases = HashMap::new();
        let mut gas = SubsData::new();
        gas.substances = vec!["O2".to_string(), "H2".to_string()];
        let mut liquid = SubsData::new();
        liquid.substances = vec!["H2O".to_string()];
        phases.insert(Some("liquid".to_string()), liquid);
        phases.insert(Some("gas".to_string()), gas);

        let layout = SystemLayout::from_phase_map(&phases);
        assert_eq!(
            layout.phases(),
            &[
                PhaseId::new(Some("gas".to_string())),
                PhaseId::new(Some("liquid".to_string()))
            ]
        );
        assert_eq!(layout.component_count(), 3);
        assert_eq!(
            layout.component_labels(),
            vec![
                "gas::O2".to_string(),
                "gas::H2".to_string(),
                "liquid::H2O".to_string()
            ]
        );

        let gas_phase = PhaseId::new(Some("gas".to_string()));
        let gas_components = layout.components_for_phase(&gas_phase).unwrap();
        assert_eq!(gas_components.len(), 2);
        assert_eq!(gas_components[0].label(), "gas::O2".to_string());
        assert_eq!(gas_components[1].label(), "gas::H2".to_string());
    }

    #[test]
    fn test_single_phase_layout_is_canonical_and_ordered() {
        let layout = SystemLayout::from_single_phase(&vec!["CO2".to_string(), "H2O".to_string()]);
        assert_eq!(layout.phases(), &[PhaseId::new(None)]);
        assert_eq!(
            layout.component_labels(),
            vec!["CO2".to_string(), "H2O".to_string()]
        );
        let phase = PhaseId::new(None);
        let range = layout.phase_component_range(&phase).unwrap();
        assert_eq!(range.start, 0);
        assert_eq!(range.end, 2);

        let component = PhaseComponentId::new(phase.clone(), "H2O");
        assert_eq!(layout.component_index(&component), Some(1));
        assert_eq!(layout.component_index_by_parts(&phase, "CO2"), Some(0));
    }
}
