//! Typed construction and repository resolution for phase systems.
//!
//! This module is the boundary between user declarations and resolved
//! `SubsData` payloads. It does lookup/I-O and never evaluates or owns
//! thermodynamic cache state.

use std::collections::{HashMap, HashSet};
use std::fmt;

use crate::Thermodynamics::User_PhaseOrSolution::{CustomSubstance, PhaseOrSolution};
use crate::Thermodynamics::User_PhaseOrSolution2::OnePhase;
use crate::Thermodynamics::User_substances::{LibraryPriority, Phases, SubsData};
use crate::Thermodynamics::User_substances_error::{SubsDataError, SubsDataResult};
use crate::Thermodynamics::phase_layout::PhaseId;
use crate::Thermodynamics::thermo_lib_api::{ThermoData, ThermoLibraryError, ThermoRepository};
use std::sync::Arc;

use super::{PhasePhysicalState, PhaseSpec, ResolvedPhaseSystem};

/// Maps substances to their respective phases while maintaining order.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SubstancePhaseMapping {
    pub all_substances: Vec<String>,
    pub substance_to_phase: Vec<Option<String>>,
}

impl SubstancePhaseMapping {
    pub fn new(substances: Vec<String>, phases: Vec<Option<String>>) -> SubsDataResult<Self> {
        if substances.len() != phases.len() {
            return Err(SubsDataError::calculation_failed(
                "multiple",
                "substance-phase mapping",
                format!(
                    "substances length {} does not match phases length {}",
                    substances.len(),
                    phases.len()
                ),
            ));
        }
        Ok(Self {
            all_substances: substances,
            substance_to_phase: phases,
        })
    }

    pub fn get_phase_for_substance(&self, index: usize) -> Option<&Option<String>> {
        self.substance_to_phase.get(index)
    }
}

/// Compatibility input for callers which have not yet declared `PhaseSpec`s.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SubstancesContainer {
    SinglePhase(Vec<String>),
    MultiPhase(HashMap<String, Vec<String>>),
}

impl SubstancesContainer {
    /// Returns all unique substance names in deterministic order.
    pub fn get_all_substances(&self) -> Vec<String> {
        let mut substances = match self {
            Self::SinglePhase(substances) => substances.clone(),
            Self::MultiPhase(phase_substances) => {
                phase_substances.values().flatten().cloned().collect()
            }
        };
        substances.sort();
        substances.dedup();
        substances
    }

    /// Preserves duplicate components when one substance occurs in two phases.
    pub fn get_ordered_component_labels(&self) -> Vec<(Option<String>, String)> {
        match self {
            Self::SinglePhase(substances) => substances
                .iter()
                .cloned()
                .map(|substance| (None, substance))
                .collect(),
            Self::MultiPhase(phase_substances) => {
                let mut phase_names = phase_substances.keys().cloned().collect::<Vec<_>>();
                phase_names.sort();
                phase_names
                    .into_iter()
                    .flat_map(|phase_name| {
                        phase_substances
                            .get(&phase_name)
                            .into_iter()
                            .flatten()
                            .cloned()
                            .map(move |substance| (Some(phase_name.clone()), substance))
                    })
                    .collect()
            }
        }
    }

    pub fn get_ordered_substances(&self) -> Vec<String> {
        self.get_ordered_component_labels()
            .into_iter()
            .map(|(_, substance)| substance)
            .collect()
    }
}

/// Typed specification for building a thermodynamic system before resolution.
#[derive(Debug, Clone)]
pub struct SubstanceSystemSpec {
    phases: Vec<PhaseSpec>,
    library_priorities: Vec<String>,
    permitted_libraries: Vec<String>,
    explicit_search_instructions: Option<HashMap<String, String>>,
    search_in_nist: bool,
}

impl SubstanceSystemSpec {
    pub fn from_phases(phases: Vec<PhaseSpec>) -> Result<Self, SubstanceSystemFactoryError> {
        let spec = Self {
            phases,
            library_priorities: Vec::new(),
            permitted_libraries: Vec::new(),
            explicit_search_instructions: None,
            search_in_nist: false,
        };
        SubstanceSystemFactory::validate_spec(&spec)?;
        Ok(spec)
    }

    pub fn phases(&self) -> &[PhaseSpec] {
        &self.phases
    }

    pub fn builder(container: SubstancesContainer) -> SubstanceSystemSpecBuilder {
        SubstanceSystemSpecBuilder::new(container)
    }

    pub fn resolve(self) -> Result<CustomSubstance, SubstanceSystemFactoryError> {
        SubstanceSystemFactory::resolve_spec(self)
    }

    /// Resolves this specification against an explicitly supplied immutable
    /// catalog. All phase payloads share that same repository handle.
    pub fn resolve_with_repository(
        self,
        repository: Arc<ThermoRepository>,
    ) -> Result<CustomSubstance, SubstanceSystemFactoryError> {
        SubstanceSystemFactory::resolve_spec_with_repository(self, repository)
    }
}

/// Fluent compatibility builder for `SubstanceSystemSpec`.
#[derive(Debug, Clone)]
pub struct SubstanceSystemSpecBuilder {
    container: SubstancesContainer,
    phase_natures: Option<HashMap<String, Phases>>,
    library_priorities: Vec<String>,
    permitted_libraries: Vec<String>,
    explicit_search_instructions: Option<HashMap<String, String>>,
    search_in_nist: bool,
}

impl SubstanceSystemSpecBuilder {
    pub fn new(container: SubstancesContainer) -> Self {
        Self {
            container,
            phase_natures: None,
            library_priorities: Vec::new(),
            permitted_libraries: Vec::new(),
            explicit_search_instructions: None,
            search_in_nist: false,
        }
    }

    /// Supplies the physical state for every named phase in a multi-phase
    /// compatibility container. `MultiPhase` has no safe implicit state: a
    /// missing entry is rejected instead of silently becoming an ideal gas.
    /// `SinglePhase` remains the explicit ideal-gas convenience form.
    pub fn with_phase_natures(mut self, phase_natures: Option<HashMap<String, Phases>>) -> Self {
        self.phase_natures = phase_natures;
        self
    }

    pub fn with_library_priorities(mut self, library_priorities: Vec<String>) -> Self {
        self.library_priorities = library_priorities;
        self
    }

    pub fn with_permitted_libraries(mut self, permitted_libraries: Vec<String>) -> Self {
        self.permitted_libraries = permitted_libraries;
        self
    }

    pub fn with_explicit_search_instructions(
        mut self,
        explicit_search_instructions: Option<HashMap<String, String>>,
    ) -> Self {
        self.explicit_search_instructions = explicit_search_instructions;
        self
    }

    pub fn with_search_in_nist(mut self, search_in_nist: bool) -> Self {
        self.search_in_nist = search_in_nist;
        self
    }

    pub fn build(self) -> Result<SubstanceSystemSpec, SubstanceSystemFactoryError> {
        let phases = Self::legacy_container_to_phase_specs(self.container, self.phase_natures)?;
        let spec = SubstanceSystemSpec {
            phases,
            library_priorities: self.library_priorities,
            permitted_libraries: self.permitted_libraries,
            explicit_search_instructions: self.explicit_search_instructions,
            search_in_nist: self.search_in_nist,
        };
        SubstanceSystemFactory::validate_spec(&spec)?;
        Ok(spec)
    }

    fn legacy_container_to_phase_specs(
        container: SubstancesContainer,
        phase_natures: Option<HashMap<String, Phases>>,
    ) -> Result<Vec<PhaseSpec>, SubstanceSystemFactoryError> {
        match container {
            SubstancesContainer::SinglePhase(components) => {
                Ok(vec![PhaseSpec::ideal_gas(PhaseId::new(None), components)?])
            }
            SubstancesContainer::MultiPhase(phases) => {
                if phases.is_empty() {
                    return Err(SubstanceSystemFactoryError::InvalidSpecification {
                        field: "phases".to_string(),
                        reason: "a phase system must contain at least one phase".to_string(),
                    });
                }
                let natures = phase_natures.ok_or_else(|| {
                    SubstanceSystemFactoryError::InvalidSpecification {
                        field: "physical phase state".to_string(),
                        reason: "named multi-phase systems require an explicit physical state for every phase"
                            .to_string(),
                    }
                })?;
                for phase_name in phases.keys() {
                    if !natures.contains_key(phase_name) {
                        return Err(SubstanceSystemFactoryError::InvalidSpecification {
                            field: "physical phase state".to_string(),
                            reason: format!(
                                "missing physical state for declared phase '{}'",
                                phase_name
                            ),
                        });
                    }
                }
                for phase_name in natures.keys() {
                    if !phases.contains_key(phase_name) {
                        return Err(SubstanceSystemFactoryError::InvalidSpecification {
                            field: "physical phase state".to_string(),
                            reason: format!(
                                "physical state declared for unknown phase '{}'",
                                phase_name
                            ),
                        });
                    }
                }
                phases
                    .into_iter()
                    .map(|(name, components)| {
                        let id = PhaseId::new(Some(name.clone()));
                        let physical_state = natures.get(&name).copied().ok_or_else(|| {
                            SubstanceSystemFactoryError::InvalidSpecification {
                                field: "physical phase state".to_string(),
                                reason: format!(
                                    "missing physical state for declared phase '{}'",
                                    name
                                ),
                            }
                        })?;
                        match PhasePhysicalState::from(physical_state) {
                            PhasePhysicalState::Gas => PhaseSpec::ideal_gas(id, components),
                            physical_state => {
                                PhaseSpec::pure_condensed(id, components, physical_state)
                            }
                        }
                    })
                    .collect()
            }
        }
    }
}

/// Typed errors for phase-system specification and resolution.
#[derive(Debug)]
pub enum SubstanceSystemFactoryError {
    InvalidSpecification { field: String, reason: String },
    Repository(ThermoLibraryError),
    Resolution(SubsDataError),
}

impl fmt::Display for SubstanceSystemFactoryError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidSpecification { field, reason } => {
                write!(f, "invalid specification in {}: {}", field, reason)
            }
            Self::Repository(err) => write!(f, "thermodynamic repository error: {}", err),
            Self::Resolution(err) => write!(f, "{}", err),
        }
    }
}

impl std::error::Error for SubstanceSystemFactoryError {}

impl From<SubsDataError> for SubstanceSystemFactoryError {
    fn from(value: SubsDataError) -> Self {
        Self::Resolution(value)
    }
}

/// Resolves typed phase specifications into ready-to-evaluate facades.
pub struct SubstanceSystemFactory;

impl SubstanceSystemFactory {
    pub(crate) fn validate_spec(
        spec: &SubstanceSystemSpec,
    ) -> Result<(), SubstanceSystemFactoryError> {
        if spec.phases.is_empty() {
            return Err(SubstanceSystemFactoryError::InvalidSpecification {
                field: "phases".to_string(),
                reason: "a phase system must contain at least one phase".to_string(),
            });
        }
        let mut seen_ids = HashSet::with_capacity(spec.phases.len());
        for phase in &spec.phases {
            phase.validate()?;
            if !seen_ids.insert(phase.id().clone()) {
                return Err(SubstanceSystemFactoryError::InvalidSpecification {
                    field: "phase id".to_string(),
                    reason: format!("phase {:?} is declared more than once", phase.id()),
                });
            }
        }
        Ok(())
    }

    fn apply_lookup_policy(
        subs_data: &mut SubsData,
        library_priorities: &[String],
        permitted_libraries: &[String],
        explicit_search_instructions: Option<&HashMap<String, String>>,
    ) {
        subs_data.set_multiple_library_priorities(
            library_priorities.to_vec(),
            LibraryPriority::Priority,
        );
        subs_data.set_multiple_library_priorities(
            permitted_libraries.to_vec(),
            LibraryPriority::Permitted,
        );
        if let Some(instructions) = explicit_search_instructions {
            subs_data.set_explicit_search_instructions(instructions.clone());
        }
    }

    fn resolve_subs_data(
        subs_data: &mut SubsData,
        search_in_nist: bool,
    ) -> Result<(), SubstanceSystemFactoryError> {
        subs_data.search_substances()?;
        if search_in_nist {
            subs_data.if_not_found_go_NIST()?;
        }
        subs_data.parse_all_thermal_coeffs()?;
        Ok(())
    }

    /// Resolves validated declarations into one immutable, layout-aligned system.
    pub fn resolve_phase_system(
        spec: SubstanceSystemSpec,
    ) -> Result<ResolvedPhaseSystem, SubstanceSystemFactoryError> {
        let repository = ThermoData::try_default_repository()
            .map_err(SubstanceSystemFactoryError::Repository)?;
        Self::resolve_phase_system_with_repository(spec, repository)
    }

    /// Resolves validated declarations against one injected immutable catalog.
    ///
    /// Every resulting `SubsData` owns query-local calculators and search state,
    /// but borrows the same read-only library payload through `Arc`.
    pub fn resolve_phase_system_with_repository(
        spec: SubstanceSystemSpec,
        repository: Arc<ThermoRepository>,
    ) -> Result<ResolvedPhaseSystem, SubstanceSystemFactoryError> {
        Self::validate_spec(&spec)?;
        let SubstanceSystemSpec {
            phases,
            library_priorities,
            permitted_libraries,
            explicit_search_instructions,
            search_in_nist,
        } = spec;
        let mut resolved_data = HashMap::with_capacity(phases.len());
        for phase in &phases {
            let phase_key = phase.id().as_option().clone();
            if phases.len() > 1 && phase_key.is_none() {
                return Err(SubstanceSystemFactoryError::InvalidSpecification {
                    field: "phase id".to_string(),
                    reason: "multi-phase systems require named phases".to_string(),
                });
            }
            let mut phase_data = SubsData::from_thermo_repository(Arc::clone(&repository));
            phase_data.substances = phase.components().to_vec();
            Self::apply_lookup_policy(
                &mut phase_data,
                &library_priorities,
                &permitted_libraries,
                explicit_search_instructions.as_ref(),
            );
            // A phase model carries physical context, but state-specific
            // library lookup remains opt-in. The compatibility projection is
            // retained for NIST fallback without restricting local records.
            phase_data.map_of_phases = phase.legacy_component_phase_map();
            Self::resolve_subs_data(&mut phase_data, search_in_nist)?;
            resolved_data.insert(phase_key, phase_data);
        }
        ResolvedPhaseSystem::new_with_nist_fallback_policy(phases, resolved_data, search_in_nist)
    }

    pub fn resolve_spec(
        spec: SubstanceSystemSpec,
    ) -> Result<CustomSubstance, SubstanceSystemFactoryError> {
        let repository = ThermoData::try_default_repository()
            .map_err(SubstanceSystemFactoryError::Repository)?;
        Self::resolve_spec_with_repository(spec, repository)
    }

    /// Resolves a facade against one explicit repository handle.
    pub fn resolve_spec_with_repository(
        spec: SubstanceSystemSpec,
        repository: Arc<ThermoRepository>,
    ) -> Result<CustomSubstance, SubstanceSystemFactoryError> {
        let resolved = Self::resolve_phase_system_with_repository(spec, repository)?;
        let is_single_phase = resolved.phase_specs().len() == 1
            && resolved.phase_specs()[0].id().as_option().is_none();
        if is_single_phase {
            let mut one_phase = OnePhase::new();
            one_phase.install_resolved_system(resolved)?;
            Ok(CustomSubstance::OnePhase(one_phase))
        } else {
            let mut phase_system = PhaseOrSolution::new();
            phase_system.install_resolved_system(resolved);
            Ok(CustomSubstance::PhaseOrSolution(phase_system))
        }
    }

    /// Compatibility adapter returning the legacy string error type.
    pub fn create_system(
        container: SubstancesContainer,
        physical_phase_natures: Option<HashMap<String, Phases>>,
        library_priorities: Vec<String>,
        permitted_libraries: Vec<String>,
        explicit_search_instructions: Option<HashMap<String, String>>,
        search_in_nist: bool,
    ) -> Result<CustomSubstance, String> {
        let spec = SubstanceSystemSpecBuilder::new(container)
            .with_phase_natures(physical_phase_natures)
            .with_library_priorities(library_priorities)
            .with_permitted_libraries(permitted_libraries)
            .with_explicit_search_instructions(explicit_search_instructions)
            .with_search_in_nist(search_in_nist)
            .build()
            .map_err(|error| error.to_string())?;
        Self::resolve_spec(spec).map_err(|error| error.to_string())
    }
}
