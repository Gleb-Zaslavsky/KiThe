//! Deterministic candidate selection for element-defined equilibrium systems.
//!
//! Element search is deliberately kept separate from solving.  A catalog can
//! tell us that a name occurs under an element, but that fact alone does not
//! decide which library record, physical state, or temperature interval is
//! suitable for a production equilibrium problem.  This module makes those
//! decisions explicit and returns an auditable report instead of mutating
//! `SubsData` or the repository.

use crate::Thermodynamics::User_PhaseOrSolution::{
    PhaseModel, PhaseSpec, SubstanceSystemFactoryError, SubstanceSystemSpec,
};
use crate::Thermodynamics::phase_layout::PhaseId;
use crate::Thermodynamics::physical_state::PhysicalState;
use crate::Thermodynamics::thermo_lib_api::{ElementSearchMode, ThermoData, ThermoRepository};
use serde_json::Value;
use std::collections::{HashMap, HashSet};
use std::fmt;
use std::sync::Arc;

/// A requested temperature interval for candidate screening.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CandidateTemperatureRange {
    lower: f64,
    upper: f64,
}

impl CandidateTemperatureRange {
    /// Creates a finite, positive, non-empty temperature interval.
    pub fn new(lower: f64, upper: f64) -> Result<Self, CandidateSelectionError> {
        if !lower.is_finite() || !upper.is_finite() || lower <= 0.0 || lower > upper {
            return Err(CandidateSelectionError::InvalidPolicy(
                "temperature interval must be finite, positive, and lower <= upper".into(),
            ));
        }
        Ok(Self { lower, upper })
    }

    pub fn lower(self) -> f64 {
        self.lower
    }

    pub fn upper(self) -> f64 {
        self.upper
    }
}

/// How much temperature metadata was available in a catalog record.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CandidateTemperatureSupport {
    /// At least one stored coefficient interval covers the requested range.
    Supported,
    /// Stored intervals exist, but none covers the requested range.
    Unsupported,
    /// The library schema does not expose a recognizable interval field.
    Unknown,
}

/// Policy controlling conversion from an element query to equilibrium
/// candidates.
#[derive(Debug, Clone, PartialEq)]
pub struct EquilibriumCandidatePolicy {
    element_mode: ElementSearchMode,
    library_preference: Vec<String>,
    physical_states: Option<Vec<PhysicalState>>,
    temperature_range: Option<CandidateTemperatureRange>,
    max_candidates: Option<usize>,
}

impl Default for EquilibriumCandidatePolicy {
    fn default() -> Self {
        Self {
            element_mode: ElementSearchMode::SubsetOf,
            library_preference: Vec::new(),
            physical_states: None,
            temperature_range: None,
            max_candidates: None,
        }
    }
}

impl EquilibriumCandidatePolicy {
    pub fn new(element_mode: ElementSearchMode) -> Self {
        Self {
            element_mode,
            ..Self::default()
        }
    }

    pub fn with_library_preference(mut self, libraries: Vec<String>) -> Self {
        self.library_preference = libraries;
        self
    }

    pub fn with_physical_states(mut self, states: Vec<PhysicalState>) -> Self {
        self.physical_states = Some(states);
        self
    }

    pub fn with_temperature_range(
        mut self,
        lower: f64,
        upper: f64,
    ) -> Result<Self, CandidateSelectionError> {
        self.temperature_range = Some(CandidateTemperatureRange::new(lower, upper)?);
        Ok(self)
    }

    pub fn with_max_candidates(mut self, limit: usize) -> Result<Self, CandidateSelectionError> {
        if limit == 0 {
            return Err(CandidateSelectionError::InvalidPolicy(
                "max_candidates must be greater than zero".into(),
            ));
        }
        self.max_candidates = Some(limit);
        Ok(self)
    }

    pub fn element_mode(&self) -> ElementSearchMode {
        self.element_mode
    }

    pub fn library_preference(&self) -> &[String] {
        &self.library_preference
    }

    pub fn physical_states(&self) -> Option<&[PhysicalState]> {
        self.physical_states.as_deref()
    }

    pub fn temperature_range(&self) -> Option<CandidateTemperatureRange> {
        self.temperature_range
    }

    pub fn max_candidates(&self) -> Option<usize> {
        self.max_candidates
    }
}

/// Why a catalog record was not selected.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CandidateRejectionReason {
    ElementSetMismatch,
    LibraryNotPreferred,
    PhysicalStateMismatch,
    TemperatureUnsupported,
    MissingRecord,
    CandidateLimit,
}

/// One rejected record retained for diagnostics.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CandidateRejection {
    substance: String,
    library: String,
    reason: CandidateRejectionReason,
}

impl CandidateRejection {
    pub fn substance(&self) -> &str {
        &self.substance
    }

    pub fn library(&self) -> &str {
        &self.library
    }

    pub fn reason(&self) -> &CandidateRejectionReason {
        &self.reason
    }
}

/// A record that passed the candidate policy and can be used to construct a
/// `PhaseSpec`/lookup request.
#[derive(Debug, Clone, PartialEq)]
pub struct EquilibriumCandidate {
    substance: String,
    library: String,
    record_key: String,
    physical_state: Option<PhysicalState>,
    elements: Vec<String>,
    temperature_support: CandidateTemperatureSupport,
    library_rank: usize,
}

impl EquilibriumCandidate {
    pub fn substance(&self) -> &str {
        &self.substance
    }

    pub fn library(&self) -> &str {
        &self.library
    }

    pub fn record_key(&self) -> &str {
        &self.record_key
    }

    pub fn physical_state(&self) -> Option<PhysicalState> {
        self.physical_state
    }

    pub fn elements(&self) -> &[String] {
        &self.elements
    }

    pub fn temperature_support(&self) -> CandidateTemperatureSupport {
        self.temperature_support
    }

    pub fn library_rank(&self) -> usize {
        self.library_rank
    }
}

/// Complete immutable evidence from one candidate-selection transaction.
#[derive(Debug, Clone, PartialEq)]
pub struct EquilibriumCandidateSelectionReport {
    requested_elements: Vec<String>,
    policy: EquilibriumCandidatePolicy,
    selected: Vec<EquilibriumCandidate>,
    rejected: Vec<CandidateRejection>,
}

impl EquilibriumCandidateSelectionReport {
    pub fn requested_elements(&self) -> &[String] {
        &self.requested_elements
    }

    pub fn policy(&self) -> &EquilibriumCandidatePolicy {
        &self.policy
    }

    pub fn selected(&self) -> &[EquilibriumCandidate] {
        &self.selected
    }

    pub fn rejected(&self) -> &[CandidateRejection] {
        &self.rejected
    }
}

/// Explicit assignment of selected records to one physical phase.
///
/// Candidate selection answers "which records satisfy the catalog policy?";
/// it must not answer "which activity model should this phase use?".  This
/// value keeps that decision visible at the production boundary.  Components
/// are identified by their exact `record_key`, not by a broadened base name,
/// so condensed polymorphs such as `Fe(a)` and `Fe(c)` cannot be silently
/// substituted for one another.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct EquilibriumCandidatePhaseAssignment {
    phase_id: PhaseId,
    physical_state: PhysicalState,
    model: PhaseModel,
    record_keys: Vec<String>,
}

impl EquilibriumCandidatePhaseAssignment {
    /// Creates one explicit phase assignment.
    pub fn new(
        phase_id: PhaseId,
        physical_state: PhysicalState,
        model: PhaseModel,
        record_keys: Vec<String>,
    ) -> Self {
        Self {
            phase_id,
            physical_state,
            model,
            record_keys,
        }
    }

    /// Convenience constructor for the currently supported ideal-gas model.
    pub fn ideal_gas(phase_id: PhaseId, record_keys: Vec<String>) -> Self {
        Self::new(
            phase_id,
            PhysicalState::Gas,
            PhaseModel::IdealGas,
            record_keys,
        )
    }

    /// Convenience constructor for a one-component pure condensed phase.
    pub fn pure_condensed(
        phase_id: PhaseId,
        physical_state: PhysicalState,
        record_keys: Vec<String>,
    ) -> Self {
        Self::new(
            phase_id,
            physical_state,
            PhaseModel::PureCondensed,
            record_keys,
        )
    }

    /// Convenience constructor for an explicitly ideal condensed solution.
    pub fn ideal_solution(
        phase_id: PhaseId,
        physical_state: PhysicalState,
        record_keys: Vec<String>,
    ) -> Self {
        Self::new(
            phase_id,
            physical_state,
            PhaseModel::IdealSolution,
            record_keys,
        )
    }

    pub fn phase_id(&self) -> &PhaseId {
        &self.phase_id
    }

    pub fn physical_state(&self) -> PhysicalState {
        self.physical_state
    }

    pub fn model(&self) -> PhaseModel {
        self.model
    }

    pub fn record_keys(&self) -> &[String] {
        &self.record_keys
    }
}

/// Physical phase plan that turns selected catalog records into a solver
/// specification.
///
/// The plan is deliberately separate from [`EquilibriumCandidatePolicy`].
/// The former is a user/model decision; the latter is a repository search
/// decision.  Keeping them separate prevents a broad element search from
/// accidentally creating an ideal-gas phase or a solution model.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct EquilibriumCandidatePhasePlan {
    assignments: Vec<EquilibriumCandidatePhaseAssignment>,
}

impl EquilibriumCandidatePhasePlan {
    pub fn new(assignments: Vec<EquilibriumCandidatePhaseAssignment>) -> Self {
        Self { assignments }
    }

    pub fn assignments(&self) -> &[EquilibriumCandidatePhaseAssignment] {
        &self.assignments
    }

    /// Builds a typed specification and pins each selected record to the
    /// library reported by the selection transaction.
    pub fn build_spec(
        &self,
        selection: &EquilibriumCandidateSelectionReport,
    ) -> Result<SubstanceSystemSpec, SubstanceSystemFactoryError> {
        let selected_by_key: HashMap<_, _> = selection
            .selected()
            .iter()
            .map(|candidate| (candidate.record_key().to_string(), candidate))
            .collect();
        if selected_by_key.len() != selection.selected().len() {
            return Err(invalid_plan(
                "candidate selection contains duplicate record keys",
            ));
        }

        let mut assigned = HashSet::new();
        let mut phase_specs = Vec::with_capacity(self.assignments.len());
        let mut explicit_libraries = HashMap::with_capacity(selected_by_key.len());
        let mut library_order = selection
            .selected()
            .iter()
            .map(|candidate| (candidate.library_rank(), candidate.library().to_string()))
            .collect::<Vec<_>>();
        library_order.sort_by_key(|(rank, _)| *rank);
        library_order.dedup_by(|left, right| left.1 == right.1);

        for assignment in &self.assignments {
            let mut components = Vec::with_capacity(assignment.record_keys.len());
            for record_key in &assignment.record_keys {
                let Some(candidate) = selected_by_key.get(record_key) else {
                    return Err(invalid_plan(format!(
                        "phase {:?} references unselected record '{}'",
                        assignment.phase_id.as_option(),
                        record_key
                    )));
                };
                if !assigned.insert(record_key.clone()) {
                    return Err(invalid_plan(format!(
                        "record '{}' is assigned to more than one phase",
                        record_key
                    )));
                }
                components.push(record_key.clone());
                explicit_libraries.insert(record_key.clone(), candidate.library().to_string());
            }
            let phase = PhaseSpec::new(
                assignment.phase_id.clone(),
                components,
                assignment.physical_state,
                assignment.model,
            )
            .map_err(|error| invalid_plan(error.to_string()))?;
            phase_specs.push(phase);
        }

        if assigned.len() != selected_by_key.len() {
            let missing = selected_by_key
                .keys()
                .filter(|record_key| !assigned.contains(*record_key))
                .cloned()
                .collect::<Vec<_>>();
            return Err(invalid_plan(format!(
                "selected records were not assigned to a phase: {missing:?}"
            )));
        }

        SubstanceSystemSpec::from_phases(phase_specs).map(|spec| {
            spec.with_lookup_policy(
                library_order
                    .iter()
                    .map(|(_, library)| library.clone())
                    .collect(),
                library_order
                    .iter()
                    .map(|(_, library)| library.clone())
                    .collect(),
                Some(explicit_libraries),
                false,
            )
        })
    }
}

fn invalid_plan(reason: impl Into<String>) -> SubstanceSystemFactoryError {
    SubstanceSystemFactoryError::InvalidSpecification {
        field: "candidate phase plan".to_string(),
        reason: reason.into(),
    }
}

/// Typed errors raised before any equilibrium problem is constructed.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CandidateSelectionError {
    InvalidPolicy(String),
    EmptyElementSet,
    UnknownLibrary(String),
}

impl fmt::Display for CandidateSelectionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidPolicy(message) => write!(f, "invalid candidate policy: {message}"),
            Self::EmptyElementSet => write!(f, "candidate selection requires at least one element"),
            Self::UnknownLibrary(library) => {
                write!(f, "unknown thermochemical library '{library}'")
            }
        }
    }
}

impl std::error::Error for CandidateSelectionError {}

/// Selects real thermochemical records without mutating the repository.
pub struct EquilibriumCandidateSelector {
    repository: Arc<ThermoRepository>,
}

impl EquilibriumCandidateSelector {
    pub fn new(repository: Arc<ThermoRepository>) -> Self {
        Self { repository }
    }

    pub fn repository(&self) -> &Arc<ThermoRepository> {
        &self.repository
    }

    /// Performs one deterministic, auditable selection transaction.
    pub fn select(
        &self,
        elements: &[String],
        policy: EquilibriumCandidatePolicy,
    ) -> Result<EquilibriumCandidateSelectionReport, CandidateSelectionError> {
        let requested_elements = normalize_elements(elements)?;
        let libraries = ordered_libraries(&self.repository, &policy)?;
        let library_ranks: HashMap<_, _> = libraries
            .iter()
            .enumerate()
            .map(|(rank, library)| (library.clone(), rank))
            .collect();
        let requested_set: HashSet<_> = requested_elements.iter().cloned().collect();

        let mut element_sets: HashMap<(String, String), HashSet<String>> = HashMap::new();
        for (element, pairs) in self.repository.ElementsData.iter() {
            for pair in pairs {
                if pair.len() < 2 {
                    continue;
                }
                let library = ThermoData::canonical_library_name(&pair[1]);
                if library_ranks.contains_key(&library) {
                    element_sets
                        .entry((library, pair[0].clone()))
                        .or_default()
                        .insert(element.clone());
                }
            }
        }

        let mut keys: Vec<_> = element_sets.keys().cloned().collect();
        keys.sort_by(|left, right| {
            library_ranks[&left.0]
                .cmp(&library_ranks[&right.0])
                .then_with(|| left.1.cmp(&right.1))
        });

        let mut selected = Vec::new();
        let mut rejected = Vec::new();
        let mut selected_substances = HashSet::new();
        for (library, substance) in keys {
            let elements = &element_sets[&(library.clone(), substance.clone())];
            if !element_match(&requested_set, elements, policy.element_mode()) {
                rejected.push(rejection(
                    &substance,
                    &library,
                    CandidateRejectionReason::ElementSetMismatch,
                ));
                continue;
            }
            let Some(records) = self.repository.LibThermoData.get(&library) else {
                rejected.push(rejection(
                    &substance,
                    &library,
                    CandidateRejectionReason::MissingRecord,
                ));
                continue;
            };
            let Some(data) = records.get(&substance) else {
                rejected.push(rejection(
                    &substance,
                    &library,
                    CandidateRejectionReason::MissingRecord,
                ));
                continue;
            };
            let query = crate::Thermodynamics::physical_state::ThermoRecordQuery::new(&substance);
            let Some(record) = self
                .repository
                .resolve_thermo_record(&library, &query)
                .map_err(|_| CandidateSelectionError::UnknownLibrary(library.clone()))?
            else {
                rejected.push(rejection(
                    &substance,
                    &library,
                    CandidateRejectionReason::MissingRecord,
                ));
                continue;
            };
            if !physical_state_match(record.physical_state, policy.physical_states()) {
                rejected.push(rejection(
                    &substance,
                    &library,
                    CandidateRejectionReason::PhysicalStateMismatch,
                ));
                continue;
            }
            let support = temperature_support(data, policy.temperature_range());
            if support == CandidateTemperatureSupport::Unsupported {
                rejected.push(rejection(
                    &substance,
                    &library,
                    CandidateRejectionReason::TemperatureUnsupported,
                ));
                continue;
            }
            if policy
                .max_candidates()
                .is_some_and(|limit| selected.len() >= limit)
            {
                rejected.push(rejection(
                    &substance,
                    &library,
                    CandidateRejectionReason::CandidateLimit,
                ));
                continue;
            }
            if selected_substances.contains(&substance) {
                rejected.push(rejection(
                    &substance,
                    &library,
                    CandidateRejectionReason::LibraryNotPreferred,
                ));
                continue;
            }
            let library_rank = library_ranks[&library];
            let mut elements = elements.iter().cloned().collect::<Vec<_>>();
            elements.sort();
            selected.push(EquilibriumCandidate {
                substance: substance.clone(),
                library,
                record_key: record.record_key,
                physical_state: record.physical_state,
                elements,
                temperature_support: support,
                library_rank,
            });
            selected_substances.insert(substance);
        }

        Ok(EquilibriumCandidateSelectionReport {
            requested_elements,
            policy,
            selected,
            rejected,
        })
    }
}

fn normalize_elements(elements: &[String]) -> Result<Vec<String>, CandidateSelectionError> {
    let mut normalized: Vec<_> = elements
        .iter()
        .map(|element| element.trim().to_string())
        .filter(|element| !element.is_empty())
        .collect();
    normalized.sort();
    normalized.dedup();
    if normalized.is_empty() {
        return Err(CandidateSelectionError::EmptyElementSet);
    }
    Ok(normalized)
}

fn ordered_libraries(
    repository: &ThermoRepository,
    policy: &EquilibriumCandidatePolicy,
) -> Result<Vec<String>, CandidateSelectionError> {
    let source = if policy.library_preference().is_empty() {
        repository.thermo_libs.as_ref().clone()
    } else {
        policy.library_preference().to_vec()
    };
    let mut libraries = Vec::new();
    for library in source {
        let canonical = ThermoData::canonical_library_name(&library);
        if !repository.LibThermoData.contains_key(&canonical) {
            return Err(CandidateSelectionError::UnknownLibrary(library));
        }
        if !libraries.contains(&canonical) {
            libraries.push(canonical);
        }
    }
    Ok(libraries)
}

fn element_match(
    requested: &HashSet<String>,
    candidate: &HashSet<String>,
    mode: ElementSearchMode,
) -> bool {
    match mode {
        ElementSearchMode::AnyRequested => !requested.is_disjoint(candidate),
        ElementSearchMode::SubsetOf => candidate.is_subset(requested),
        ElementSearchMode::ExactSet => candidate == requested,
    }
}

fn physical_state_match(
    observed: Option<PhysicalState>,
    requested: Option<&[PhysicalState]>,
) -> bool {
    requested.is_none_or(|states| {
        observed.is_some_and(|observed| states.iter().any(|state| state.accepts(observed)))
    })
}

fn rejection(
    substance: &str,
    library: &str,
    reason: CandidateRejectionReason,
) -> CandidateRejection {
    CandidateRejection {
        substance: substance.to_string(),
        library: library.to_string(),
        reason,
    }
}

fn temperature_support(
    data: &Value,
    requested: Option<CandidateTemperatureRange>,
) -> CandidateTemperatureSupport {
    let Some(requested) = requested else {
        return CandidateTemperatureSupport::Unknown;
    };
    let ranges = temperature_ranges(data);
    if ranges.is_empty() {
        return CandidateTemperatureSupport::Unknown;
    }
    if ranges
        .iter()
        .any(|(lower, upper)| *lower <= requested.lower() && *upper >= requested.upper())
    {
        CandidateTemperatureSupport::Supported
    } else {
        CandidateTemperatureSupport::Unsupported
    }
}

/// Extracts only interval-like fields, avoiding accidental interpretation of
/// arbitrary coefficient arrays as temperature ranges.
fn temperature_ranges(value: &Value) -> Vec<(f64, f64)> {
    let mut ranges = Vec::new();
    if let Value::Object(object) = value {
        for (key, value) in object {
            let key_lower = key.to_ascii_lowercase();
            let hinted =
                key == "T" || key_lower.contains("temperature") || key_lower.contains("range");
            if hinted {
                collect_ranges(value, &mut ranges);
            }
        }
    }
    ranges.sort_by(|left, right| left.partial_cmp(right).unwrap_or(std::cmp::Ordering::Equal));
    ranges.dedup();
    ranges
}

fn collect_ranges(value: &Value, ranges: &mut Vec<(f64, f64)>) {
    match value {
        Value::Array(values) if values.len() == 2 => {
            if let (Some(lower), Some(upper)) = (values[0].as_f64(), values[1].as_f64()) {
                if lower.is_finite() && upper.is_finite() && lower > 0.0 && lower <= upper {
                    ranges.push((lower, upper));
                    return;
                }
            }
            for value in values {
                collect_ranges(value, ranges);
            }
        }
        Value::Array(values) => {
            for value in values {
                collect_ranges(value, ranges);
            }
        }
        _ => {}
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::json;
    use std::collections::HashMap;

    #[test]
    fn policy_validates_temperature_and_candidate_limits() {
        assert!(CandidateTemperatureRange::new(300.0, 1000.0).is_ok());
        assert!(CandidateTemperatureRange::new(0.0, 1000.0).is_err());
        assert!(
            EquilibriumCandidatePolicy::default()
                .with_max_candidates(0)
                .is_err()
        );
    }

    #[test]
    fn temperature_metadata_is_conservative() {
        let data = serde_json::json!({"T": [[200.0, 1000.0], [1000.0, 6000.0]], "Cp": [1.0]});
        assert_eq!(
            temperature_support(
                &data,
                Some(CandidateTemperatureRange::new(300.0, 500.0).unwrap())
            ),
            CandidateTemperatureSupport::Supported
        );
        assert_eq!(
            temperature_support(
                &data,
                Some(CandidateTemperatureRange::new(7000.0, 8000.0).unwrap())
            ),
            CandidateTemperatureSupport::Unsupported
        );
        assert_eq!(
            temperature_support(
                &serde_json::json!({"Cp": [1.0]}),
                Some(CandidateTemperatureRange::new(300.0, 500.0).unwrap())
            ),
            CandidateTemperatureSupport::Unknown
        );
    }

    #[test]
    fn selector_is_exact_state_aware_and_does_not_mutate_repository() {
        let repository = Arc::new(ThermoRepository::from_parts(
            vec![
                ("NASA_gas".into(), "CO".into()),
                ("NASA_gas".into(), "CO2".into()),
                ("NASA_cond".into(), "CO2(s)".into()),
            ],
            HashMap::from([
                (
                    "NASA_gas".into(),
                    HashMap::from([
                        (
                            "CO".into(),
                            json!({"T": [[200.0, 6000.0]], "model": "NASA"}),
                        ),
                        (
                            "CO2".into(),
                            json!({"T": [[200.0, 6000.0]], "model": "NASA"}),
                        ),
                    ]),
                ),
                (
                    "NASA_cond".into(),
                    HashMap::from([(
                        "CO2(s)".into(),
                        json!({"T": [[200.0, 1000.0]], "model": "NASA"}),
                    )]),
                ),
            ]),
            HashMap::from([
                (
                    "C".into(),
                    vec![
                        vec!["CO".into(), "NASA_gas".into()],
                        vec!["CO2".into(), "NASA_gas".into()],
                        vec!["CO2(s)".into(), "NASA_cond".into()],
                    ],
                ),
                (
                    "O".into(),
                    vec![
                        vec!["CO".into(), "NASA_gas".into()],
                        vec!["CO2".into(), "NASA_gas".into()],
                        vec!["CO2(s)".into(), "NASA_cond".into()],
                    ],
                ),
            ]),
            vec!["NASA_gas".into(), "NASA_cond".into()],
            HashMap::new(),
            HashMap::new(),
            vec!["NASA_gas".into(), "NASA_cond".into()],
            Vec::new(),
        ));
        let before = repository.ElementsData.as_ref().clone();
        let selector = EquilibriumCandidateSelector::new(Arc::clone(&repository));
        let policy = EquilibriumCandidatePolicy::new(ElementSearchMode::ExactSet)
            .with_library_preference(vec!["NASA_cond".into(), "NASA_gas".into()])
            .with_physical_states(vec![PhysicalState::Solid])
            .with_temperature_range(300.0, 900.0)
            .unwrap();
        let report = selector
            .select(&["O".into(), "C".into(), "C".into()], policy)
            .unwrap();

        assert_eq!(report.selected().len(), 1);
        assert_eq!(report.selected()[0].substance(), "CO2(s)");
        assert_eq!(report.selected()[0].library(), "NASA_cond");
        assert_eq!(
            report.selected()[0].temperature_support(),
            CandidateTemperatureSupport::Supported
        );
        assert_eq!(
            report.selected()[0].elements(),
            &["C".to_string(), "O".to_string()]
        );
        assert_eq!(repository.ElementsData.as_ref(), &before);
    }

    #[test]
    fn phase_plan_preserves_exact_record_and_library_provenance() {
        let repository = Arc::new(ThermoRepository::from_parts(
            vec![("NASA_cond".into(), "CO2(s)".into())],
            HashMap::from([(
                "NASA_cond".into(),
                HashMap::from([("CO2(s)".into(), json!({"T": [[200.0, 1000.0]]}))]),
            )]),
            HashMap::from([
                (
                    "C".into(),
                    vec![vec!["CO2(s)".to_string(), "NASA_cond".to_string()]],
                ),
                (
                    "O".into(),
                    vec![vec!["CO2(s)".to_string(), "NASA_cond".to_string()]],
                ),
            ]),
            vec!["NASA_cond".into()],
            HashMap::new(),
            HashMap::new(),
            vec!["NASA_cond".into()],
            Vec::new(),
        ));
        let selection = EquilibriumCandidateSelector::new(Arc::clone(&repository))
            .select(
                &["C".to_string(), "O".to_string()],
                EquilibriumCandidatePolicy::new(ElementSearchMode::ExactSet)
                    .with_physical_states(vec![PhysicalState::Solid]),
            )
            .unwrap();
        assert_eq!(selection.selected().len(), 1);

        let plan = EquilibriumCandidatePhasePlan::new(vec![
            EquilibriumCandidatePhaseAssignment::pure_condensed(
                PhaseId::new(Some("solid".to_string())),
                PhysicalState::Solid,
                vec!["CO2(s)".into()],
            ),
        ]);
        let spec = plan.build_spec(&selection).unwrap();

        assert_eq!(spec.phases()[0].components(), &["CO2(s)".to_string()]);
        assert_eq!(spec.library_priorities(), &["NASA_cond".to_string()]);
        assert_eq!(spec.permitted_libraries(), &["NASA_cond".to_string()]);
        assert_eq!(
            spec.explicit_search_instructions()
                .and_then(|instructions| instructions.get("CO2(s)")),
            Some(&"NASA_cond".to_string())
        );
        assert!(!spec.search_in_nist());
    }

    #[test]
    fn phase_plan_rejects_unassigned_duplicate_and_unknown_records() {
        let selection = EquilibriumCandidateSelectionReport {
            requested_elements: vec!["C".into()],
            policy: EquilibriumCandidatePolicy::default(),
            selected: vec![EquilibriumCandidate {
                substance: "C(gr)".into(),
                library: "NASA_cond".into(),
                record_key: "C(gr)".into(),
                physical_state: Some(PhysicalState::Solid),
                elements: vec!["C".into()],
                temperature_support: CandidateTemperatureSupport::Unknown,
                library_rank: 0,
            }],
            rejected: Vec::new(),
        };

        let omitted = EquilibriumCandidatePhasePlan::new(Vec::new());
        assert!(omitted.build_spec(&selection).is_err());

        let unknown = EquilibriumCandidatePhasePlan::new(vec![
            EquilibriumCandidatePhaseAssignment::pure_condensed(
                PhaseId::new(Some("solid".to_string())),
                PhysicalState::Solid,
                vec!["missing".into()],
            ),
        ]);
        assert!(unknown.build_spec(&selection).is_err());

        let duplicate = EquilibriumCandidatePhasePlan::new(vec![
            EquilibriumCandidatePhaseAssignment::pure_condensed(
                PhaseId::new(Some("first".to_string())),
                PhysicalState::Solid,
                vec!["C(gr)".into()],
            ),
            EquilibriumCandidatePhaseAssignment::pure_condensed(
                PhaseId::new(Some("second".to_string())),
                PhysicalState::Solid,
                vec!["C(gr)".into()],
            ),
        ]);
        assert!(duplicate.build_spec(&selection).is_err());
    }
}
