//! Exportable immutable evidence for reproducing an accepted equilibrium run.
//!
//! A solver result intentionally does not retain a mutable request or a live
//! repository handle. This module therefore captures the exact resolved record
//! identities and the *effective* numerical policy beside the accepted result.
//! It does not claim to hash arbitrary JSON payloads: a reviewed data-release
//! label is optional until the repository owns a versioned payload manifest.

use serde::{Deserialize, Serialize};
use std::fmt;

use crate::Thermodynamics::ChemEquilibrium::equilibrium_candidate_selection::EquilibriumCandidateSelectionReport;
use crate::Thermodynamics::ChemEquilibrium::phase_equilibrium_problem::PhaseEquilibriumInputKind;
use crate::Thermodynamics::ChemEquilibrium::phase_equilibrium_workflow::{
    EquilibriumSolveOptions, EquilibriumSolveOptionsSnapshot, ResolvedPhaseEquilibriumOutcome,
};
use crate::Thermodynamics::thermo_lib_api::ThermoCatalogConsistencyReport;

/// Current JSON-compatible capsule schema.
pub const EQUILIBRIUM_REPRODUCIBILITY_SCHEMA_VERSION: u32 = 2;

/// Phase-stability mathematics represented by a reproducibility capsule.
///
/// The tag prevents a historical `driving_force` artifact from being read as
/// evidence for the canonical constrained TPD workflow.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum PhaseStabilitySemantics {
    CanonicalTpdV1,
}

/// Stable symbolic representation of the physical input origin.
///
/// This is intentionally a capsule-local type so the JSON contract does not
/// depend on the internal bridge enum layout or derives.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum EquilibriumReproducibilityInputKind {
    ExplicitComposition,
    ElementInventory,
}

/// Typed failure while loading a reproducibility artifact.
#[derive(Debug)]
pub enum ReproducibilityCapsuleError {
    Json(serde_json::Error),
    UnsupportedSchema { found: Option<u64>, expected: u32 },
    ObsoleteStabilityField { field: String },
}

impl fmt::Display for ReproducibilityCapsuleError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Json(error) => write!(formatter, "invalid reproducibility JSON: {error}"),
            Self::UnsupportedSchema { found, expected } => write!(
                formatter,
                "unsupported equilibrium reproducibility schema {:?}; expected {expected}",
                found
            ),
            Self::ObsoleteStabilityField { field } => write!(
                formatter,
                "reproducibility capsule contains obsolete phase-stability field '{field}'"
            ),
        }
    }
}

impl std::error::Error for ReproducibilityCapsuleError {}

/// Immutable component-level identity of one selected thermochemistry record.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct EquilibriumRecordIdentity {
    /// Collision-free `phase::substance` component label.
    pub component: String,
    /// Semantic phase name.
    pub phase: String,
    /// Bare substance name used for record lookup.
    pub substance: String,
    /// Thermochemistry library that supplied the record.
    pub library: String,
    /// Exact record key selected in that library.
    pub record_key: String,
    /// Lookup-priority label describing how the record was resolved.
    pub lookup_priority: String,
}

/// Immutable declaration of one resolved physical phase.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct EquilibriumPhaseSpecSnapshot {
    /// Semantic phase name.
    pub phase: String,
    /// Physical state (`Gas`, `Liquid`, ...) selected for the phase.
    pub physical_state: String,
    /// Domain-level phase model label.
    pub model: String,
    /// Component labels belonging to this phase, in declared order.
    pub components: Vec<String>,
}

/// Read-only structural evidence for the repository used by a run.
///
/// This is a catalog-structure fingerprint, not a payload-content hash. A
/// production data release should additionally supply `data_release_label`.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct ThermoCatalogSnapshot {
    /// Stable structural fingerprint over all indexed/payload record identities.
    pub structure_fingerprint: u64,
    /// Number of index entries in the catalog.
    pub indexed_pair_count: usize,
    /// Number of unique index entries.
    pub unique_indexed_pair_count: usize,
    /// Number of payload entries.
    pub payload_pair_count: usize,
    /// Number of duplicated index pairs.
    pub duplicate_index_pair_count: usize,
    /// Number of indexed records without a matching payload.
    pub indexed_without_payload_count: usize,
    /// Number of payload records without a matching index.
    pub payload_without_index_count: usize,
    /// Whether the catalog is structurally consistent.
    pub consistent: bool,
}

/// Candidate-selection evidence when a top-level element query built the run.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct EquilibriumCandidateSelectionSnapshot {
    /// Elements requested by the selection transaction.
    pub requested_elements: Vec<String>,
    /// Element-search mode label.
    pub element_mode: String,
    /// Library preference order used by selection.
    pub library_preference: Vec<String>,
    /// Optional physical-state filter applied during selection.
    pub physical_states: Option<Vec<String>>,
    /// Optional temperature-domain filter in K.
    pub temperature_range_kelvin: Option<(f64, f64)>,
    /// Optional cap on the number of selected candidates.
    pub max_candidates: Option<usize>,
    /// Selected candidate records in deterministic order.
    pub selected_records: Vec<EquilibriumCandidateRecordSnapshot>,
    /// Number of candidates rejected during selection.
    pub rejected_record_count: usize,
    /// Full deterministic rejection evidence retained for replay/audit.
    pub rejected_records: Vec<EquilibriumCandidateRejectionSnapshot>,
    /// Whether the candidate limit excluded otherwise eligible records.
    pub truncated: bool,
}

/// Portable rejection evidence for one candidate-selection record.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct EquilibriumCandidateRejectionSnapshot {
    /// Bare substance name considered by the selector.
    pub substance: String,
    /// Library containing the considered record.
    pub library: String,
    /// Stable selector reason label.
    pub reason: String,
}

/// Candidate record identity retained before phase-plan construction.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct EquilibriumCandidateRecordSnapshot {
    /// Bare substance name.
    pub substance: String,
    /// Thermochemistry library that supplied the record.
    pub library: String,
    /// Exact record key selected in that library.
    pub record_key: String,
    /// Optional physical state selected for this record.
    pub physical_state: Option<String>,
    /// Elements present in this record's composition.
    pub elements: Vec<String>,
    /// Temperature-support interval label reported by the library.
    pub temperature_support: String,
    /// Zero-based priority rank within the library preference order.
    pub library_rank: usize,
}

/// Portable metadata for one accepted fixed-`P,T` equilibrium result.
///
/// The capsule is intentionally descriptive, never executable: it has no
/// write path, no live closures, and no mutable repository. Consumers may
/// serialize it alongside a result, then reconstruct a fresh request through
/// the public phase/candidate APIs under the recorded data release.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct EquilibriumReproducibilityCapsule {
    pub schema_version: u32,
    pub phase_stability_semantics: PhaseStabilitySemantics,
    /// Physical input origin. Optional for compatibility with older schema v2
    /// artifacts written before input provenance was exported.
    #[serde(default)]
    pub input_kind: Option<EquilibriumReproducibilityInputKind>,
    /// Full canonical element order paired by index with `canonical_b`.
    #[serde(default)]
    pub canonical_element_labels: Vec<String>,
    /// Physical conserved elemental totals, before numerical trace seeding.
    #[serde(default)]
    pub canonical_b: Vec<f64>,
    pub temperature_kelvin: f64,
    pub pressure_pa: f64,
    pub reference_pressure_pa: f64,
    pub layout_fingerprint: u64,
    /// Stable FNV-1a identity hash over the exact selected component records.
    pub selected_record_identity_fingerprint: u64,
    /// Effective NIST fallback policy label used by the run.
    pub nist_fallback_policy: String,
    /// Declared physical phases in canonical order.
    pub phases: Vec<EquilibriumPhaseSpecSnapshot>,
    /// Exact selected record identities consumed by the run.
    pub selected_records: Vec<EquilibriumRecordIdentity>,
    /// Effective numerical options snapshot submitted to the pipeline.
    pub solve_options: EquilibriumSolveOptionsSnapshot,
    /// Backend that accepted the final candidate.
    pub accepted_backend: String,
    /// Accepted residual L2 norm.
    pub residual_l2_norm: f64,
    /// Maximum absolute elemental-balance error of the accepted result.
    pub max_abs_element_balance_error: f64,
    /// Optional user/release-supplied payload manifest identifier.
    pub data_release_label: Option<String>,
    /// Repository structural evidence, when a live repository was used.
    pub catalog: Option<ThermoCatalogSnapshot>,
    /// Candidate-selection evidence, when an element query built the run.
    pub candidate_selection: Option<EquilibriumCandidateSelectionSnapshot>,
}

impl EquilibriumReproducibilityCapsule {
    /// Captures one accepted outcome and the exact options clone submitted to
    /// the pipeline. Pipeline requests consume their options by design, so
    /// callers who need a capsule should keep an inexpensive clone.
    pub fn from_outcome(
        outcome: &ResolvedPhaseEquilibriumOutcome,
        solve_options: &EquilibriumSolveOptions,
    ) -> Self {
        let solution = outcome.solution();
        let presentation = crate::Thermodynamics::ChemEquilibrium::equilibrium_presentation::
            EquilibriumPresentationReport::from_solution(solution);
        let selected_records = presentation
            .components
            .into_iter()
            .map(|component| EquilibriumRecordIdentity {
                component: component.component,
                phase: component.phase,
                substance: component.substance,
                library: component.library,
                record_key: component.record_key,
                lookup_priority: component.lookup_priority,
            })
            .collect::<Vec<_>>();
        let phases = outcome
            .resolved()
            .phase_specs()
            .iter()
            .map(|phase| EquilibriumPhaseSpecSnapshot {
                phase: phase
                    .id()
                    .as_option()
                    .clone()
                    .unwrap_or_else(|| "single".to_string()),
                physical_state: format!("{:?}", phase.physical_state()),
                model: format!("{:?}", phase.model()),
                components: phase.components().to_vec(),
            })
            .collect();
        let validation = solution.accepted_solution().validation();
        let conditions = solution.conditions();
        let build_report = solution.build_report();
        let input_kind = Some(match build_report.input_kind() {
            PhaseEquilibriumInputKind::ExplicitComposition => {
                EquilibriumReproducibilityInputKind::ExplicitComposition
            }
            PhaseEquilibriumInputKind::ElementInventory => {
                EquilibriumReproducibilityInputKind::ElementInventory
            }
        });
        Self {
            schema_version: EQUILIBRIUM_REPRODUCIBILITY_SCHEMA_VERSION,
            phase_stability_semantics: PhaseStabilitySemantics::CanonicalTpdV1,
            input_kind,
            canonical_element_labels: build_report.element_labels().to_vec(),
            canonical_b: build_report.element_totals().to_vec(),
            temperature_kelvin: conditions.temperature(),
            pressure_pa: conditions.pressure(),
            reference_pressure_pa: conditions.reference_pressure(),
            layout_fingerprint: solution.metadata().layout_fingerprint(),
            selected_record_identity_fingerprint: selected_record_fingerprint(&selected_records),
            nist_fallback_policy: format!("{:?}", outcome.lookup_report().nist_fallback_policy()),
            phases,
            selected_records,
            solve_options: solve_options.reproducibility_snapshot(),
            accepted_backend: format!("{:?}", solution.solve_report().accepted_backend),
            residual_l2_norm: validation.residual_l2_norm,
            max_abs_element_balance_error: validation.max_abs_element_balance_error,
            data_release_label: None,
            catalog: None,
            candidate_selection: build_report
                .candidate_selection()
                .map(EquilibriumCandidateSelectionSnapshot::from_report),
        }
    }

    /// Attaches a reviewed immutable data-release label chosen by the caller.
    pub fn with_data_release_label(mut self, label: impl Into<String>) -> Self {
        let label = label.into();
        self.data_release_label = (!label.trim().is_empty()).then_some(label);
        self
    }

    /// Attaches repository structural evidence without performing I/O.
    pub fn with_catalog_consistency(mut self, report: &ThermoCatalogConsistencyReport) -> Self {
        self.catalog = Some(ThermoCatalogSnapshot::from_report(report));
        self
    }

    /// Attaches the selection evidence when the input came from an element
    /// query. Explicit `PhaseSpec` callers simply leave this field absent.
    pub fn with_candidate_selection(
        mut self,
        selection: &EquilibriumCandidateSelectionReport,
    ) -> Self {
        self.candidate_selection = Some(EquilibriumCandidateSelectionSnapshot::from_report(
            selection,
        ));
        self
    }

    /// Produces stable, human-readable JSON for a sidecar result artifact.
    pub fn to_pretty_json(&self) -> Result<String, serde_json::Error> {
        serde_json::to_string_pretty(self)
    }

    /// Loads only the current capsule contract.
    ///
    /// Replaying an old physical score under the new TPD name is worse than a
    /// hard error: it would create a plausible but physically ambiguous audit
    /// trail. Migrations must therefore be explicit application-level work.
    pub fn from_json(json: &str) -> Result<Self, ReproducibilityCapsuleError> {
        let value: serde_json::Value =
            serde_json::from_str(json).map_err(ReproducibilityCapsuleError::Json)?;
        let found = value
            .get("schema_version")
            .and_then(serde_json::Value::as_u64);
        if found != Some(u64::from(EQUILIBRIUM_REPRODUCIBILITY_SCHEMA_VERSION)) {
            return Err(ReproducibilityCapsuleError::UnsupportedSchema {
                found,
                expected: EQUILIBRIUM_REPRODUCIBILITY_SCHEMA_VERSION,
            });
        }
        if let Some(field) = find_obsolete_stability_field(&value) {
            return Err(ReproducibilityCapsuleError::ObsoleteStabilityField { field });
        }
        serde_json::from_value(value).map_err(ReproducibilityCapsuleError::Json)
    }
}

fn find_obsolete_stability_field(value: &serde_json::Value) -> Option<String> {
    const OBSOLETE: &[&str] = &[
        "driving_force",
        "phase_stability_model",
        "PhaseStabilityModel",
    ];
    match value {
        serde_json::Value::Object(map) => map.iter().find_map(|(key, value)| {
            OBSOLETE
                .contains(&key.as_str())
                .then(|| key.clone())
                .or_else(|| find_obsolete_stability_field(value))
        }),
        serde_json::Value::Array(values) => values.iter().find_map(find_obsolete_stability_field),
        _ => None,
    }
}

impl ThermoCatalogSnapshot {
    /// Projects a repository consistency report into an immutable snapshot.
    ///
    /// Derives a stable structural fingerprint from the indexed/unique/payload
    /// counts plus every duplicate, missing, and orphan record identity, then
    /// retains the discrete counts and the consistency flag. No repository I/O
    /// is performed; the snapshot only captures already-collected evidence.
    fn from_report(report: &ThermoCatalogConsistencyReport) -> Self {
        let mut identities = Vec::new();
        identities.push(format!("indexed={}", report.indexed_pair_count()));
        identities.push(format!("unique={}", report.unique_indexed_pair_count()));
        identities.push(format!("payload={}", report.payload_pair_count()));
        identities.extend(
            report
                .duplicate_index_pairs()
                .iter()
                .map(|(library, substance)| format!("duplicate:{library}:{substance}")),
        );
        identities.extend(
            report
                .indexed_without_payload()
                .iter()
                .map(|(library, substance)| format!("missing:{library}:{substance}")),
        );
        identities.extend(
            report
                .payload_without_index()
                .iter()
                .map(|(library, substance)| format!("orphan:{library}:{substance}")),
        );
        Self {
            structure_fingerprint: stable_fingerprint(identities.iter().map(String::as_str)),
            indexed_pair_count: report.indexed_pair_count(),
            unique_indexed_pair_count: report.unique_indexed_pair_count(),
            payload_pair_count: report.payload_pair_count(),
            duplicate_index_pair_count: report.duplicate_index_pairs().len(),
            indexed_without_payload_count: report.indexed_without_payload().len(),
            payload_without_index_count: report.payload_without_index().len(),
            consistent: report.is_consistent(),
        }
    }
}

impl EquilibriumCandidateSelectionSnapshot {
    /// Projects an auditable candidate-selection report into an immutable
    /// snapshot for reproducibility.
    ///
    /// Captures the effective selection policy (element mode, library
    /// preference, physical states, temperature range, candidate cap) and the
    /// selected record identities together with the rejected-record count. This
    /// preserves exactly how records were chosen without holding a live
    /// repository handle.
    fn from_report(report: &EquilibriumCandidateSelectionReport) -> Self {
        let policy = report.policy();
        Self {
            requested_elements: report.requested_elements().to_vec(),
            element_mode: format!("{:?}", policy.element_mode()),
            library_preference: policy.library_preference().to_vec(),
            physical_states: policy
                .physical_states()
                .map(|states| states.iter().map(|state| format!("{state:?}")).collect()),
            temperature_range_kelvin: policy
                .temperature_range()
                .map(|range| (range.lower(), range.upper())),
            max_candidates: policy.max_candidates(),
            selected_records: report
                .selected()
                .iter()
                .map(|candidate| EquilibriumCandidateRecordSnapshot {
                    substance: candidate.substance().to_string(),
                    library: candidate.library().to_string(),
                    record_key: candidate.record_key().to_string(),
                    physical_state: candidate.physical_state().map(|state| format!("{state:?}")),
                    elements: candidate.elements().to_vec(),
                    temperature_support: format!("{:?}", candidate.temperature_support()),
                    library_rank: candidate.library_rank(),
                })
                .collect(),
            rejected_record_count: report.rejected().len(),
            rejected_records: report
                .rejected()
                .iter()
                .map(|rejection| EquilibriumCandidateRejectionSnapshot {
                    substance: rejection.substance().to_string(),
                    library: rejection.library().to_string(),
                    reason: format!("{:?}", rejection.reason()),
                })
                .collect(),
            truncated: report.is_truncated(),
        }
    }
}

fn selected_record_fingerprint(records: &[EquilibriumRecordIdentity]) -> u64 {
    stable_fingerprint(records.iter().flat_map(|record| {
        [
            record.component.as_str(),
            record.phase.as_str(),
            record.substance.as_str(),
            record.library.as_str(),
            record.record_key.as_str(),
            record.lookup_priority.as_str(),
        ]
    }))
}

fn stable_fingerprint<'a>(parts: impl IntoIterator<Item = &'a str>) -> u64 {
    parts
        .into_iter()
        .flat_map(|part| part.bytes().chain(std::iter::once(0xff)))
        .fold(0xcbf2_9ce4_8422_2325_u64, |hash, byte| {
            (hash ^ u64::from(byte)).wrapping_mul(0x0000_0100_0000_01b3)
        })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_element_inventory::ElementInventory;
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_log_moles::Solvers;
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_problem::EquilibriumConditions;
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_solver_policy::{
        SolverBackend, SolverPolicy,
    };
    use crate::Thermodynamics::ChemEquilibrium::phase_equilibrium_workflow::PhaseEquilibriumPipelineRequest;
    use crate::Thermodynamics::User_PhaseOrSolution::{
        SubstanceSystemSpecBuilder, SubstancesContainer,
    };
    use crate::Thermodynamics::thermo_lib_api::ThermoData;

    #[test]
    fn local_outcome_exports_stable_policy_provenance_and_catalog_evidence() {
        let repository = ThermoData::try_default_repository().unwrap();
        let spec = SubstanceSystemSpecBuilder::new(SubstancesContainer::SinglePhase(vec![
            "N2".to_string(),
            "O2".to_string(),
        ]))
        .with_library_priorities(vec!["NASA_gas".to_string()])
        .with_search_in_nist(false)
        .build()
        .unwrap();
        let options = EquilibriumSolveOptions::new()
            .with_solver_policy(SolverPolicy::Single(SolverBackend::Legacy(Solvers::NR)))
            .unwrap();
        let outcome = PhaseEquilibriumPipelineRequest::new(
            spec,
            vec![0.79, 0.21],
            EquilibriumConditions::new(500.0, 101_325.0, 101_325.0).unwrap(),
        )
        .with_repository(repository.clone())
        .with_solve_options(options.clone())
        .solve()
        .unwrap();

        let capsule = EquilibriumReproducibilityCapsule::from_outcome(&outcome, &options)
            .with_data_release_label("bundled-local-test-data")
            .with_catalog_consistency(&repository.consistency_report());
        let repeated = EquilibriumReproducibilityCapsule::from_outcome(&outcome, &options)
            .with_data_release_label("bundled-local-test-data")
            .with_catalog_consistency(&repository.consistency_report());
        assert_eq!(capsule, repeated);
        assert_eq!(
            capsule.schema_version,
            EQUILIBRIUM_REPRODUCIBILITY_SCHEMA_VERSION
        );
        assert_eq!(capsule.selected_records.len(), 2);
        assert_eq!(
            capsule.input_kind,
            Some(EquilibriumReproducibilityInputKind::ExplicitComposition)
        );
        assert_eq!(capsule.canonical_element_labels, ["N", "O"]);
        assert_eq!(capsule.canonical_b, [1.58, 0.42]);
        assert_eq!(
            capsule.solve_options.effective_backend_order,
            vec!["Legacy(NR)"]
        );
        assert!(capsule.selected_record_identity_fingerprint != 0);
        assert!(capsule.catalog.is_some());
        let json = capsule.to_pretty_json().unwrap();
        assert!(json.contains("NASA_gas"));
        assert!(json.contains("bundled-local-test-data"));
        assert_eq!(
            EquilibriumReproducibilityCapsule::from_json(&json).unwrap(),
            capsule
        );

        // Phase models are stored by their explicit symbolic names, not by a
        // fragile enum ordinal. Adding `IdealSolution` therefore does not
        // reinterpret existing `IdealGas`/`PureCondensed` artifacts or force
        // a schema bump by itself.
        let mut ideal_solution_named = capsule;
        ideal_solution_named.phases[0].model = "IdealSolution".to_string();
        let round_trip = EquilibriumReproducibilityCapsule::from_json(
            &ideal_solution_named.to_pretty_json().unwrap(),
        )
        .unwrap();
        assert_eq!(round_trip.phases[0].model, "IdealSolution");
        assert_eq!(
            round_trip.phase_stability_semantics,
            PhaseStabilitySemantics::CanonicalTpdV1
        );
    }

    #[test]
    fn elemental_outcome_exports_input_kind_and_canonical_b() {
        let repository = ThermoData::try_default_repository().unwrap();
        let spec = SubstanceSystemSpecBuilder::new(SubstancesContainer::SinglePhase(vec![
            "H2".to_string(),
            "O2".to_string(),
            "H2O".to_string(),
        ]))
        .with_library_priorities(vec!["NASA_gas".to_string()])
        .with_search_in_nist(false)
        .build()
        .unwrap();
        let options = EquilibriumSolveOptions::new()
            .with_solver_policy(SolverPolicy::Single(SolverBackend::Legacy(Solvers::NR)))
            .unwrap();
        let outcome = PhaseEquilibriumPipelineRequest::from_element_inventory(
            spec,
            ElementInventory::from_amounts([("H", 4.0), ("O", 2.0)]).unwrap(),
            EquilibriumConditions::new(1_200.0, 101_325.0, 101_325.0).unwrap(),
        )
        .with_repository(repository)
        .with_solve_options(options.clone())
        .solve()
        .unwrap();

        let capsule = EquilibriumReproducibilityCapsule::from_outcome(&outcome, &options);
        assert_eq!(
            capsule.input_kind,
            Some(EquilibriumReproducibilityInputKind::ElementInventory)
        );
        assert_eq!(capsule.canonical_element_labels, ["H", "O"]);
        assert_eq!(capsule.canonical_b, [4.0, 2.0]);
        assert!(capsule.candidate_selection.is_none());

        let json = capsule.to_pretty_json().unwrap();
        assert!(json.contains("ElementInventory"));
        assert!(json.contains("canonical_b"));
        assert_eq!(
            EquilibriumReproducibilityCapsule::from_json(&json).unwrap(),
            capsule
        );

        let mut legacy_value = serde_json::to_value(&capsule).unwrap();
        let legacy_object = legacy_value.as_object_mut().unwrap();
        legacy_object.remove("input_kind");
        legacy_object.remove("canonical_element_labels");
        legacy_object.remove("canonical_b");
        let legacy = EquilibriumReproducibilityCapsule::from_json(
            &serde_json::to_string(&legacy_value).unwrap(),
        )
        .unwrap();
        assert_eq!(legacy.input_kind, None);
        assert!(legacy.canonical_element_labels.is_empty());
        assert!(legacy.canonical_b.is_empty());
    }

    #[test]
    fn reproducibility_loader_rejects_old_schema_and_obsolete_stability_fields() {
        let old_schema = serde_json::json!({ "schema_version": 1 });
        assert!(matches!(
            EquilibriumReproducibilityCapsule::from_json(&old_schema.to_string()),
            Err(ReproducibilityCapsuleError::UnsupportedSchema { found: Some(1), .. })
        ));

        let obsolete_field = serde_json::json!({
            "schema_version": EQUILIBRIUM_REPRODUCIBILITY_SCHEMA_VERSION,
            "phase_stability": { "driving_force": -1.0 }
        });
        assert!(matches!(
            EquilibriumReproducibilityCapsule::from_json(&obsolete_field.to_string()),
            Err(ReproducibilityCapsuleError::ObsoleteStabilityField { field }) if field == "driving_force"
        ));
    }

    #[test]
    fn options_snapshot_expands_the_implicit_production_cascade() {
        let snapshot = EquilibriumSolveOptions::new().reproducibility_snapshot();
        assert!(
            snapshot
                .effective_backend_order
                .iter()
                .any(|backend| backend.contains("RustedSciThe"))
        );
        assert!(
            snapshot
                .effective_backend_order
                .iter()
                .any(|backend| backend.contains("Legacy"))
        );
        assert_eq!(snapshot.timing_mode, "Disabled");
        assert!(!snapshot.execution_control_attached);
    }
}
