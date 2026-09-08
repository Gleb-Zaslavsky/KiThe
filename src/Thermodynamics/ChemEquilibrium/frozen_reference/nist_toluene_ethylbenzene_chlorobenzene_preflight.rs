//! Capability and source audit for the first ternary VLE candidate.
//!
//! This is deliberately not a TPD benchmark.  It records whether KiThe can
//! represent the six exact gas/liquid standard states required by the
//! toluene/ethylbenzene/chlorobenzene ThermoML experiment.  A missing state is
//! a typed `ValidationNotApplicable` result, never permission to substitute a
//! chemically related record or to query NIST at test time.

use crate::Thermodynamics::ChemEquilibrium::equilibrium_nonlinear::ReactionExtentError;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_ph_thermochemistry::ResolvedThermochemistry;
use crate::Thermodynamics::User_PhaseOrSolution::{
    PhaseModel, PhaseSpec, SubstanceSystemFactory, SubstanceSystemSpec,
};
use crate::Thermodynamics::phase_layout::PhaseId;
use crate::Thermodynamics::physical_state::PhysicalState;

const NIST_THERMOML_DOI: &str = "10.1021/je020186c";
const NIST_THERMOML_URL: &str = "https://trc.nist.gov/ThermoML/10.1021/je020186c.html";

/// Exact local identity required for one state-qualified ternary component.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct TernaryVleRecordRequest {
    /// Compound name.
    pub(crate) compound: &'static str,
    /// Molecular formula.
    pub(crate) formula: &'static str,
    /// Library expected to hold the exact record.
    pub(crate) library: &'static str,
    /// Exact repository record key.
    pub(crate) record_key: &'static str,
    /// Physical state of the standard state.
    pub(crate) state: PhysicalState,
    /// Activity model required for the phase.
    pub(crate) model: PhaseModel,
}

/// Read-only result of resolving one exact local state.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct TernaryVleLocalRecordAudit {
    /// The exact record requested by the preflight.
    pub(crate) request: TernaryVleRecordRequest,
    /// Library actually selected after resolution, when available.
    pub(crate) selected_library: Option<String>,
    /// Record key actually selected after resolution, when available.
    pub(crate) selected_record_key: Option<String>,
    /// Resolved physical-state label, when available.
    pub(crate) selected_state: Option<String>,
    /// Whether the resolved provenance matches the requested library/record.
    pub(crate) phase_model_compatible: bool,
    /// Native temperature interval, K, when the record resolved.
    pub(crate) temperature_interval_k: Option<(f64, f64)>,
    /// Whether the record provides standard Gibbs capability.
    pub(crate) supports_standard_gibbs: bool,
    /// Whether the record provides enthalpy capability.
    pub(crate) supports_enthalpy: bool,
    /// Standard-state pressure provenance, when available.
    pub(crate) standard_state_pressure_provenance: Option<String>,
    /// Reason the record is unavailable, or `None` when it resolved.
    pub(crate) unavailable_reason: Option<String>,
}

impl TernaryVleLocalRecordAudit {
    pub(crate) fn is_available(&self) -> bool {
        self.unavailable_reason.is_none()
    }
}

/// Machine-audited facts from the official NIST ThermoML payload.
///
/// Dataset 8 of the payload contains 48 ternary rows.  The property is a
/// boiling temperature; the variables are liquid mole fractions of toluene
/// and ethylbenzene plus pressure.  Chlorobenzene is the implied third liquid
/// simplex coordinate.  No experimental vapour composition is published.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct NistTernaryThermoMlAudit {
    /// DOI of the source publication.
    pub(crate) doi: &'static str,
    /// URL of the official ThermoML record.
    pub(crate) source_url: &'static str,
    /// Component order in the ThermoML dataset.
    pub(crate) component_order: [&'static str; 3],
    /// Measured property (boiling temperature at pressure).
    pub(crate) property: &'static str,
    /// Published pressure levels, kPa.
    pub(crate) pressure_levels_kpa: [f64; 4],
    /// Number of ternary rows in the payload.
    pub(crate) ternary_row_count: usize,
    /// Whether liquid composition was published.
    pub(crate) liquid_composition_published: bool,
    /// Whether experimental vapor composition was published.
    pub(crate) vapor_composition_published: bool,
    /// Whether the property uncertainty was published.
    pub(crate) property_uncertainty_published: bool,
    /// Whether the same payload also contains pure saturation support.
    pub(crate) same_payload_has_pure_saturation_support: bool,
}

/// Final preflight decision.  `Ready` is intentionally impossible until all
/// six local records have both numeric thermochemistry capabilities and one
/// common native temperature interval.
#[derive(Debug, Clone, PartialEq)]
pub(crate) enum TernaryVlePreflightDecision {
    /// All six exact records resolved with one common native temperature interval.
    Ready {
        /// Common (lower, upper) temperature interval over all records, K.
        common_temperature_interval_k: (f64, f64),
        /// Number of liquid components in the candidate.
        liquid_component_count: usize,
        /// Dimension of the liquid composition simplex.
        liquid_simplex_dimension: usize,
    },
    /// One or more exact records are missing or incompatible.
    ValidationNotApplicable {
        /// Description of the missing/incompatible records.
        missing_or_incompatible_records: Vec<String>,
        /// Human-readable reason for the negative result.
        reason: String,
    },
}

/// Full local-plus-external capability audit for the proposed ternary story.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct TernaryVlePreflightReport {
    /// Per-record local capability audits.
    pub(crate) local_records: Vec<TernaryVleLocalRecordAudit>,
    /// Machine-audited facts from the official NIST ThermoML payload.
    pub(crate) source: NistTernaryThermoMlAudit,
    /// Final preflight decision.
    pub(crate) decision: TernaryVlePreflightDecision,
}

impl TernaryVlePreflightReport {
    /// Human-readable audit summary for an explicitly requested diagnostic
    /// run. Normal tests keep it captured; a maintainer can use `--nocapture`
    /// while selecting the next local ternary candidate.
    pub(crate) fn summary(&self) -> String {
        let mut lines = vec![format!(
            "NIST ternary VLE preflight: DOI={} rows={} property={} vapor_x={} pure_psat={}",
            self.source.doi,
            self.source.ternary_row_count,
            self.source.property,
            self.source.vapor_composition_published,
            self.source.same_payload_has_pure_saturation_support,
        )];
        for record in &self.local_records {
            let availability = if let Some(reason) = &record.unavailable_reason {
                format!("unavailable ({reason})")
            } else {
                format!(
                    "{}:{} state={} T={:?} G={} H={} p0={}",
                    record.selected_library.as_deref().unwrap_or("?"),
                    record.selected_record_key.as_deref().unwrap_or("?"),
                    record.selected_state.as_deref().unwrap_or("?"),
                    record.temperature_interval_k,
                    record.supports_standard_gibbs,
                    record.supports_enthalpy,
                    record
                        .standard_state_pressure_provenance
                        .as_deref()
                        .unwrap_or("?"),
                )
            };
            lines.push(format!(
                "  {} {:?}: {} [{}]",
                record.request.compound,
                record.request.state,
                availability,
                record.request.record_key
            ));
        }
        lines.push(format!("decision={:?}", self.decision));
        lines.join("\n")
    }

    /// Converts a negative preflight result into the normal typed validation
    /// boundary used by the equilibrium suite.
    pub(crate) fn require_executable(&self) -> Result<(), ReactionExtentError> {
        match &self.decision {
            TernaryVlePreflightDecision::Ready { .. } => Ok(()),
            TernaryVlePreflightDecision::ValidationNotApplicable { reason, .. } => {
                Err(ReactionExtentError::ValidationNotApplicable {
                    path: "nist_toluene_ethylbenzene_chlorobenzene_ternary_vle",
                    message: reason.clone(),
                })
            }
        }
    }
}

const PRIMARY_RECORDS: [TernaryVleRecordRequest; 6] = [
    TernaryVleRecordRequest {
        compound: "toluene",
        formula: "C7H8",
        library: "NASA_gas",
        record_key: "C7H8",
        state: PhysicalState::Gas,
        model: PhaseModel::IdealGas,
    },
    TernaryVleRecordRequest {
        compound: "ethylbenzene",
        formula: "C8H10",
        library: "NASA_gas",
        record_key: "C8H10,ethylbenz",
        state: PhysicalState::Gas,
        model: PhaseModel::IdealGas,
    },
    TernaryVleRecordRequest {
        compound: "chlorobenzene",
        formula: "C6H5Cl",
        library: "NASA_gas",
        record_key: "C6H5Cl",
        state: PhysicalState::Gas,
        model: PhaseModel::IdealGas,
    },
    TernaryVleRecordRequest {
        compound: "toluene",
        formula: "C7H8",
        library: "NASA_cond",
        record_key: "C7H8(L)",
        state: PhysicalState::Liquid,
        model: PhaseModel::IdealSolution,
    },
    TernaryVleRecordRequest {
        compound: "ethylbenzene",
        formula: "C8H10",
        library: "NASA_cond",
        record_key: "C8H10(L),ethylbenz",
        state: PhysicalState::Liquid,
        model: PhaseModel::IdealSolution,
    },
    TernaryVleRecordRequest {
        compound: "chlorobenzene",
        formula: "C6H5Cl",
        library: "NASA_cond",
        record_key: "C6H5Cl(L)",
        state: PhysicalState::Liquid,
        model: PhaseModel::IdealSolution,
    },
];

/// Returns the exact, manually reviewed ThermoML source shape.
pub(crate) fn nist_ternary_thermoml_audit() -> NistTernaryThermoMlAudit {
    NistTernaryThermoMlAudit {
        doi: NIST_THERMOML_DOI,
        source_url: NIST_THERMOML_URL,
        component_order: ["toluene", "ethylbenzene", "chlorobenzene"],
        property: "boiling temperature at pressure P, K",
        pressure_levels_kpa: [26.66, 53.33, 79.99, 101.32],
        ternary_row_count: 48,
        liquid_composition_published: true,
        vapor_composition_published: false,
        property_uncertainty_published: true,
        same_payload_has_pure_saturation_support: false,
    }
}

/// Resolves only one state-qualified record with the requested local library.
///
/// A one-component `IdealSolution` here is a capability probe, not a claim
/// that a future ternary liquid may be reduced to a pure-condensed phase.
/// Resolves one state-qualified record with the requested library and probes its
/// G(T)/H(T) capabilities, returning an audit even on failure.
fn audit_local_record(request: TernaryVleRecordRequest) -> TernaryVleLocalRecordAudit {
    let phase = match request.model {
        PhaseModel::IdealGas => PhaseSpec::ideal_gas(
            PhaseId::new(Some(format!("preflight:{}:gas", request.compound))),
            vec![request.record_key.to_owned()],
        ),
        PhaseModel::IdealSolution => PhaseSpec::ideal_solution(
            PhaseId::new(Some(format!("preflight:{}:liquid", request.compound))),
            vec![request.record_key.to_owned()],
            request.state,
        ),
        _ => unreachable!("the ternary VLE preflight only declares ideal phases"),
    };
    let phase = match phase {
        Ok(phase) => phase,
        Err(error) => {
            return unavailable(
                request,
                format!("invalid static phase declaration: {error}"),
            );
        }
    };
    let spec = match SubstanceSystemSpec::from_phases(vec![phase]) {
        Ok(spec) => spec.with_lookup_policy(
            vec![request.library.to_owned()],
            vec![request.library.to_owned()],
            None,
            false,
        ),
        Err(error) => return unavailable(request, format!("invalid local lookup spec: {error}")),
    };
    let resolved = match SubstanceSystemFactory::resolve_spec(spec) {
        Ok(resolved) => resolved,
        Err(error) => {
            return unavailable(
                request,
                format!(
                    "exact local {:?} record '{}' in {} did not resolve with NIST fallback disabled: {error}",
                    request.state, request.record_key, request.library
                ),
            );
        }
    };
    let thermochemistry = match ResolvedThermochemistry::from_resolved_system(&resolved) {
        Ok(thermochemistry) => thermochemistry,
        Err(error) => {
            return unavailable(
                request,
                format!(
                    "exact local {:?} record '{}' resolved but lacks usable G(T)/H(T) capabilities: {error}",
                    request.state, request.record_key
                ),
            );
        }
    };
    let Some(provenance) = thermochemistry.provenance().first() else {
        return unavailable(request, "resolved record produced no provenance".to_owned());
    };
    let bounds = thermochemistry.temperature_bounds();
    TernaryVleLocalRecordAudit {
        request,
        selected_library: Some(provenance.library().to_owned()),
        selected_record_key: Some(provenance.record_key().to_owned()),
        selected_state: Some(provenance.state().to_owned()),
        phase_model_compatible: provenance.library() == request.library
            && provenance.record_key() == request.record_key,
        temperature_interval_k: Some((bounds.lower(), bounds.upper())),
        supports_standard_gibbs: true,
        supports_enthalpy: true,
        standard_state_pressure_provenance: Some(provenance.standard_state_pressure().to_string()),
        unavailable_reason: None,
    }
}

/// Builds an unavailable audit row for a failed record request.
fn unavailable(
    request: TernaryVleRecordRequest,
    unavailable_reason: String,
) -> TernaryVleLocalRecordAudit {
    TernaryVleLocalRecordAudit {
        request,
        selected_library: None,
        selected_record_key: None,
        selected_state: None,
        phase_model_compatible: false,
        temperature_interval_k: None,
        supports_standard_gibbs: false,
        supports_enthalpy: false,
        standard_state_pressure_provenance: None,
        unavailable_reason: Some(unavailable_reason),
    }
}

/// Performs the complete preflight without contacting NIST or mutating a
/// library.  This must remain a cheap capability gate before any future 2D
/// TPD/minimization fixture is constructed.
pub(crate) fn preflight_primary_ternary_vle() -> TernaryVlePreflightReport {
    let local_records = PRIMARY_RECORDS
        .into_iter()
        .map(audit_local_record)
        .collect::<Vec<_>>();
    let missing_or_incompatible_records = local_records
        .iter()
        .filter(|record| !record.is_available() || !record.phase_model_compatible)
        .map(|record| {
            format!(
                "{} {:?} '{}'/{}",
                record.request.compound,
                record.request.state,
                record.request.library,
                record.request.record_key
            )
        })
        .collect::<Vec<_>>();
    let decision = if missing_or_incompatible_records.is_empty() {
        let lower = local_records
            .iter()
            .filter_map(|record| record.temperature_interval_k.map(|bounds| bounds.0))
            .fold(f64::NEG_INFINITY, f64::max);
        let upper = local_records
            .iter()
            .filter_map(|record| record.temperature_interval_k.map(|bounds| bounds.1))
            .fold(f64::INFINITY, f64::min);
        if lower.is_finite() && upper.is_finite() && lower < upper {
            TernaryVlePreflightDecision::Ready {
                common_temperature_interval_k: (lower, upper),
                liquid_component_count: 3,
                liquid_simplex_dimension: 2,
            }
        } else {
            TernaryVlePreflightDecision::ValidationNotApplicable {
                missing_or_incompatible_records: Vec::new(),
                reason: format!(
                    "six exact records have no common native temperature interval: [{lower}, {upper}] K"
                ),
            }
        }
    } else {
        TernaryVlePreflightDecision::ValidationNotApplicable {
            missing_or_incompatible_records,
            reason: "the exact local six-state inventory is incomplete; do not build a ternary TPD fixture or substitute a related compound".to_owned(),
        }
    };
    TernaryVlePreflightReport {
        local_records,
        source: nist_ternary_thermoml_audit(),
        decision,
    }
}
