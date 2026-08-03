//! Read-only GUI projection of the engine's element-candidate report.
//!
//! Candidate discovery is a repository transaction, not editable document
//! state. This projection keeps the audit columns needed by the editor while
//! deliberately preserving the engine's selected/rejected decision.

use crate::Thermodynamics::ChemEquilibrium::prelude::{
    CandidateRejectionReason, CandidateTemperatureSupport, EquilibriumCandidateSelectionReport,
};
use crate::Thermodynamics::physical_state::PhysicalState;

/// One auditable candidate row shown by the equilibrium editor.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct EquilibriumGuiCandidateRow {
    substance: String,
    record_key: Option<String>,
    library: String,
    physical_state: Option<PhysicalState>,
    temperature_support: CandidateTemperatureSupport,
    included: bool,
    rejection: Option<CandidateRejectionReason>,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Thermodynamics::ChemEquilibrium::prelude::{
        ElementSearchMode, EquilibriumCandidatePolicy, EquilibriumCandidateSelector,
        ThermoRepository,
    };
    use serde_json::json;
    use std::collections::HashMap;
    use std::sync::Arc;

    #[test]
    fn projection_retains_selected_and_rejected_audit_rows() {
        let repository = Arc::new(ThermoRepository::from_parts(
            vec![
                ("NASA_gas".into(), "CO".into()),
                ("NASA_gas".into(), "CO2".into()),
            ],
            HashMap::from([(
                "NASA_gas".into(),
                HashMap::from([
                    ("CO".into(), json!({"T": [[200.0, 6000.0]]})),
                    ("CO2".into(), json!({"T": [[200.0, 6000.0]]})),
                ]),
            )]),
            HashMap::from([
                (
                    "C".into(),
                    vec![
                        vec!["CO".into(), "NASA_gas".into()],
                        vec!["CO2".into(), "NASA_gas".into()],
                    ],
                ),
                (
                    "O".into(),
                    vec![
                        vec!["CO".into(), "NASA_gas".into()],
                        vec!["CO2".into(), "NASA_gas".into()],
                    ],
                ),
            ]),
            vec!["NASA_gas".into()],
            HashMap::new(),
            HashMap::new(),
            vec!["NASA_gas".into()],
            Vec::new(),
        ));
        let report = EquilibriumCandidateSelector::new(repository)
            .select(
                &["C".into(), "O".into()],
                EquilibriumCandidatePolicy::new(ElementSearchMode::ExactSet)
                    .with_max_candidates(1)
                    .expect("positive candidate limit"),
            )
            .expect("synthetic repository selection succeeds");
        let preview = EquilibriumGuiCandidatePreview::from_report(&report);

        assert_eq!(
            preview.requested_elements(),
            &["C".to_string(), "O".to_string()]
        );
        assert_eq!(preview.selected_count(), 1);
        assert_eq!(preview.rejected_count(), 1);
        assert!(preview.rows()[0].included());
        assert!(preview.rows()[0].record_key().is_some());
        assert_eq!(preview.rows()[1].record_key(), None);
        assert!(preview.rows()[1].rejection().is_some());
    }
}

impl EquilibriumGuiCandidateRow {
    pub fn substance(&self) -> &str {
        &self.substance
    }

    pub fn record_key(&self) -> Option<&str> {
        self.record_key.as_deref()
    }

    pub fn library(&self) -> &str {
        &self.library
    }

    pub fn physical_state(&self) -> Option<PhysicalState> {
        self.physical_state
    }

    pub fn temperature_support(&self) -> CandidateTemperatureSupport {
        self.temperature_support
    }

    pub fn included(&self) -> bool {
        self.included
    }

    pub fn rejection(&self) -> Option<&CandidateRejectionReason> {
        self.rejection.as_ref()
    }
}

/// Immutable candidate preview tied to one validated element query.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct EquilibriumGuiCandidatePreview {
    requested_elements: Vec<String>,
    rows: Vec<EquilibriumGuiCandidateRow>,
}

impl EquilibriumGuiCandidatePreview {
    /// Converts the complete engine report without dropping rejected evidence.
    pub fn from_report(report: &EquilibriumCandidateSelectionReport) -> Self {
        let mut rows = Vec::with_capacity(report.selected().len() + report.rejected().len());
        rows.extend(
            report
                .selected()
                .iter()
                .map(|candidate| EquilibriumGuiCandidateRow {
                    substance: candidate.substance().to_string(),
                    record_key: Some(candidate.record_key().to_string()),
                    library: candidate.library().to_string(),
                    physical_state: candidate.physical_state(),
                    temperature_support: candidate.temperature_support(),
                    included: true,
                    rejection: None,
                }),
        );
        rows.extend(report.rejected().iter().map(|rejection| {
            EquilibriumGuiCandidateRow {
                substance: rejection.substance().to_string(),
                // The engine currently has no record key for rejected rows;
                // the GUI must not invent a misleading identity.
                record_key: None,
                library: rejection.library().to_string(),
                physical_state: None,
                temperature_support: CandidateTemperatureSupport::Unknown,
                included: false,
                rejection: Some(rejection.reason().clone()),
            }
        }));
        Self {
            requested_elements: report.requested_elements().to_vec(),
            rows,
        }
    }

    pub fn requested_elements(&self) -> &[String] {
        &self.requested_elements
    }

    pub fn rows(&self) -> &[EquilibriumGuiCandidateRow] {
        &self.rows
    }

    pub fn selected_count(&self) -> usize {
        self.rows.iter().filter(|row| row.included).count()
    }

    pub fn rejected_count(&self) -> usize {
        self.rows
            .iter()
            .filter(|row| row.rejection.is_some())
            .count()
    }
}
