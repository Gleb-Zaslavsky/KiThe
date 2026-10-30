//! Non-transactional diagnostics for reaction-library maintenance.
//!
//! Production APIs remain fail-fast and transactional. These audit functions
//! instead retain one outcome per requested reaction so a damaged catalog can
//! be inspected without losing later valid records.
//!
//! A clean report confirms that selected records exist, satisfy the current
//! payload schema, and have parseable equations. It does not replace runtime
//! model checks performed by the production calculation workflow.
//!
//! A typical caller audits first, decides whether the requested set is clean,
//! and then independently uses the ordinary production facade:
//!
//! ```no_run
//! use KiThe::Kinetics::library_audit::{
//!     audit_library_payload_range, audit_parsed_equations,
//! };
//! use KiThe::Kinetics::User_reactions::KinData;
//! use KiThe::Kinetics::error::{KineticsError, KineticsResult};
//!
//! fn prepare_mechanism() -> KineticsResult<KinData> {
//!     let payloads = audit_library_payload_range("Cantera", 1..=3)?;
//!     let equations = audit_parsed_equations(&payloads);
//!     if !payloads.is_clean() || !equations.is_clean() {
//!         return Err(KineticsError::InvalidReactionData(
//!             "selected Cantera reactions failed library audit".into(),
//!         ));
//!     }
//!
//!     let mut kinetics = KinData::new();
//!     kinetics.set_reactions_from_shortcut_range("C1..3".into())?;
//!     kinetics.get_reactions_from_shortcuts()?;
//!     kinetics.reactdata_parsing()?;
//!     kinetics.analyze_reactions()?;
//!     Ok(kinetics)
//! }
//! ```

use std::ops::RangeInclusive;
use std::path::Path;

use serde::{Deserialize, Serialize};
use serde_json::Value;

use crate::Kinetics::error::{KineticsError, KineticsResult};
use crate::Kinetics::kinetics_lib_api::KineticData;
use crate::Kinetics::mechfinder_api::ReactionData;
use crate::Kinetics::stoichiometry_analyzer::StoichAnalyzer;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum AuditStage {
    CatalogLookup,
    PayloadDeserialize,
    PayloadValidation,
    EquationAnalysis,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum AuditErrorKind {
    MissingReaction,
    Json,
    InvalidPayload,
    ReactionSyntax,
    UnknownSpecies,
    Other,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct AuditFailure {
    pub stage: AuditStage,
    pub kind: AuditErrorKind,
    pub message: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PayloadAuditRecord {
    pub reaction_id: String,
    pub range_offset: u64,
    pub raw_payload: Option<Value>,
    pub result: Result<ReactionData, AuditFailure>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PayloadAuditReport {
    pub library: String,
    pub requested_reaction_ids: Vec<String>,
    pub requested_range: Option<RangeInclusive<u64>>,
    pub records: Vec<PayloadAuditRecord>,
}

impl PayloadAuditReport {
    pub fn successful_count(&self) -> usize {
        self.records
            .iter()
            .filter(|record| record.result.is_ok())
            .count()
    }

    pub fn failed_count(&self) -> usize {
        self.records.len() - self.successful_count()
    }

    pub fn is_clean(&self) -> bool {
        self.failed_count() == 0
    }

    pub fn write_json(&self, path: impl AsRef<Path>) -> KineticsResult<()> {
        write_json(path, self)
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EquationParseSummary {
    pub equation: String,
    pub species: Vec<String>,
    pub stoichiometric_coefficients: Vec<f64>,
    pub reactant_coefficients: Vec<f64>,
    pub product_coefficients: Vec<f64>,
    pub reactant_orders: Vec<f64>,
    pub product_orders: Vec<f64>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EquationAuditRecord {
    pub reaction_id: String,
    pub equation: String,
    pub result: Result<EquationParseSummary, AuditFailure>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EquationAuditReport {
    pub library: String,
    pub requested_reaction_ids: Vec<String>,
    pub requested_range: Option<RangeInclusive<u64>>,
    pub records: Vec<EquationAuditRecord>,
}

impl EquationAuditReport {
    pub fn successful_count(&self) -> usize {
        self.records
            .iter()
            .filter(|record| record.result.is_ok())
            .count()
    }

    pub fn failed_count(&self) -> usize {
        self.records.len() - self.successful_count()
    }

    pub fn is_clean(&self) -> bool {
        self.failed_count() == 0
    }

    pub fn write_json(&self, path: impl AsRef<Path>) -> KineticsResult<()> {
        write_json(path, self)
    }
}

/// Audit every numeric reaction ID in an inclusive range, continuing after
/// missing IDs, JSON-shape errors, and domain validation failures.
pub fn audit_library_payload_range(
    library: &str,
    range: RangeInclusive<u64>,
) -> KineticsResult<PayloadAuditReport> {
    if range.is_empty() {
        return Err(KineticsError::InvalidShortcutRange(format!(
            "{}..={}",
            range.start(),
            range.end()
        )));
    }

    let reaction_ids: Vec<_> = range.clone().map(|id| id.to_string()).collect();
    audit_library_reactions_inner(library, &reaction_ids, Some(range))
}

/// Audit an explicit list of canonical reaction IDs from one library.
///
/// Each requested ID receives an outcome, including IDs absent from the
/// library. This returns diagnostics only and never mutates production state.
pub fn audit_library_reactions(
    library: &str,
    reaction_ids: &[String],
) -> KineticsResult<PayloadAuditReport> {
    if reaction_ids.is_empty() {
        return Err(KineticsError::InvalidReactionData(
            "library audit requires at least one reaction ID".to_string(),
        ));
    }
    audit_library_reactions_inner(library, reaction_ids, None)
}

fn audit_library_reactions_inner(
    library: &str,
    reaction_ids: &[String],
    requested_range: Option<RangeInclusive<u64>>,
) -> KineticsResult<PayloadAuditReport> {
    let mut catalog = KineticData::new();
    catalog.open_json_files(library)?;
    let records = reaction_ids
        .iter()
        .enumerate()
        .map(|(offset, reaction_id)| {
            audit_catalog_entry(
                library,
                reaction_id.clone(),
                offset as u64,
                catalog.LibKineticData.get(reaction_id).cloned(),
            )
        })
        .collect();

    Ok(PayloadAuditReport {
        library: library.to_string(),
        requested_reaction_ids: reaction_ids.to_vec(),
        requested_range,
        records,
    })
}

fn audit_catalog_entry(
    library: &str,
    reaction_id: String,
    range_offset: u64,
    raw_payload: Option<Value>,
) -> PayloadAuditRecord {
    let result = match raw_payload.as_ref() {
        None => Err(AuditFailure {
            stage: AuditStage::CatalogLookup,
            kind: AuditErrorKind::MissingReaction,
            message: format!("reaction `{reaction_id}` is absent from library `{library}`"),
        }),
        Some(payload) => audit_payload(payload.clone(), &reaction_id),
    };
    PayloadAuditRecord {
        reaction_id,
        range_offset,
        raw_payload,
        result,
    }
}

/// Parse equations from successful payload records independently. Failed
/// payload records are intentionally excluded; each equation failure is kept
/// beside its source reaction ID while later equations continue to be checked.
pub fn audit_parsed_equations(payloads: &PayloadAuditReport) -> EquationAuditReport {
    let records = payloads
        .records
        .iter()
        .filter_map(|record| {
            let reaction = record.result.as_ref().ok()?;
            Some(EquationAuditRecord {
                reaction_id: record.reaction_id.clone(),
                equation: reaction.eq.clone(),
                result: audit_equation(reaction.eq.clone()),
            })
        })
        .collect();

    EquationAuditReport {
        library: payloads.library.clone(),
        requested_reaction_ids: payloads.requested_reaction_ids.clone(),
        requested_range: payloads.requested_range.clone(),
        records,
    }
}

fn audit_payload(payload: Value, reaction_id: &str) -> Result<ReactionData, AuditFailure> {
    let reaction =
        serde_json::from_value::<ReactionData>(payload).map_err(|error| AuditFailure {
            stage: AuditStage::PayloadDeserialize,
            kind: AuditErrorKind::Json,
            message: error.to_string(),
        })?;
    reaction.validate().map_err(|error| AuditFailure {
        stage: AuditStage::PayloadValidation,
        kind: classify_error(&error),
        message: format!("reaction `{reaction_id}`: {error}"),
    })?;
    Ok(reaction)
}

fn audit_equation(equation: String) -> Result<EquationParseSummary, AuditFailure> {
    let mut analyzer = StoichAnalyzer::new();
    analyzer.reactions = vec![equation.clone()];
    analyzer
        .search_substances()
        .and_then(|()| analyzer.analyse_reactions())
        .map_err(|error| AuditFailure {
            stage: AuditStage::EquationAnalysis,
            kind: classify_error(&error),
            message: error.to_string(),
        })?;

    Ok(EquationParseSummary {
        equation,
        species: analyzer.substances,
        stoichiometric_coefficients: analyzer.stecheo_matrx.remove(0),
        reactant_coefficients: analyzer.stecheo_reags.remove(0),
        product_coefficients: analyzer.stecheo_prods.remove(0),
        reactant_orders: analyzer.G_matrix_reag.remove(0),
        product_orders: analyzer.G_matrix_prod.remove(0),
    })
}

fn classify_error(error: &KineticsError) -> AuditErrorKind {
    match error {
        KineticsError::Context { source, .. } => classify_error(source),
        KineticsError::Json(_) => AuditErrorKind::Json,
        KineticsError::InvalidReactionPayload { .. } | KineticsError::InvalidReactionData(_) => {
            AuditErrorKind::InvalidPayload
        }
        KineticsError::ReactionSyntax { .. } | KineticsError::InvalidShortcutRange(_) => {
            AuditErrorKind::ReactionSyntax
        }
        KineticsError::UnknownReactionSpecies { .. } => AuditErrorKind::UnknownSpecies,
        _ => AuditErrorKind::Other,
    }
}

fn write_json(path: impl AsRef<Path>, report: &impl Serialize) -> KineticsResult<()> {
    let path = path.as_ref();
    let contents = serde_json::to_vec_pretty(report).map_err(KineticsError::from)?;
    std::fs::write(path, contents).map_err(|error| {
        KineticsError::Io(error).with_context(format!(
            "writing library audit report to `{}`",
            path.display()
        ))
    })
}

#[cfg(test)]
mod tests {
    use serde_json::json;

    use super::*;

    fn payload(eq: &str, arrhenius: Value) -> Value {
        json!({ "type": "elem", "eq": eq, "Arrenius": arrhenius })
    }

    #[test]
    fn non_transactional_audits_keep_good_records_around_middle_failures() {
        let values = [
            payload("A=>B", json!([1.0, 0.0, 0.0])),
            payload("malformed equation", json!([1.0, 0.0, 0.0])),
            payload("C=>D", json!([1.0, 0.0])),
            payload("E=>F", json!([1.0, 0.0, 0.0])),
        ];

        let records = values
            .iter()
            .enumerate()
            .map(|(index, value)| {
                let id = (index + 10).to_string();
                let result = audit_payload(value.clone(), &id);
                PayloadAuditRecord {
                    reaction_id: id,
                    range_offset: index as u64,
                    raw_payload: Some(value.clone()),
                    result,
                }
            })
            .collect();
        let payload_report = PayloadAuditReport {
            library: "fixture".into(),
            requested_reaction_ids: vec!["10".into(), "11".into(), "12".into(), "13".into()],
            requested_range: Some(10..=13),
            records,
        };

        assert_eq!(payload_report.successful_count(), 3);
        assert_eq!(payload_report.failed_count(), 1);
        assert!(!payload_report.is_clean());
        assert_eq!(payload_report.records[2].reaction_id, "12");
        assert_eq!(
            payload_report.records[2].result.as_ref().unwrap_err().stage,
            AuditStage::PayloadValidation
        );

        let equation_report = audit_parsed_equations(&payload_report);
        assert_eq!(equation_report.records.len(), 3);
        assert_eq!(equation_report.successful_count(), 2);
        assert_eq!(equation_report.failed_count(), 1);
        assert_eq!(equation_report.records[0].reaction_id, "10");
        assert_eq!(equation_report.records[1].reaction_id, "11");
        assert_eq!(
            equation_report.records[1]
                .result
                .as_ref()
                .unwrap_err()
                .stage,
            AuditStage::EquationAnalysis
        );
        assert_eq!(equation_report.records[2].reaction_id, "13");
    }

    #[test]
    fn audit_reports_serialize_as_json() {
        let record = PayloadAuditRecord {
            reaction_id: "1".into(),
            range_offset: 0,
            raw_payload: Some(payload("A=>B", json!([1.0, 0.0, 0.0]))),
            result: audit_payload(payload("A=>B", json!([1.0, 0.0, 0.0])), "1"),
        };
        let report = PayloadAuditReport {
            library: "fixture".into(),
            requested_reaction_ids: vec!["1".into()],
            requested_range: Some(1..=1),
            records: vec![record],
        };
        let json = serde_json::to_value(report).expect("audit report should serialize");
        assert_eq!(json["records"][0]["reaction_id"], "1");
        assert!(json["records"][0]["result"]["Ok"].is_object());
    }

    #[test]
    fn payload_audit_records_missing_ids_instead_of_skipping_them() {
        let record = audit_catalog_entry("fixture", "12".into(), 2, None);

        assert!(record.raw_payload.is_none());
        assert!(matches!(
            record.result,
            Err(AuditFailure {
                stage: AuditStage::CatalogLookup,
                kind: AuditErrorKind::MissingReaction,
                ..
            })
        ));
    }

    #[test]
    #[ignore = "bundled catalog audit story; included in the Kinetics campaign"]
    fn story_cantera_range_audit_writes_payload_and_equation_reports() {
        let payload_report = audit_library_payload_range("Cantera", 0..=100)
            .expect("Cantera catalog should open for a bounded audit");
        assert_eq!(payload_report.records.len(), 101);
        assert!(payload_report.is_clean());

        let equation_report = audit_parsed_equations(&payload_report);
        assert!(equation_report.is_clean());
        assert_eq!(
            equation_report.records.len(),
            payload_report.successful_count()
        );

        let output = tempfile::tempdir().expect("temporary audit report directory");
        let payload_path = output.path().join("payload-audit.json");
        let equation_path = output.path().join("equation-audit.json");
        payload_report
            .write_json(&payload_path)
            .expect("payload audit report should be persisted");
        equation_report
            .write_json(&equation_path)
            .expect("equation audit report should be persisted");

        let payload_json: Value = serde_json::from_slice(
            &std::fs::read(&payload_path).expect("payload report file should exist"),
        )
        .expect("payload report should be valid JSON");
        let equation_json: Value = serde_json::from_slice(
            &std::fs::read(&equation_path).expect("equation report file should exist"),
        )
        .expect("equation report should be valid JSON");
        assert_eq!(payload_json["records"].as_array().unwrap().len(), 101);
        assert_eq!(
            equation_json["records"].as_array().unwrap().len(),
            payload_report.successful_count()
        );

        let explicit_ids = vec!["1".to_string(), "__missing_audit_fixture_id__".to_string()];
        let explicit_report = audit_library_reactions("Cantera", &explicit_ids)
            .expect("explicit reaction IDs should be auditable");
        assert_eq!(explicit_report.records.len(), 2);
        assert_eq!(explicit_report.successful_count(), 1);
        assert_eq!(explicit_report.failed_count(), 1);
        assert!(!explicit_report.is_clean());
        assert_eq!(
            explicit_report.records[1].result.as_ref().unwrap_err().kind,
            AuditErrorKind::MissingReaction
        );

        println!(
            "Cantera audit 0..=100: payloads_ok={}, payloads_failed={}, equations_ok={}, equations_failed={}",
            payload_report.successful_count(),
            payload_report.failed_count(),
            equation_report.successful_count(),
            equation_report.failed_count()
        );
        for record in payload_report
            .records
            .iter()
            .filter_map(|record| record.result.as_ref().err().map(|error| (record, error)))
        {
            println!(
                "payload issue id={} stage={:?} kind={:?}: {}",
                record.0.reaction_id, record.1.stage, record.1.kind, record.1.message
            );
        }
        for record in equation_report
            .records
            .iter()
            .filter_map(|record| record.result.as_ref().err().map(|error| (record, error)))
        {
            println!(
                "equation issue id={} equation={:?} kind={:?}: {}",
                record.0.reaction_id, record.0.equation, record.1.kind, record.1.message
            );
        }
    }
}
