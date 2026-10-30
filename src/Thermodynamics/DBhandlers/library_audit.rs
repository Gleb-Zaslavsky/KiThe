//! Record-by-record diagnostics for thermodynamic and transport libraries.
//!
//! Production handlers remain fail-fast. This audit API is intentionally
//! separate: every selected record gets an outcome, so one damaged row does not
//! hide later records. Reports do not mutate or populate production state.
//!
//! ```no_run
//! use KiThe::Thermodynamics::DBhandlers::library_audit::{
//!     audit_entire_library_with_telemetry, DatabaseAuditApiError,
//! };
//! use KiThe::Thermodynamics::thermo_lib_api::ThermoData;
//!
//! fn inspect_nasa() -> Result<(), DatabaseAuditApiError> {
//!     let mut database = ThermoData::open()?;
//!     let (report, timings) = audit_entire_library_with_telemetry(&database, "NASA_gas")?;
//!     report.write_json("nasa-gas-audit.json")?;
//!     println!("clean={} timings={timings:?}", report.is_clean());
//!     database.close();
//!     Ok(())
//! }
//! ```

use std::path::Path;
use std::time::{Duration, Instant};

use serde::{Deserialize, Serialize};
use serde_json::Value;
use thiserror::Error;

use crate::Thermodynamics::DBhandlers::thermo_api::{
    ThermoCalculator, ThermoError, create_thermal_by_name,
};
use crate::Thermodynamics::DBhandlers::transport_api::{
    TransportCalculator, TransportError, create_transport_calculator_by_name,
};
use crate::Thermodynamics::thermo_lib_api::{LibraryCapability, ThermoData, ThermoLibraryError};

#[derive(Debug, Error)]
pub enum DatabaseAuditApiError {
    #[error(transparent)]
    Library(#[from] ThermoLibraryError),
    #[error("no thermo/transport handler is registered for library `{library}`")]
    UnsupportedLibrary { library: String },
    #[error("library `{library}` is not present in the open repository snapshot")]
    LibraryNotOpen { library: String },
    #[error("failed to serialize database audit report: {0}")]
    Serialize(#[from] serde_json::Error),
    #[error("failed to write database audit report to `{path}`: {source}")]
    WriteReport {
        path: String,
        #[source]
        source: std::io::Error,
    },
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum DatabaseAuditStage {
    RecordLookup,
    PayloadDeserialize,
    CoefficientParse,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum DatabaseAuditErrorKind {
    MissingRecord,
    Json,
    Schema,
    CoefficientParse,
    Other,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct DatabaseAuditFailure {
    pub library: String,
    pub record_key: String,
    pub stage: DatabaseAuditStage,
    pub kind: DatabaseAuditErrorKind,
    pub message: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DatabaseAuditRecord {
    pub record_key: String,
    pub raw_payload: Option<Value>,
    pub result: Result<(), DatabaseAuditFailure>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DatabaseAuditReport {
    pub library: String,
    pub records: Vec<DatabaseAuditRecord>,
}

impl DatabaseAuditReport {
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

    pub fn write_json(&self, path: impl AsRef<Path>) -> Result<(), DatabaseAuditApiError> {
        let path = path.as_ref();
        let bytes = serde_json::to_vec_pretty(self)?;
        std::fs::write(path, bytes).map_err(|source| DatabaseAuditApiError::WriteReport {
            path: path.display().to_string(),
            source,
        })
    }
}

#[derive(Debug, Clone, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct DatabaseAuditTelemetry {
    pub total: Duration,
    pub catalog_selection: Duration,
    pub record_lookup: Duration,
    pub payload_deserialize: Duration,
    pub coefficient_parse: Duration,
    pub records_requested: usize,
    pub records_found: usize,
    pub records_valid: usize,
    pub records_failed: usize,
}

/// Inspect explicitly selected record keys and continue after every row error.
pub fn audit_library_records(
    database: &ThermoData,
    library: &str,
    record_keys: &[String],
) -> Result<DatabaseAuditReport, DatabaseAuditApiError> {
    audit_library_records_inner(database, library, record_keys, None).map(|(report, _)| report)
}

/// Inspect explicitly selected record keys and return opt-in stage timings.
pub fn audit_library_records_with_telemetry(
    database: &ThermoData,
    library: &str,
    record_keys: &[String],
) -> Result<(DatabaseAuditReport, DatabaseAuditTelemetry), DatabaseAuditApiError> {
    let (report, telemetry) =
        audit_library_records_inner(database, library, record_keys, Some(Instant::now()))?;
    Ok((report, telemetry.expect("telemetry requested")))
}

/// Audit every record in one already-open library snapshot.
pub fn audit_entire_library(
    database: &ThermoData,
    library: &str,
) -> Result<DatabaseAuditReport, DatabaseAuditApiError> {
    audit_entire_library_inner(database, library, false).map(|(report, _)| report)
}

/// Audit every record in one already-open library snapshot with stage timings.
pub fn audit_entire_library_with_telemetry(
    database: &ThermoData,
    library: &str,
) -> Result<(DatabaseAuditReport, DatabaseAuditTelemetry), DatabaseAuditApiError> {
    audit_entire_library_inner(database, library, true)
        .map(|(report, telemetry)| (report, telemetry.expect("telemetry requested")))
}

fn audit_entire_library_inner(
    database: &ThermoData,
    library: &str,
    telemetry: bool,
) -> Result<(DatabaseAuditReport, Option<DatabaseAuditTelemetry>), DatabaseAuditApiError> {
    let canonical = ThermoData::canonical_library_name(library);
    let records = database.LibThermoData.get(&canonical).ok_or_else(|| {
        DatabaseAuditApiError::LibraryNotOpen {
            library: canonical.clone(),
        }
    })?;
    let mut keys: Vec<_> = records.keys().cloned().collect();
    keys.sort();
    audit_library_records_inner(database, &canonical, &keys, telemetry.then(Instant::now))
}

fn audit_library_records_inner(
    database: &ThermoData,
    library: &str,
    record_keys: &[String],
    started: Option<Instant>,
) -> Result<(DatabaseAuditReport, Option<DatabaseAuditTelemetry>), DatabaseAuditApiError> {
    let selection_started = started.map(|_| Instant::now());
    let canonical = ThermoData::canonical_library_name(library);
    let capability = ThermoData::library_capability(&canonical).ok_or_else(|| {
        DatabaseAuditApiError::UnsupportedLibrary {
            library: canonical.clone(),
        }
    })?;
    let records = database.LibThermoData.get(&canonical).ok_or_else(|| {
        DatabaseAuditApiError::LibraryNotOpen {
            library: canonical.clone(),
        }
    })?;

    let mut telemetry = started.map(|_| DatabaseAuditTelemetry {
        records_requested: record_keys.len(),
        ..DatabaseAuditTelemetry::default()
    });
    if let (Some(stats), Some(started)) = (telemetry.as_mut(), selection_started) {
        stats.catalog_selection = started.elapsed();
    }
    let mut audited = Vec::with_capacity(record_keys.len());
    for key in record_keys {
        let lookup_started = telemetry.as_ref().map(|_| Instant::now());
        let raw_payload = records.get(key).cloned();
        if let Some(stats) = telemetry.as_mut() {
            stats.record_lookup += lookup_started.unwrap().elapsed();
            stats.records_found += usize::from(raw_payload.is_some());
        }

        let result = match raw_payload.as_ref() {
            None => Err(failure(
                &canonical,
                key,
                DatabaseAuditStage::RecordLookup,
                DatabaseAuditErrorKind::MissingRecord,
                "record key is absent from the selected library",
            )),
            Some(payload) => {
                audit_payload(&canonical, key, payload, capability, telemetry.as_mut())
            }
        };
        if let Some(stats) = telemetry.as_mut() {
            if result.is_ok() {
                stats.records_valid += 1;
            } else {
                stats.records_failed += 1;
            }
        }
        audited.push(DatabaseAuditRecord {
            record_key: key.clone(),
            raw_payload,
            result,
        });
    }

    let report = DatabaseAuditReport {
        library: canonical,
        records: audited,
    };
    if let Some(stats) = telemetry.as_mut() {
        stats.total = started.unwrap().elapsed();
    }
    Ok((report, telemetry))
}

fn audit_payload(
    library: &str,
    record_key: &str,
    payload: &Value,
    capability: LibraryCapability,
    mut telemetry: Option<&mut DatabaseAuditTelemetry>,
) -> Result<(), DatabaseAuditFailure> {
    match capability {
        LibraryCapability::Thermo => {
            let deserialize_started = telemetry.as_ref().map(|_| Instant::now());
            let mut calculator = create_thermal_by_name(library);
            let deserialize_result = calculator.from_serde(payload.clone());
            if let (Some(stats), Some(started)) = (telemetry.as_deref_mut(), deserialize_started) {
                stats.payload_deserialize += started.elapsed();
            }
            deserialize_result.map_err(|error| {
                failure_from_thermo(
                    library,
                    record_key,
                    DatabaseAuditStage::PayloadDeserialize,
                    error,
                )
            })?;

            let parse_started = deserialize_started.map(|_| Instant::now());
            let parse_result = calculator.parse_coefficients();
            if let (Some(stats), Some(started)) = (telemetry.as_deref_mut(), parse_started) {
                stats.coefficient_parse += started.elapsed();
            }
            parse_result.map_err(|error| {
                failure_from_thermo(
                    library,
                    record_key,
                    DatabaseAuditStage::CoefficientParse,
                    error,
                )
            })
        }
        LibraryCapability::Transport => {
            let deserialize_started = telemetry.as_ref().map(|_| Instant::now());
            let mut calculator = create_transport_calculator_by_name(library).map_err(|error| {
                failure_from_transport(
                    library,
                    record_key,
                    DatabaseAuditStage::PayloadDeserialize,
                    error,
                )
            })?;
            let result = calculator.from_serde(payload.clone());
            if let (Some(stats), Some(started)) = (telemetry.as_deref_mut(), deserialize_started) {
                stats.payload_deserialize += started.elapsed();
            }
            result.map_err(|error| {
                failure_from_transport(
                    library,
                    record_key,
                    DatabaseAuditStage::PayloadDeserialize,
                    error,
                )
            })
        }
    }
}

fn failure_from_thermo(
    library: &str,
    record_key: &str,
    stage: DatabaseAuditStage,
    error: ThermoError,
) -> DatabaseAuditFailure {
    let kind = match &error {
        ThermoError::SerdeError(_) => DatabaseAuditErrorKind::Json,
        ThermoError::NISTSchema(_) | ThermoError::InvalidTemperatureRange => {
            DatabaseAuditErrorKind::Schema
        }
        ThermoError::NoCoefficientsFound { .. } | ThermoError::NISTFit(_) => {
            DatabaseAuditErrorKind::CoefficientParse
        }
        _ => DatabaseAuditErrorKind::Other,
    };
    failure(library, record_key, stage, kind, error.to_string())
}

fn failure_from_transport(
    library: &str,
    record_key: &str,
    stage: DatabaseAuditStage,
    error: TransportError,
) -> DatabaseAuditFailure {
    let kind = match &error {
        TransportError::SerdeError(_) => DatabaseAuditErrorKind::Json,
        TransportError::MissingData(_) | TransportError::MissingCoefficients(_) => {
            DatabaseAuditErrorKind::Schema
        }
        TransportError::ParseError(_) => DatabaseAuditErrorKind::CoefficientParse,
        _ => DatabaseAuditErrorKind::Other,
    };
    failure(library, record_key, stage, kind, error.to_string())
}

fn failure(
    library: &str,
    record_key: &str,
    stage: DatabaseAuditStage,
    kind: DatabaseAuditErrorKind,
    message: impl Into<String>,
) -> DatabaseAuditFailure {
    DatabaseAuditFailure {
        library: library.to_string(),
        record_key: record_key.to_string(),
        stage,
        kind,
        message: message.into(),
    }
}

#[cfg(test)]
mod tests {
    use std::collections::HashMap;
    use std::sync::Arc;

    use super::*;

    #[test]
    fn record_audit_keeps_location_and_continues_after_middle_payload_failure() {
        let mut database = ThermoData::try_new().expect("bundled thermo repository loads");
        let valid_payload = database.LibThermoData["NASA_gas"]["CO"].clone();
        let records = HashMap::from([
            ("first".to_string(), valid_payload.clone()),
            (
                "broken".to_string(),
                serde_json::json!({ "malformed": true }),
            ),
            ("last".to_string(), valid_payload),
        ]);
        Arc::make_mut(&mut database.LibThermoData).insert("NASA_gas".into(), records);

        let keys = vec!["first".into(), "broken".into(), "last".into()];
        let (report, telemetry) =
            audit_library_records_with_telemetry(&database, "NASA_gas", &keys)
                .expect("registered database library should audit");

        assert_eq!(report.records.len(), 3);
        assert_eq!(report.successful_count(), 2);
        assert_eq!(report.failed_count(), 1);
        assert_eq!(report.records[1].record_key, "broken");
        let failure = report.records[1].result.as_ref().unwrap_err();
        assert_eq!(failure.library, "NASA_gas");
        assert_eq!(failure.record_key, "broken");
        assert_eq!(failure.stage, DatabaseAuditStage::PayloadDeserialize);
        assert_eq!(failure.kind, DatabaseAuditErrorKind::Json);
        assert!(report.records[1].raw_payload.is_some());

        assert_eq!(telemetry.records_requested, 3);
        assert_eq!(telemetry.records_found, 3);
        assert_eq!(telemetry.records_valid, 2);
        assert_eq!(telemetry.records_failed, 1);
        assert!(telemetry.payload_deserialize > Duration::ZERO);
        assert!(telemetry.coefficient_parse > Duration::ZERO);
    }

    #[test]
    fn record_audit_reports_missing_keys_without_aborting() {
        let database = ThermoData::try_new().expect("bundled thermo repository loads");
        let keys = vec!["CO".into(), "__missing_record__".into()];
        let report = audit_library_records(&database, "NASA_gas", &keys)
            .expect("registered database library should audit");

        assert_eq!(report.successful_count(), 1);
        let failure = report.records[1].result.as_ref().unwrap_err();
        assert_eq!(failure.stage, DatabaseAuditStage::RecordLookup);
        assert_eq!(failure.kind, DatabaseAuditErrorKind::MissingRecord);
        assert_eq!(failure.library, "NASA_gas");
    }

    #[test]
    #[ignore = "full NASA gas catalog stress story; included in DBhandlers campaign"]
    fn story_nasa_gas_library_audit_scans_and_persists_every_record() {
        let database = ThermoData::try_new().expect("bundled thermo repository loads");
        let total = database.LibThermoData["NASA_gas"].len();
        let (report, telemetry) = audit_entire_library_with_telemetry(&database, "NASA_gas")
            .expect("NASA gas library should be auditable record by record");
        assert_eq!(report.records.len(), total);
        assert_eq!(
            telemetry.records_requested, total,
            "telemetry accounts for every catalog row"
        );
        assert_eq!(telemetry.records_valid + telemetry.records_failed, total);

        let output = tempfile::tempdir().expect("temporary database audit report directory");
        let path = output.path().join("nasa-gas-audit.json");
        report
            .write_json(&path)
            .expect("full audit report should persist as JSON");
        let persisted: DatabaseAuditReport =
            serde_json::from_slice(&std::fs::read(path).expect("audit output should exist"))
                .expect("persisted audit should deserialize");
        assert_eq!(persisted.records.len(), total);

        println!(
            "NASA_gas full audit: records={} valid={} failed={} total_ms={:.3} lookup_ms={:.3} deserialize_ms={:.3} coefficient_parse_ms={:.3}",
            total,
            telemetry.records_valid,
            telemetry.records_failed,
            telemetry.total.as_secs_f64() * 1e3,
            telemetry.record_lookup.as_secs_f64() * 1e3,
            telemetry.payload_deserialize.as_secs_f64() * 1e3,
            telemetry.coefficient_parse.as_secs_f64() * 1e3,
        );
    }

    #[test]
    #[ignore = "remaining full-catalog stress stories; included in DBhandlers campaign"]
    fn story_remaining_supported_catalog_families_are_auditable_without_aborting() {
        let database = ThermoData::try_new().expect("bundled thermo repository loads");
        let output = tempfile::tempdir().expect("temporary audit report directory");
        let mut available_libraries: Vec<_> = database.LibThermoData.keys().cloned().collect();
        available_libraries.sort();
        for library in ["nuig_thermo", "NASA_cond", "CEA", "Aramco_transport"] {
            let expected = database
                .LibThermoData
                .get(library)
                .unwrap_or_else(|| panic!("missing {library}; available={available_libraries:?}"))
                .len();
            let (report, telemetry) = audit_entire_library_with_telemetry(&database, library)
                .unwrap_or_else(|error| panic!("failed to start audit for {library}: {error}"));
            assert_eq!(report.records.len(), expected, "library {library}");
            assert_eq!(
                telemetry.records_valid + telemetry.records_failed,
                expected,
                "telemetry should account for all {library} rows"
            );

            let path = output.path().join(format!("{library}-audit.json"));
            report
                .write_json(&path)
                .unwrap_or_else(|error| panic!("failed to persist {library} report: {error}"));
            println!(
                "{library}: records={expected} valid={} failed={} total_ms={:.3} decode_ms={:.3} coefficient_parse_ms={:.3}",
                telemetry.records_valid,
                telemetry.records_failed,
                telemetry.total.as_secs_f64() * 1e3,
                telemetry.payload_deserialize.as_secs_f64() * 1e3,
                telemetry.coefficient_parse.as_secs_f64() * 1e3,
            );
            for record in report
                .records
                .iter()
                .filter(|record| record.result.is_err())
            {
                let failure = record.result.as_ref().expect_err("filtered failed record");
                println!(
                    "{library} invalid record `{}` at {:?} ({:?}): {}",
                    failure.record_key, failure.stage, failure.kind, failure.message
                );
            }
        }
    }
}
