#[cfg(test)]
mod tests {
    use super::super::frozen_reference::{
        FrozenReferenceDataset, FrozenReferenceError, FrozenReferenceEvidenceKind,
        NistBenzeneTolueneVleReference, TemperaturePressureReference,
    };
    use serde_json::{Value, json};
    use std::fs;
    use std::path::PathBuf;

    const DATA_FILE: &str = "temperature_pressure.rows.json";

    fn fixture_dir() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("src")
            .join("Thermodynamics")
            .join("ChemEquilibrium")
            .join("frozen_reference")
            .join("data")
            .join("synthetic")
    }

    fn valid_metadata() -> Value {
        json!({
            "dataset_format_version": 1,
            "dataset_id": "synthetic.temperature_pressure.v1",
            "title": "Synthetic temperature-pressure loader fixture",
            "evidence_kind": "synthetic_infrastructure",
            "data_file": DATA_FILE,
            "source": {
                "organization": "KiThe test suite",
                "name": "synthetic infrastructure fixture",
                "version": "1",
                "citation": "Not physical evidence; hand-authored for parser tests",
                "source_table": "three-row example",
                "stable_identifier": "kithe.synthetic.temperature-pressure.v1"
            },
            "transcription": "Hand-authored values used only to test frozen-reference infrastructure.",
            "quantities": [
                {
                    "column": "temperature_k",
                    "meaning": "absolute temperature",
                    "unit": "K",
                    "source_precision": "1 K"
                },
                {
                    "column": "pressure_pa",
                    "meaning": "absolute pressure",
                    "unit": "Pa",
                    "uncertainty": { "value": 0.0, "unit": "Pa" }
                }
            ]
        })
    }

    fn valid_rows() -> Value {
        json!({
            "dataset_id": "synthetic.temperature_pressure.v1",
            "rows": [
                { "temperature_k": 300.0, "pressure_pa": 1000.0 },
                { "temperature_k": 400.0, "pressure_pa": 2000.0 },
                { "temperature_k": 500.0, "pressure_pa": 3000.0 }
            ]
        })
    }

    fn parse(
        metadata: &Value,
        rows: &Value,
    ) -> Result<FrozenReferenceDataset<TemperaturePressureReference>, FrozenReferenceError> {
        FrozenReferenceDataset::from_json_strs(
            &serde_json::to_string(metadata).unwrap(),
            &serde_json::to_string(rows).unwrap(),
            DATA_FILE,
        )
    }

    #[test]
    fn synthetic_fixture_loads_with_visible_provenance_but_is_not_i5_evidence() {
        let directory = fixture_dir();
        let dataset = FrozenReferenceDataset::<TemperaturePressureReference>::load(
            directory.join("temperature_pressure.metadata.json"),
            directory.join(DATA_FILE),
        )
        .unwrap();

        assert_eq!(dataset.rows().len(), 3);
        assert_eq!(dataset.rows()[1].temperature_k, 400.0);
        assert_eq!(dataset.rows()[1].pressure_pa, 2000.0);
        assert_eq!(
            dataset.metadata().evidence_kind,
            FrozenReferenceEvidenceKind::SyntheticInfrastructure
        );
        assert!(!dataset.is_external_evidence());
        let context = dataset.provenance_context(1);
        assert!(context.contains("synthetic.temperature_pressure.v1"));
        assert!(context.contains("row=1"));
        assert!(context.contains("synthetic infrastructure fixture"));
    }

    #[test]
    fn nist_benzene_toluene_vle_subset_is_external_read_only_px_evidence() {
        let directory = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("src")
            .join("Thermodynamics")
            .join("ChemEquilibrium")
            .join("frozen_reference")
            .join("data")
            .join("nist_thermoml");
        let metadata_path = directory.join("benzene_toluene_vle_353_15k.metadata.json");
        let rows_path = directory.join("benzene_toluene_vle_353_15k.rows.json");
        let metadata_before = fs::read(&metadata_path).unwrap();
        let rows_before = fs::read(&rows_path).unwrap();

        let dataset = FrozenReferenceDataset::<NistBenzeneTolueneVleReference>::load(
            &metadata_path,
            &rows_path,
        )
        .expect("reviewed NIST ThermoML benzene/toluene P-x subset must load");

        assert!(dataset.is_external_evidence());
        assert_eq!(dataset.rows().len(), 5);
        assert_eq!(
            dataset.rows().first().unwrap().liquid_benzene_mole_fraction,
            0.0
        );
        assert_eq!(
            dataset.rows().last().unwrap().liquid_benzene_mole_fraction,
            1.0
        );
        assert!(
            dataset
                .provenance_context(2)
                .contains("nist.thermoml.benzene_toluene.vle.353_15k.v1")
        );
        assert_eq!(fs::read(metadata_path).unwrap(), metadata_before);
        assert_eq!(fs::read(rows_path).unwrap(), rows_before);
    }

    #[test]
    fn loading_is_byte_for_byte_read_only() {
        let directory = fixture_dir();
        let metadata_path = directory.join("temperature_pressure.metadata.json");
        let data_path = directory.join(DATA_FILE);
        let metadata_before = fs::read(&metadata_path).unwrap();
        let data_before = fs::read(&data_path).unwrap();

        FrozenReferenceDataset::<TemperaturePressureReference>::load(&metadata_path, &data_path)
            .unwrap();

        assert_eq!(fs::read(metadata_path).unwrap(), metadata_before);
        assert_eq!(fs::read(data_path).unwrap(), data_before);
    }

    #[test]
    fn missing_file_is_reported_as_a_typed_read_error() {
        let directory = fixture_dir();
        let error = FrozenReferenceDataset::<TemperaturePressureReference>::load(
            directory.join("missing.metadata.json"),
            directory.join(DATA_FILE),
        )
        .unwrap_err();
        assert!(matches!(error, FrozenReferenceError::Read { .. }));
    }

    #[test]
    fn row_validation_rejects_empty_duplicate_unordered_and_nonpositive_data() {
        let cases = [
            (json!([]), "at least one row"),
            (
                json!([
                    { "temperature_k": 300.0, "pressure_pa": 1000.0 },
                    { "temperature_k": 300.0, "pressure_pa": 2000.0 }
                ]),
                "duplicate temperature",
            ),
            (
                json!([
                    { "temperature_k": 400.0, "pressure_pa": 1000.0 },
                    { "temperature_k": 300.0, "pressure_pa": 2000.0 }
                ]),
                "strictly increasing",
            ),
            (
                json!([{ "temperature_k": 300.0, "pressure_pa": 0.0 }]),
                "finite and positive",
            ),
        ];

        for (rows_value, expected) in cases {
            let mut rows = valid_rows();
            rows["rows"] = rows_value;
            let error = parse(&valid_metadata(), &rows).unwrap_err();
            assert!(error.to_string().contains(expected), "{error}");
        }
    }

    #[test]
    fn metadata_validation_rejects_empty_fields_duplicate_columns_and_bad_uncertainty() {
        let mut missing_format_version = valid_metadata();
        missing_format_version
            .as_object_mut()
            .unwrap()
            .remove("dataset_format_version");
        assert!(matches!(
            parse(&missing_format_version, &valid_rows()),
            Err(FrozenReferenceError::Parse { .. })
        ));

        let mut zero_format_version = valid_metadata();
        zero_format_version["dataset_format_version"] = json!(0);
        assert!(
            parse(&zero_format_version, &valid_rows())
                .unwrap_err()
                .to_string()
                .contains("unsupported format version 0")
        );

        let mut future_format_version = valid_metadata();
        future_format_version["dataset_format_version"] = json!(2);
        assert!(
            parse(&future_format_version, &valid_rows())
                .unwrap_err()
                .to_string()
                .contains("unsupported format version 2")
        );

        let mut empty_title = valid_metadata();
        empty_title["title"] = json!("  ");
        assert!(matches!(
            parse(&empty_title, &valid_rows()).unwrap_err(),
            FrozenReferenceError::InvalidMetadata { ref field, .. } if field == "title"
        ));

        let mut duplicate = valid_metadata();
        duplicate["quantities"][1]["column"] = json!("temperature_k");
        assert!(matches!(
            parse(&duplicate, &valid_rows()).unwrap_err(),
            FrozenReferenceError::InvalidMetadata { ref field, .. }
                if field == "quantities.column"
        ));

        let mut uncertainty = valid_metadata();
        uncertainty["quantities"][1]["uncertainty"]["value"] = json!(-1.0);
        assert!(matches!(
            parse(&uncertainty, &valid_rows()).unwrap_err(),
            FrozenReferenceError::InvalidMetadata { ref field, .. }
                if field == "quantities.uncertainty.value"
        ));
    }

    #[test]
    fn schema_validation_rejects_missing_unknown_and_mismatched_units() {
        let mut wrong_unit = valid_metadata();
        wrong_unit["quantities"][1]["unit"] = json!("kPa");
        assert!(matches!(
            parse(&wrong_unit, &valid_rows()).unwrap_err(),
            FrozenReferenceError::SchemaMismatch { ref column, .. }
                if column == "pressure_pa"
        ));

        let mut unknown_column = valid_metadata();
        unknown_column["quantities"][1]["column"] = json!("pressure_bar");
        assert!(matches!(
            parse(&unknown_column, &valid_rows()).unwrap_err(),
            FrozenReferenceError::SchemaMismatch { ref column, .. }
                if column == "pressure_pa"
        ));

        let mut missing_column = valid_metadata();
        missing_column["quantities"].as_array_mut().unwrap().pop();
        assert!(matches!(
            parse(&missing_column, &valid_rows()).unwrap_err(),
            FrozenReferenceError::SchemaMismatch { .. }
        ));
    }

    #[test]
    fn pair_identity_rejects_dataset_id_and_filename_mismatches() {
        let mut wrong_id = valid_rows();
        wrong_id["dataset_id"] = json!("another.dataset");
        assert!(matches!(
            parse(&valid_metadata(), &wrong_id).unwrap_err(),
            FrozenReferenceError::DatasetIdMismatch { .. }
        ));

        let mut wrong_file = valid_metadata();
        wrong_file["data_file"] = json!("other.rows.json");
        assert!(matches!(
            parse(&wrong_file, &valid_rows()).unwrap_err(),
            FrozenReferenceError::DataFileMismatch { .. }
        ));
    }

    #[test]
    fn strict_parser_rejects_missing_extra_malformed_and_nonfinite_rows() {
        let metadata_json = serde_json::to_string(&valid_metadata()).unwrap();
        let malformed_cases = [
            r#"{"dataset_id":"synthetic.temperature_pressure.v1","rows":[{"temperature_k":300.0}]}"#,
            r#"{"dataset_id":"synthetic.temperature_pressure.v1","rows":[{"temperature_k":300.0,"pressure_pa":1000.0,"extra":1.0}]}"#,
            r#"{"dataset_id":"synthetic.temperature_pressure.v1","rows":[{"temperature_k":"NaN","pressure_pa":1000.0}]}"#,
            r#"{"dataset_id":"synthetic.temperature_pressure.v1","rows":[{"temperature_k":1e999,"pressure_pa":1000.0}]}"#,
        ];

        for rows_json in malformed_cases {
            let error = FrozenReferenceDataset::<TemperaturePressureReference>::from_json_strs(
                &metadata_json,
                rows_json,
                DATA_FILE,
            )
            .unwrap_err();
            assert!(matches!(error, FrozenReferenceError::Parse { .. }));
        }

        let mut missing_metadata = valid_metadata();
        missing_metadata["source"]
            .as_object_mut()
            .unwrap()
            .remove("citation");
        assert!(matches!(
            parse(&missing_metadata, &valid_rows()).unwrap_err(),
            FrozenReferenceError::Parse { .. }
        ));

        let mut extra_metadata = valid_metadata();
        extra_metadata["unreviewed_note"] = json!("must not be silently ignored");
        assert!(matches!(
            parse(&extra_metadata, &valid_rows()).unwrap_err(),
            FrozenReferenceError::Parse { .. }
        ));
    }
}
