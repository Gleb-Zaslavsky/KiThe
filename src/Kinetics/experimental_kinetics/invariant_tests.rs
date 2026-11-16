use crate::Kinetics::experimental_kinetics::experiment_series_main::{TGAExperiment, TGASeries};
use crate::Kinetics::experimental_kinetics::one_experiment_dataset::{
    ColumnNature, ColumnOrigin, TGADataset, Unit,
};
use crate::Kinetics::experimental_kinetics::testing_mod::VirtualTGA;
use proptest::prelude::*;
use std::collections::HashSet;
use std::sync::Arc;

fn synthetic_dataset() -> TGADataset {
    let time: Vec<f64> = (0..12).map(|index| index as f64).collect();
    TGADataset::create_from_synthetic_data(&VirtualTGA {
        temperature: time.iter().map(|t| 300.0 + 2.0 * t).collect(),
        mass: time.iter().map(|t| 10.0 - 0.1 * t).collect(),
        time,
    })
    .expect("valid synthetic test dataset")
}

proptest! {
    #![proptest_config(ProptestConfig {
        cases: 64,
        max_shrink_iters: 1_000,
        .. ProptestConfig::default()
    })]

    #[test]
    fn generated_series_lifecycle_preserves_id_index(
        operations in prop::collection::vec((0_u8..4, any::<u8>()), 0..64)
    ) {
        let dataset = synthetic_dataset();
        let mut series = TGASeries::new();
        let mut next_name = 0_usize;

        for (operation, selector) in operations {
            match operation {
                0 => {
                    let id = format!("experiment_{next_name}");
                    next_name += 1;
                    series
                        .push(TGAExperiment::new(dataset.clone()).with_id(id))
                        .expect("generated IDs are unique");
                }
                1 if !series.is_empty() => {
                    let index = selector as usize % series.len();
                    let new_id = format!("renamed_{next_name}");
                    next_name += 1;
                    series
                        .set_experiment_id(index, &new_id)
                        .expect("generated rename target is unique");
                }
                2 if !series.is_empty() => {
                    let index = selector as usize % series.len();
                    series.drop_experiment(index).expect("valid generated index");
                }
                3 if let Some(existing_id) = series.ids().first().cloned() => {
                    let before = series.ids();
                    let duplicate = TGAExperiment::new(dataset.clone()).with_id(existing_id);
                    prop_assert!(series.push(duplicate).is_err());
                    prop_assert_eq!(series.ids(), before);
                }
                _ => {}
            }

            let ids = series.ids();
            let stored_ids: Vec<String> = series
                .experiments()
                .iter()
                .map(|experiment| experiment.meta.id.clone())
                .collect();
            prop_assert_eq!(&stored_ids, &ids);
            prop_assert_eq!(ids.iter().collect::<HashSet<_>>().len(), ids.len());

            for (index, id) in ids.iter().enumerate() {
                prop_assert_eq!(series.index_by_id(id).unwrap(), index);
                prop_assert_eq!(&series.get_experiment(index).unwrap().meta.id, id);
                prop_assert_eq!(&series.get_experiment_by_id(id).unwrap().meta.id, id);
            }
        }
    }

    #[test]
    fn generated_dataset_column_lifecycle_preserves_frame_schema_and_history(
        operations in prop::collection::vec((0_u8..4, any::<u8>()), 0..40)
    ) {
        let mut dataset = synthetic_dataset();
        let initial_rows = dataset.frame.clone().collect().unwrap().height();
        let mut generated_columns = Vec::new();
        let mut next_name = 0_usize;

        for (operation, selector) in operations {
            match operation {
                0 => {
                    let name = format!("generated_{next_name}");
                    next_name += 1;
                    let values = (0..initial_rows)
                        .map(|row| row as f64 + selector as f64 / 255.0)
                        .collect();
                    dataset = dataset
                        .add_numeric_column(
                            &name,
                            Unit::Dimensionless,
                            Arc::new(values),
                            ColumnNature::Unknown,
                        )
                        .expect("generated output name is unused");
                    generated_columns.push(name);
                }
                1 if !generated_columns.is_empty() => {
                    let old_name = generated_columns.pop().unwrap();
                    let new_name = format!("generated_{next_name}");
                    next_name += 1;
                    dataset = dataset
                        .rename_column(&old_name, &new_name)
                        .expect("generated rename target is unused");
                    generated_columns.push(new_name);
                }
                2 if !generated_columns.is_empty() => {
                    let name = generated_columns.pop().unwrap();
                    dataset = dataset.drop_column(&name).expect("generated column exists");
                }
                3 if dataset.can_undo() => {
                    dataset.undo_last().expect("snapshot is available");
                    generated_columns = dataset
                        .schema
                        .columns
                        .keys()
                        .filter(|name| name.starts_with("generated_"))
                        .cloned()
                        .collect();
                }
                _ => {}
            }

            let frame = dataset.frame.clone().collect().unwrap();
            let frame_names: HashSet<String> = frame
                .get_column_names()
                .into_iter()
                .map(|name| name.to_string())
                .collect();
            let schema_names: HashSet<String> = dataset.schema.columns.keys().cloned().collect();
            prop_assert_eq!(&frame_names, &schema_names);
            prop_assert_eq!(frame.height(), initial_rows);
            for name in frame.get_column_names() {
                prop_assert_eq!(frame.column(name.as_str()).unwrap().len(), initial_rows);
            }
            let changed_names: HashSet<&str> = dataset
                .history_of_operations
                .columns_has_changed
                .iter()
                .flatten()
                .map(String::as_str)
                .collect();
            prop_assert_eq!(
                changed_names.len(),
                dataset.history_of_operations.columns_has_changed.len()
            );
            for (index, record) in dataset
                .history_of_operations
                .vector_of_operations
                .iter()
                .enumerate()
            {
                prop_assert_eq!(record.timestamp, index);
            }

            for (name, metadata) in &dataset.schema.columns {
                prop_assert_eq!(&metadata.name, name);
                let provenance_matches_name = metadata
                    .provenance
                    .steps
                    .last()
                    .is_some_and(|step| step.output_column == *name);
                prop_assert!(provenance_matches_name);
            }

            for binding in [
                dataset.schema.time.as_ref(),
                dataset.schema.temperature.as_ref(),
                dataset.schema.mass.as_ref(),
                dataset.schema.alpha.as_ref(),
                dataset.schema.dm_dt.as_ref(),
                dataset.schema.eta.as_ref(),
                dataset.schema.deta_dt.as_ref(),
                dataset.schema.dalpha_dt.as_ref(),
                dataset.schema.dT_dt.as_ref(),
                dataset.schema.E.as_ref(),
                dataset.schema.R2.as_ref(),
            ]
            .into_iter()
            .flatten()
            {
                prop_assert!(frame_names.contains(binding));
                prop_assert!(schema_names.contains(binding));
            }
        }
    }

    #[test]
    fn generated_time_sort_is_strict_and_keeps_rows_aligned(
        increments in prop::collection::vec(1_u16..1000, 2..40),
        rotation in any::<usize>(),
    ) {
        let mut time = Vec::with_capacity(increments.len() + 1);
        time.push(0.0);
        for increment in increments {
            time.push(time.last().unwrap() + f64::from(increment));
        }
        let mut mass: Vec<f64> = time.iter().map(|value| 10_000.0 - value).collect();
        let mut temperature: Vec<f64> = time.iter().map(|value| 300.0 + value).collect();
        let shift = rotation % (time.len() - 1) + 1;
        time.rotate_left(shift);
        mass.rotate_left(shift);
        temperature.rotate_left(shift);

        let mut dataset = TGADataset::create_from_synthetic_data(&VirtualTGA {
            time,
            temperature,
            mass,
        }).expect("generated finite equal-length columns");
        prop_assert!(!dataset.time_monotonicity_report().unwrap().is_strictly_increasing);

        dataset.sort_by_bound_time_inplace().unwrap();
        prop_assert!(dataset.time_monotonicity_report().unwrap().is_strictly_increasing);
        let block = dataset.materialize("time", &["mass", "temperature"]).unwrap();
        prop_assert_eq!(block.grid.len(), block.columns["mass"].len());
        prop_assert_eq!(block.grid.len(), block.columns["temperature"].len());
        for ((time, mass), temperature) in block.grid.iter()
            .zip(block.columns["mass"].iter())
            .zip(block.columns["temperature"].iter())
        {
            prop_assert!((mass - (10_000.0 - time)).abs() < 1e-9);
            prop_assert!((temperature - (300.0 + time)).abs() < 1e-9);
        }
    }

    #[test]
    fn generated_series_csv_roundtrip_preserves_metadata_and_free_text(
        id in "[a-zA-Z0-9_-]{1,16}",
        comment in any::<String>(),
        unit_index in 0_u8..4,
        nature_index in 0_u8..6,
        origin_index in 0_u8..4,
    ) {
        let mut dataset = synthetic_dataset()
            .add_numeric_column(
                "auxiliary",
                Unit::Dimensionless,
                Arc::new((0..12).map(|row| row as f64).collect()),
                ColumnNature::Unknown,
            )
            .unwrap();
        let units = [Unit::Unknown, Unit::Dimensionless, Unit::Minute, Unit::Kelvin];
        let natures = [
            ColumnNature::Unknown,
            ColumnNature::Time,
            ColumnNature::Temperature,
            ColumnNature::Mass,
            ColumnNature::Conversion,
            ColumnNature::MassRate,
        ];
        let origins = [
            ColumnOrigin::Raw,
            ColumnOrigin::PolarsDerived,
            ColumnOrigin::NumericDerived,
            ColumnOrigin::Imported,
        ];
        let metadata = dataset.schema.columns.get_mut("auxiliary").unwrap();
        metadata.unit = units[unit_index as usize];
        metadata.nature = natures[nature_index as usize];
        metadata.origin = origins[origin_index as usize];

        let expected_unit = metadata.unit;
        let expected_nature = metadata.nature;
        let expected_origin = metadata.origin;
        let expected_roots = metadata.provenance.root_columns.clone();
        let expected_steps = metadata.provenance.steps.len();
        let mut experiment = TGAExperiment::new(dataset).with_id(id.clone());
        experiment.meta.comment = Some(comment.clone());
        let mut series = TGASeries::new();
        series.push(experiment).unwrap();

        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("roundtrip.csv");
        series.to_csv_series(&path).unwrap();
        let restored = TGASeries::from_csv_series(&path).unwrap();
        let experiment = restored.get_experiment_by_id(&id).unwrap();
        let restored_metadata = &experiment.dataset.schema.columns["auxiliary"];
        prop_assert_eq!(experiment.meta.comment.as_deref(), Some(comment.as_str()));
        prop_assert_eq!(restored_metadata.unit, expected_unit);
        prop_assert_eq!(restored_metadata.nature, expected_nature);
        prop_assert_eq!(restored_metadata.origin, expected_origin);
        prop_assert_eq!(&restored_metadata.provenance.root_columns, &expected_roots);
        prop_assert_eq!(restored_metadata.provenance.steps.len(), expected_steps);
    }
}
