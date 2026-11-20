//! MODULE FOR I/O OPERATIONS WITH KINETIC METHODS
use super::kinetic_methods::KineticDataView;
use crate::Kinetics::experimental_kinetics::experiment_series_main::{TGAExperiment, TGASeries};
use crate::Kinetics::experimental_kinetics::experiment_series2::UnitedDataset;
use crate::Kinetics::experimental_kinetics::kinetic_methods::integral_isoconversion::IsoconversionalResult;
use crate::Kinetics::experimental_kinetics::kinetic_methods::isoconversion::IsoconversionalMethod;
use crate::Kinetics::experimental_kinetics::one_experiment_dataset::{
    ColumnNature, TGADomainError, Unit,
};
use log::info;
use polars::prelude::*;
impl TGASeries {
    /// Resolve experiment references for union operations.
    ///
    /// If `what_exp_to_take` is `None`, all experiments in the series are used.
    /// If `Some(&[id...])`, only the listed ids are used (in the same order).
    fn select_experiments<'a>(
        &'a self,
        what_exp_to_take: Option<&[&str]>,
    ) -> Result<Vec<&'a TGAExperiment>, TGADomainError> {
        match what_exp_to_take {
            Some(ids) => ids.iter().map(|id| self.get_experiment_by_id(id)).collect(),
            None => Ok(self.experiments().iter().collect()),
        }
    }

    /// Build a united dataset in "vertical stack" layout.
    ///
    /// Source data:
    /// - each row comes from original `TGADataset` rows in selected experiments
    /// - only columns with requested `ColumnNature` are taken
    ///
    /// Output layout:
    /// - rows from all experiments are appended one under another
    /// - every row gets `exp_id` showing experiment ownership
    /// - selected columns are normalized to canonical names
    ///
    /// This is a long/tidy table representation suitable for filters,
    /// plotting pipelines, and group-by operations over `exp_id`.
    pub fn concat_into_vertical_stack(
        &self,
        what_exp_to_take: Option<&[&str]>,
        what_cols_take: Vec<ColumnNature>,
    ) -> Result<UnitedDataset, TGADomainError> {
        if what_cols_take.is_empty() {
            return Err(TGADomainError::InvalidOperation(
                "No column natures were requested".to_string(),
            ));
        }

        let experiments = self.select_experiments(what_exp_to_take)?;
        if experiments.is_empty() {
            return Ok(UnitedDataset::empty());
        }

        let mut frames = Vec::with_capacity(experiments.len());
        let mut meta = Vec::with_capacity(experiments.len());
        let mut units = std::collections::HashMap::with_capacity(experiments.len());
        let mut common_units: Vec<(ColumnNature, Unit)> = Vec::with_capacity(what_cols_take.len());

        for exp in experiments {
            let mut select_exprs = Vec::with_capacity(what_cols_take.len());

            for nature in &what_cols_take {
                let source_col = exp.dataset.get_column_by_nature(*nature)?.ok_or_else(|| {
                    TGADomainError::InvalidOperation(format!(
                        "Experiment '{}' has no column with nature {:?}",
                        exp.meta.id, nature
                    ))
                })?;
                let unit = exp
                    .dataset
                    .schema
                    .columns
                    .get(&source_col)
                    .map(|column| column.unit)
                    .ok_or_else(|| {
                        TGADomainError::InvalidOperation(format!(
                            "Experiment '{}' has no unit metadata for column '{}'",
                            exp.meta.id, source_col
                        ))
                    })?;

                if let Some((_, expected)) = common_units.iter().find(|(role, _)| role == nature) {
                    if *expected != unit {
                        return Err(TGADomainError::IncompatibleUnits {
                            operation: "concat_into_vertical_stack".to_string(),
                            expected: format!("{} for {:?}", expected, nature),
                            actual: unit,
                        });
                    }
                } else {
                    common_units.push((*nature, unit));
                }

                select_exprs
                    .push(col(source_col).alias(UnitedDataset::canonical_column_name(*nature)));
            }

            let complete_case = all_horizontal(
                what_cols_take
                    .iter()
                    .map(|nature| col(UnitedDataset::canonical_column_name(*nature)).is_not_null())
                    .collect::<Vec<_>>(),
            );
            let selected = exp
                .dataset
                .frame
                .clone()
                .select(select_exprs)
                .filter(complete_case?)
                .with_column(lit(exp.meta.id.clone()).alias("exp_id"));

            frames.push(selected);
            meta.push(exp.meta.clone());
            units.insert(
                exp.meta.id.clone(),
                super::experiment_series2::KineticColumnUnits::from_dataset(&exp.dataset),
            );
        }

        let combined = concat(frames, UnionArgs::default())?;
        Ok(UnitedDataset::new(combined, meta).with_units(units))
    }

    /// Build a united dataset in "Polars Struct" layout.
    ///
    /// Source data:
    /// - starts from vertical-stack representation of selected experiments
    ///
    /// Output layout:
    /// - one row per `exp_id`
    /// - `data` column stores `List<Struct>` where each struct is one original point
    ///   with fields like `time`, `mass`, etc. (canonical names).
    ///
    /// This packs per-experiment trajectories into one cell, which is useful when
    /// you want to pass complete experiment curves as one grouped object.
    pub fn concat_into_polars_struct(
        &self,
        what_exp_to_take: Option<&[&str]>,
        what_cols_take: &[ColumnNature],
    ) -> Result<UnitedDataset, TGADomainError> {
        if what_cols_take.is_empty() {
            return Err(TGADomainError::InvalidOperation(
                "No column natures were requested".to_string(),
            ));
        }

        let stacked = self.concat_into_vertical_stack(what_exp_to_take, what_cols_take.to_vec())?;

        let struct_fields: Vec<Expr> = what_cols_take
            .iter()
            .map(|nature| col(UnitedDataset::canonical_column_name(*nature)))
            .collect();

        let aggregated = stacked
            .frame
            .clone()
            .group_by([col("exp_id")])
            .agg([as_struct(struct_fields).alias("data")]);

        Ok(UnitedDataset::new(aggregated, stacked.meta).with_units(stacked.units))
    }
    //==============================================================================================
    /// OUTPUT
    /// TAKING DATA FOR KINETIC METHODS
    pub fn create_kinetic_data_view(
        &self,
        what_exp_to_take: Option<&[&str]>,
        what_cols_take: Vec<ColumnNature>,
    ) -> Result<KineticDataView, TGADomainError> {
        let united = self.concat_into_vertical_stack(what_exp_to_take, what_cols_take.clone())?;
        info!("United dataset created");
        KineticDataView::from_united_dataset_by_nature(&united, what_cols_take)
    }

    pub fn create_kinetic_data_view_for_method(
        &self,
        what_exp_to_take: Option<&[&str]>,
        method: &IsoconversionalMethod,
    ) -> Result<KineticDataView, TGADomainError> {
        use crate::Kinetics::experimental_kinetics::experiment_series_main::ExperimentThermalMode;

        let experiments = self.select_experiments(what_exp_to_take)?;
        let isothermal_method = matches!(method, IsoconversionalMethod::FriedmanIntegral);
        for experiment in &experiments {
            let isothermal =
                experiment.meta.resolved_thermal_mode() == ExperimentThermalMode::Isothermal;
            if isothermal != isothermal_method {
                return Err(TGADomainError::InvalidOperation(format!(
                    "Method '{:?}' is {}, but experiment '{}' is declared {}",
                    method,
                    if isothermal_method {
                        "isothermal"
                    } else {
                        "non-isothermal"
                    },
                    experiment.meta.id,
                    if isothermal {
                        "isothermal"
                    } else {
                        "not isothermal"
                    },
                )));
            }
        }
        let what_cols_take = method.required_columns_by_nature();
        self.create_kinetic_data_view(what_exp_to_take, what_cols_take)
    }

    //==============================================================================================
    /// INPUT
    /// TAKING ISOCONVERSIONAL RESULT AND PUSHING IT INTO THE SERIES
    /// Create a new experiment from an isoconversional result and push it into the series.
    /// The experiment id is provided by the caller; all other meta fields are left empty.
    pub fn push_isoconversional_result(
        &mut self,
        result: &IsoconversionalResult,
        id: &str,
    ) -> Result<(), TGADomainError> {
        let dataset = result.to_tga_dataset()?;
        let exp = TGAExperiment::new(dataset).with_id(id);
        self.push(exp)?;
        info!("Isoconversional method results added!");
        Ok(())
    }

    /// Add a numeric column from Vec<f64> to a specific experiment by id.
    ///
    /// This is a thin series-level wrapper around `TGADataset::add_column_from_vec`.
    pub fn add_column_from_vec(
        &mut self,
        id: &str,
        name: &str,
        unit: Unit,
        nature: ColumnNature,
        data: Vec<f64>,
    ) -> Result<(), TGADomainError> {
        self.try_transform_by_id(id, move |exp| {
            let dataset = exp.dataset.add_column_from_vec(name, unit, nature, data)?;
            Ok(TGAExperiment {
                dataset,
                meta: exp.meta,
            })
        })?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Kinetics::experimental_kinetics::experiment_series_main::TGAExperiment;
    use crate::Kinetics::experimental_kinetics::testing_mod::VirtualTGA;

    fn build_series() -> TGASeries {
        let mut series = TGASeries::new();

        let v1 = VirtualTGA {
            time: vec![0.0, 1.0, 2.0],
            temperature: vec![300.0, 301.0, 302.0],
            mass: vec![10.0, 9.8, 9.6],
        };
        let v2 = VirtualTGA {
            time: vec![0.0, 1.0, 2.0],
            temperature: vec![310.0, 311.0, 312.0],
            mass: vec![8.0, 7.9, 7.8],
        };

        let d1 = crate::Kinetics::experimental_kinetics::one_experiment_dataset::TGADataset::create_from_synthetic_data(&v1).unwrap();
        let d2 = crate::Kinetics::experimental_kinetics::one_experiment_dataset::TGADataset::create_from_synthetic_data(&v2).unwrap();

        series
            .push(TGAExperiment::new(d1).with_id("exp_1"))
            .unwrap();
        series
            .push(TGAExperiment::new(d2).with_id("exp_2"))
            .unwrap();
        series
    }

    #[test]
    fn nonisothermal_method_is_rejected_for_explicit_isothermal_experiment() {
        let mut series = build_series();
        series
            .try_mutate_by_id("exp_1", |experiment| {
                experiment.meta.isothermal_temperature = Some(500.0);
                experiment.meta.thermal_mode =
                    crate::Kinetics::experimental_kinetics::experiment_series_main::ExperimentThermalMode::Isothermal;
                Ok(())
            })
            .unwrap();

        let error = series
            .create_kinetic_data_view_for_method(None, &IsoconversionalMethod::KAS)
            .unwrap_err();
        assert!(format!("{error}").contains("declared isothermal"));
    }

    #[test]
    fn friedman_integral_builds_view_for_isothermal_series_without_t_columns() {
        use crate::Kinetics::experimental_kinetics::kinetic_methods::Friedman::FriedmanIntegral;
        use crate::Kinetics::experimental_kinetics::kinetic_methods::KineticMethod;

        let mut series = build_series();
        let third = VirtualTGA {
            time: vec![0.0, 1.0, 2.0],
            temperature: vec![320.0, 321.0, 322.0],
            mass: vec![7.0, 6.9, 6.8],
        };
        let dataset = crate::Kinetics::experimental_kinetics::one_experiment_dataset::TGADataset::create_from_synthetic_data(&third).unwrap();
        series
            .push(TGAExperiment::new(dataset).with_id("exp_3"))
            .unwrap();

        for (index, id) in series.ids().into_iter().enumerate() {
            series
                .try_transform_by_id(&id, |experiment| {
                    let mut dataset = experiment
                        .dataset
                        .add_column_from_vec(
                            "eta",
                            Unit::Dimensionless,
                            ColumnNature::Conversion,
                            vec![0.1, 0.5, 0.9],
                        )?
                        .add_column_from_vec(
                            "deta_dt",
                            Unit::PerSecond,
                            ColumnNature::ConversionRate,
                            vec![0.1, 0.08, 0.04],
                        )?
                        .drop_column("temperature")?;
                    dataset.schema.eta = Some("eta".into());
                    dataset.schema.deta_dt = Some("deta_dt".into());
                    let mut meta = experiment.meta;
                    meta.isothermal_temperature = Some(500.0 + index as f64 * 20.0);
                    meta.thermal_mode = crate::Kinetics::experimental_kinetics::experiment_series_main::ExperimentThermalMode::Isothermal;
                    Ok(TGAExperiment {
                        dataset,
                        meta,
                    })
                })
                .unwrap();
        }

        let view = series
            .create_kinetic_data_view_for_method(None, &IsoconversionalMethod::FriedmanIntegral)
            .unwrap();
        assert_eq!(view.experiments.len(), 3);
        for experiment in &view.experiments {
            assert!(
                !series
                    .get_experiment_by_id(&experiment.meta.id)
                    .unwrap()
                    .dataset
                    .list_of_columns()
                    .contains(&"temperature".to_string())
            );
            assert_eq!(experiment.temperature.len(), experiment.time.len());
            assert!(
                experiment
                    .temperature
                    .iter()
                    .all(|temperature| *temperature
                        == experiment.meta.isothermal_temperature.unwrap())
            );
        }
        FriedmanIntegral.check_input(&view).unwrap();
    }

    #[test]
    fn concat_into_vertical_stack_builds_united_dataset() {
        let series = build_series();
        let united = series
            .concat_into_vertical_stack(None, vec![ColumnNature::Time, ColumnNature::Mass])
            .unwrap();

        let df = united.frame.collect().unwrap();
        assert_eq!(df.height(), 6);
        assert!(df.column("time").is_ok());
        assert!(df.column("mass").is_ok());
        assert!(df.column("exp_id").is_ok());
        assert_eq!(united.meta.len(), 2);
    }

    #[test]
    fn vertical_stack_allows_different_experiment_lengths() {
        let mut series = TGASeries::new();
        let long = VirtualTGA {
            time: vec![0.0, 1.0, 2.0],
            temperature: vec![300.0, 301.0, 302.0],
            mass: vec![10.0, 9.8, 9.6],
        };
        let short = VirtualTGA {
            time: vec![0.0, 1.0],
            temperature: vec![310.0, 311.0],
            mass: vec![8.0, 7.9],
        };
        for (id, data) in [("long", &long), ("short", &short)] {
            let dataset = crate::Kinetics::experimental_kinetics::one_experiment_dataset::TGADataset::create_from_synthetic_data(data).unwrap();
            series
                .push(TGAExperiment::new(dataset).with_id(id))
                .unwrap();
        }

        let united = series
            .concat_into_vertical_stack(None, vec![ColumnNature::Time, ColumnNature::Mass])
            .unwrap();
        assert_eq!(united.frame.collect().unwrap().height(), 5);
    }

    #[test]
    fn vertical_stack_rejects_mixed_units_for_requested_role() {
        let mut series = build_series();
        series.bind_time("exp_2", "time", Unit::Minute).unwrap();

        let error = match series
            .concat_into_vertical_stack(None, vec![ColumnNature::Time, ColumnNature::Mass])
        {
            Err(error) => error,
            Ok(_) => panic!("mixed time units must be rejected"),
        };

        assert!(matches!(
            error,
            TGADomainError::IncompatibleUnits {
                operation,
                expected,
                actual: Unit::Minute,
            } if operation == "concat_into_vertical_stack" && expected == "s for Time"
        ));
    }

    #[test]
    fn vertical_stack_uses_each_experiments_active_mass_binding() {
        let mut series = build_series();
        series
            .add_column_from_vec(
                "exp_1",
                "mass_smooth",
                Unit::Milligram,
                ColumnNature::Mass,
                vec![7.0, 6.0, 5.0],
            )
            .unwrap();
        series
            .bind_mass("exp_1", "mass_smooth", Unit::Milligram)
            .unwrap();

        let united = series
            .concat_into_vertical_stack(None, vec![ColumnNature::Time, ColumnNature::Mass])
            .unwrap();
        let mass = united
            .materialize_vertical_column_by_nature("exp_1", ColumnNature::Mass)
            .unwrap();

        assert_eq!(mass, vec![7.0, 6.0, 5.0]);
    }

    #[test]
    fn vertical_stack_applies_one_complete_case_mask_to_rate_edges() {
        let mut series = build_series();
        series
            .derive_rate("exp_1", "mass", "mass_rate", ColumnNature::MassRate)
            .unwrap();

        let united = series
            .concat_into_vertical_stack(
                Some(&["exp_1"]),
                vec![
                    ColumnNature::Time,
                    ColumnNature::Mass,
                    ColumnNature::MassRate,
                ],
            )
            .unwrap();
        let view = KineticDataView::from_united_dataset_by_nature(
            &united,
            vec![
                ColumnNature::Time,
                ColumnNature::Mass,
                ColumnNature::MassRate,
            ],
        )
        .unwrap();
        let experiment = &view.experiments[0];

        assert_eq!(experiment.time, vec![1.0]);
        assert_eq!(experiment.mass.as_deref(), Some(&[9.8][..]));
        assert_eq!(experiment.mass_rate.as_deref().unwrap().len(), 1);
        assert!(experiment.mass_rate.as_deref().unwrap()[0].is_finite());
    }

    #[test]
    fn concat_into_polars_struct_builds_struct_column() {
        let series = build_series();
        let united = series
            .concat_into_polars_struct(None, &[ColumnNature::Time, ColumnNature::Mass])
            .unwrap();

        let df = united.frame.collect().unwrap();
        assert_eq!(df.height(), 2);
        assert!(df.column("exp_id").is_ok());
        assert!(df.column("data").is_ok());

        let dtype = df.column("data").unwrap().dtype();
        assert!(
            matches!(dtype, DataType::List(inner) if matches!(inner.as_ref(), DataType::Struct(_)))
        );
    }

    #[test]
    fn materialize_vertical_column_by_nature_returns_values_for_one_id() {
        let series = build_series();
        let united = series
            .concat_into_vertical_stack(None, vec![ColumnNature::Time, ColumnNature::Mass])
            .unwrap();

        let mass = united
            .materialize_vertical_column_by_nature("exp_1", ColumnNature::Mass)
            .unwrap();
        assert_eq!(mass, vec![10.0, 9.8, 9.6]);
    }

    #[test]
    fn materialize_struct_field_by_nature_returns_values_for_one_id() {
        let series = build_series();
        let united = series
            .concat_into_polars_struct(None, &[ColumnNature::Time, ColumnNature::Mass])
            .unwrap();

        let time = united
            .materialize_struct_field_by_nature("exp_2", ColumnNature::Time)
            .unwrap();
        assert_eq!(time, vec![0.0, 1.0, 2.0]);
    }
}
