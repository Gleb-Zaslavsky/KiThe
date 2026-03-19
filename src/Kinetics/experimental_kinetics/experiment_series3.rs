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
            Some(ids) => ids
                .iter()
                .map(|id| {
                    let idx = self.exp_map.get(*id).ok_or_else(|| {
                        TGADomainError::InvalidOperation(format!(
                            "Experiment with id '{}' was not found",
                            id
                        ))
                    })?;
                    Ok(&self.experiments[*idx])
                })
                .collect(),
            None => Ok(self.experiments.iter().collect()),
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

        for exp in experiments {
            let mut select_exprs = Vec::with_capacity(what_cols_take.len());

            for nature in &what_cols_take {
                let source_col = exp.dataset.get_column_by_nature(*nature).ok_or_else(|| {
                    TGADomainError::InvalidOperation(format!(
                        "Experiment '{}' has no column with nature {:?}",
                        exp.meta.id, nature
                    ))
                })?;

                select_exprs
                    .push(col(source_col).alias(UnitedDataset::canonical_column_name(*nature)));
            }

            let selected = exp
                .dataset
                .frame
                .clone()
                .select(select_exprs)
                .with_column(lit(exp.meta.id.clone()).alias("exp_id"));

            frames.push(selected);
            meta.push(exp.meta.clone());
        }

        let combined = concat(frames, UnionArgs::default())?;
        Ok(UnitedDataset::new(combined, meta))
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

        Ok(UnitedDataset::new(aggregated, stacked.meta))
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
        self.push(exp);
        self.rebuild_index();
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

        series.push(TGAExperiment::new(d1).with_id("exp_1"));
        series.push(TGAExperiment::new(d2).with_id("exp_2"));
        series
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
