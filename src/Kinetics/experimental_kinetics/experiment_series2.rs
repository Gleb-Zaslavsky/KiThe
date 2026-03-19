use crate::Kinetics::experimental_kinetics::exp_engine_api::GoldenPipelineConfig;
use crate::Kinetics::experimental_kinetics::experiment_series_main::{
    ExperimentMeta, TGAExperiment, TGASeries,
};
use crate::Kinetics::experimental_kinetics::one_experiment_dataset::{
    ColumnNature, TGADataset, TGADomainError, TGASchema,
};
use crate::Kinetics::experimental_kinetics::testing_mod::{AdvancedTGAConfig, VirtualTGA};
use polars::prelude::*;
use std::collections::HashMap;
//=================================================================================
/// INTERFACE FOR KINETIC METHODS
#[derive(Clone)]
pub struct UnitedDataset {
    pub frame: LazyFrame,
    pub meta: Vec<ExperimentMeta>,
}

impl UnitedDataset {
    /// Canonical output column names used by united datasets created from
    /// `ColumnNature`-based selection in series-level APIs.
    pub fn canonical_column_name(nature: ColumnNature) -> &'static str {
        match nature {
            ColumnNature::Time => "time",
            ColumnNature::Temperature => "temperature",
            ColumnNature::Mass => "mass",
            ColumnNature::Conversion => "conversion",
            ColumnNature::DimensionlessMass => "dimensionless_mass",
            ColumnNature::MassRate => "mass_rate",
            ColumnNature::ConversionRate => "conversion_rate",
            ColumnNature::TemperatureRate => "temperature_rate",
            ColumnNature::DimensionlessMassRate => "dimensionless_mass_rate",
            ColumnNature::Unknown => "unknown",
            ColumnNature::ActivationEnergy => "activation_energy",
            ColumnNature::PredexFactor => "predex_factor",
            ColumnNature::R2 => "r2",
        }
    }

    pub fn empty() -> Self {
        Self {
            frame: DataFrame::default().lazy(),
            meta: Vec::new(),
        }
    }

    pub fn new(frame: LazyFrame, meta: Vec<ExperimentMeta>) -> Self {
        Self { frame, meta }
    }

    pub fn frame(&self) -> &LazyFrame {
        &self.frame
    }

    pub fn frame_mut(&mut self) -> &mut LazyFrame {
        &mut self.frame
    }

    pub fn set_frame(&mut self, frame: LazyFrame) {
        self.frame = frame;
    }

    pub fn meta(&self) -> &[ExperimentMeta] {
        &self.meta
    }

    pub fn set_meta(&mut self, meta: Vec<ExperimentMeta>) {
        self.meta = meta;
    }

    pub fn push_meta(&mut self, meta: ExperimentMeta) {
        self.meta.push(meta);
    }

    pub fn write_united_dataset_to_csv(&self, path: &str) -> PolarsResult<()> {
        let df = &mut self.frame.clone().collect()?;
        let mut file = std::fs::File::create(path)?;
        CsvWriter::new(&mut file).include_header(true).finish(df)?;
        Ok(())
    }
    pub fn write_united_dataset_to_parquet(&self, path: &str) -> PolarsResult<()> {
        let df = &mut self.frame.clone().collect()?;
        let mut file = std::fs::File::create(path)?;
        ParquetWriter::new(&mut file).finish(df)?;
        Ok(())
    }

    /// Materialize one numeric column from a vertical-stack united dataset
    /// for one experiment (`exp_id`).
    pub fn materialize_vertical_column_by_name(
        &self,
        id: &str,
        column_name: &str,
    ) -> Result<Vec<f64>, TGADomainError> {
        let df = self
            .frame
            .clone()
            .filter(col("exp_id").eq(lit(id)))
            .select([col(column_name)])
            .collect()?;

        if df.height() == 0 {
            return Err(TGADomainError::InvalidOperation(format!(
                "No rows for experiment id '{}'",
                id
            )));
        }

        let out = df
            .column(column_name)?
            .f64()?
            .into_no_null_iter()
            .collect::<Vec<f64>>();
        Ok(out)
    }

    /// Materialize one numeric column from a vertical-stack united dataset
    /// by semantic `ColumnNature`.
    pub fn materialize_vertical_column_by_nature(
        &self,
        id: &str,
        nature: ColumnNature,
    ) -> Result<Vec<f64>, TGADomainError> {
        self.materialize_vertical_column_by_name(id, Self::canonical_column_name(nature))
    }

    /// Materialize one numeric field from a struct-based united dataset
    /// for one experiment (`exp_id`), where points are stored in `data: List<Struct>`.
    pub fn materialize_struct_field_by_name(
        &self,
        id: &str,
        field_name: &str,
    ) -> Result<Vec<f64>, TGADomainError> {
        let df = self
            .frame
            .clone()
            .filter(col("exp_id").eq(lit(id)))
            .select([col("data")])
            .collect()?;

        if df.height() == 0 {
            return Err(TGADomainError::InvalidOperation(format!(
                "No rows for experiment id '{}'",
                id
            )));
        }

        let list_col = df.column("data")?.list()?;
        let struct_series = list_col.get_as_series(0).ok_or_else(|| {
            TGADomainError::InvalidOperation(
                "Struct dataset contains no list payload for selected experiment".to_string(),
            )
        })?;

        let field = struct_series.struct_()?.field_by_name(field_name)?;
        let out = field.f64()?.into_no_null_iter().collect::<Vec<f64>>();
        Ok(out)
    }

    /// Materialize one numeric field from a struct-based united dataset
    /// by semantic `ColumnNature`.
    pub fn materialize_struct_field_by_nature(
        &self,
        id: &str,
        nature: ColumnNature,
    ) -> Result<Vec<f64>, TGADomainError> {
        self.materialize_struct_field_by_name(id, Self::canonical_column_name(nature))
    }
}
pub fn write_lazy_frame_eagerly(lf: LazyFrame, path: &str) -> PolarsResult<()> {
    let df = &mut lf.clone().collect()?;
    let mut file = std::fs::File::create(path)?;
    CsvWriter::new(&mut file).include_header(true).finish(df)?;
    Ok(())
}
//====================================================================================
//=========================================================
pub struct SampledColumns {
    cols: Vec<Vec<f64>>,
    names: Vec<String>,
}
impl SampledColumns {
    pub fn new(cols: Vec<Vec<f64>>, names: Vec<String>) -> Self {
        Self { cols, names }
    }

    pub fn new_empty() -> Self {
        Self {
            cols: Vec::new(),
            names: Vec::new(),
        }
    }

    pub fn push_column(&mut self, col: Vec<f64>, name: String) {
        self.cols.push(col);
        self.names.push(name);
    }

    pub fn column(&self, name: &str) -> Option<Vec<f64>> {
        self.names
            .iter()
            .position(|n| n == name)
            .map(|idx| self.cols[idx].clone())
    }

    pub fn column_names(&self) -> Vec<String> {
        self.names.clone()
    }

    pub fn len(&self) -> usize {
        self.cols.first().map(|c| c.len()).unwrap_or(0)
    }
}

impl TGAExperiment {
    /// Thin wrapper for dataset monotonicity diagnostics on bound time column.
    pub fn monotony_of_time_check(&self) -> Result<Vec<f64>, TGADomainError> {
        self.dataset.monotony_of_time_check()
    }
}

impl TGASeries {
    /// Thin wrapper for per-experiment time monotonicity diagnostics.
    pub fn monotony_of_time_check(&self, id: &str) -> Result<Vec<f64>, TGADomainError> {
        self.try_apply_by_id(id, |exp| exp.monotony_of_time_check())
    }

    pub fn unite_datasets(&self) -> Result<UnitedDataset, TGADomainError> {
        if self.experiments.is_empty() {
            return Ok(UnitedDataset::empty());
        }

        let mut united_df = DataFrame::default();
        let mut expected_height: Option<usize> = None;
        let mut metas: Vec<ExperimentMeta> = Vec::with_capacity(self.experiments.len());

        for (idx, exp) in self.experiments.iter().enumerate() {
            let df = exp.dataset.frame.clone().collect()?;

            if let Some(height) = expected_height {
                if df.height() != height {
                    return Err(TGADomainError::InvalidOperation(format!(
                        "Cannot unite datasets with different row counts: '{}' has {}, expected {}",
                        exp.meta.id,
                        df.height(),
                        height
                    )));
                }
            } else {
                expected_height = Some(df.height());
            }

            let prefix = if exp.meta.id.is_empty() {
                format!("exp{}", idx)
            } else {
                exp.meta.id.clone()
            };

            let mut renamed_columns = Vec::with_capacity(df.width());
            for old_name in df.get_column_names() {
                let mut column = df.column(old_name.as_str())?.clone();
                let new_name = format!("{}_{}", prefix, old_name.as_str());
                column.rename(new_name.into());
                renamed_columns.push(column);
            }

            if united_df.width() == 0 {
                let height = renamed_columns.iter().map(|c| c.len()).max().unwrap_or(0);
                united_df = DataFrame::new(height, renamed_columns)?;
            } else {
                united_df.hstack_mut(&renamed_columns)?;
            }

            metas.push(exp.meta.clone());
        }

        Ok(UnitedDataset::new(united_df.lazy(), metas))
    }

    /// Create a new experiment from an existing one by selecting specific columns.
    ///
    /// # What happens:
    /// 1. Selects specified columns from parent experiment's LazyFrame using Polars `select()`
    /// 2. Filters schema to include only metadata for selected columns
    /// 3. Updates semantic fields (time, temperature, mass, etc.) - keeps only if column is selected
    /// 4. Clones parent metadata (heating_rate, isothermal_temperature, comment) but assigns new ID
    /// 5. Creates fresh history (no operations carried over)
    /// 6. Pushes new experiment to the series
    ///
    /// # Use case:
    /// When you need a subset of columns from an experiment as a separate entity,
    /// e.g., extracting only time and smoothed mass for further processing.
    ///
    /// # Arguments:
    /// * `parent_idx` - Index of parent experiment in the series
    /// * `new_id` - Identifier for the new experiment
    /// * `columns` - Slice of column names to extract
    ///
    /// # Errors:
    /// Returns error if parent index is out of bounds, column list is empty, or column not found.
    pub fn create_experiment_from_columns(
        &mut self,
        parent_idx: usize,
        new_id: String,
        columns: &[&str],
    ) -> Result<(), TGADomainError> {
        let parent = self.experiments.get(parent_idx).ok_or_else(|| {
            TGADomainError::InvalidOperation(format!(
                "Experiment index {} out of bounds",
                parent_idx
            ))
        })?;

        let new_exp = Self::extract_columns_to_new_experiment(parent, new_id, columns)?;
        self.push(new_exp);
        Ok(())
    }

    /// Core extraction logic: creates new TGAExperiment from parent by selecting columns.
    ///
    /// # Implementation details:
    /// - Uses Polars `select()` to create new LazyFrame with only requested columns
    /// - Clones ColumnMeta for each selected column from parent schema
    /// - Semantic fields are filtered: preserved only if their column is in selection
    /// - New experiment gets empty operation history
    /// - Parent metadata cloned except for ID which is replaced
    fn extract_columns_to_new_experiment(
        parent: &TGAExperiment,
        new_id: String,
        columns: &[&str],
    ) -> Result<TGAExperiment, TGADomainError> {
        if columns.is_empty() {
            return Err(TGADomainError::InvalidOperation(
                "Column list cannot be empty".to_string(),
            ));
        }

        let select_exprs: Vec<Expr> = columns.iter().map(|&c| col(c)).collect();
        let new_frame = parent.dataset.frame.clone().select(select_exprs);

        let mut new_schema_columns = HashMap::new();
        for &col_name in columns {
            if let Some(meta) = parent.dataset.schema.columns.get(col_name) {
                new_schema_columns.insert(col_name.to_string(), meta.clone());
            } else {
                return Err(TGADomainError::ColumnNotFound(col_name.to_string()));
            }
        }

        let new_schema = TGASchema {
            columns: new_schema_columns,
            time: Self::filter_semantic_field(&parent.dataset.schema.time, columns),
            temperature: Self::filter_semantic_field(&parent.dataset.schema.temperature, columns),
            mass: Self::filter_semantic_field(&parent.dataset.schema.mass, columns),
            alpha: Self::filter_semantic_field(&parent.dataset.schema.alpha, columns),
            dm_dt: Self::filter_semantic_field(&parent.dataset.schema.dm_dt, columns),
            eta: Self::filter_semantic_field(&parent.dataset.schema.eta, columns),
            deta_dt: Self::filter_semantic_field(&parent.dataset.schema.deta_dt, columns),
            dalpha_dt: Self::filter_semantic_field(&parent.dataset.schema.dalpha_dt, columns),
            dT_dt: Self::filter_semantic_field(&parent.dataset.schema.dT_dt, columns),
        };

        let new_dataset = TGADataset {
            frame: new_frame,
            schema: new_schema,
            oneframeplot: None,
            history_of_operations:
                crate::Kinetics::experimental_kinetics::one_experiment_dataset::History::new(),
        };

        let mut new_meta = parent.meta.clone();
        new_meta.id = new_id;

        Ok(TGAExperiment {
            dataset: new_dataset,
            meta: new_meta,
        })
    }

    /// Helper to filter semantic field references.
    /// Returns Some(name) only if the field exists and its column is in the selection.
    fn filter_semantic_field(field: &Option<String>, columns: &[&str]) -> Option<String> {
        field.as_ref().and_then(|name| {
            if columns.contains(&name.as_str()) {
                Some(name.clone())
            } else {
                None
            }
        })
    }

    //=======================================================================================
    //  GOLDEN PIPELINE
    //================================================================================
    pub fn drop_nulls(&mut self, id: &str) -> Result<(), TGADomainError> {
        self.try_mutate_by_id(id, |exp| exp.drop_nulls())
    }

    pub fn apply_golden_pipeline_inner(
        &mut self,
        id: &str,
        config: GoldenPipelineConfig,
    ) -> Result<Vec<String>, TGADomainError> {
        let mut new_columns = None;
        self.try_transform_by_id(id, |exp| {
            let (new_exp, cols) = exp.apply_golden_pipeline(config)?;
            new_columns = Some(cols);
            Ok(new_exp)
        })?;
        Ok(new_columns.unwrap_or_default())
    }
    pub fn apply_golden_pipeline(
        &mut self,
        id: &str,
        config: GoldenPipelineConfig,
    ) -> Result<(), TGADomainError> {
        let do_we_need_new_exp = config.save_to_new_experiment;
        let do_we_need_drop_old_exp = config.del_old_experiment;
        let new_columns = self.apply_golden_pipeline_inner(id, config)?;

        if do_we_need_new_exp {
            let idx = self.index_by_id(id)?;
            let new_id = format!("proceeded_{}", id);
            let new_columns: Vec<&str> = new_columns.iter().map(|x| x.as_str()).collect();
            self.create_experiment_from_columns(idx, new_id.clone(), &new_columns)?;
            // self.trim_null_edges(&new_id)?;
            self.drop_nulls(&new_id)?;
            self.assert_no_nulls(&new_id, "final")?;
            if do_we_need_drop_old_exp {
                self.drop_experiment(idx);
            }
            self.rebuild_index();
        }
        Ok(())
    }

    pub fn column_samples_for_all_experiment_for_plotting(
        &mut self,
        n_points: usize,
    ) -> Result<(HashMap<String, Vec<f64>>, SampledColumns), TGADomainError> {
        let mut data_map: HashMap<String, Vec<f64>> = HashMap::new();
        let mut sampled_cols = SampledColumns::new_empty();

        for id in self.ids() {
            let idx = self.index_by_id(&id)?;
            let frame = self.experiments[idx].dataset.frame.clone().collect()?;
            for col in self.list_of_columns(&id)? {
                let short_name = format!("{}_{}", id, col);
                let eager_col: Vec<f64> = frame.column(&col)?.f64()?.into_no_null_iter().collect();
                let sampled_col = if n_points < eager_col.len() {
                    self.sample_column(&id, &col, None, n_points)?
                } else {
                    eager_col
                };
                sampled_cols.push_column(sampled_col.clone(), short_name.clone());
                data_map.insert(short_name, sampled_col);
            }
        }
        Ok((data_map, sampled_cols))
    }
}

/// Standalone function to create a new TGAExperiment from an existing one by selecting columns.
///
/// # What happens:
/// 1. Validates column list is non-empty
/// 2. Creates new LazyFrame with selected columns using Polars `select()`
/// 3. Builds new schema containing only metadata for selected columns
/// 4. Filters semantic bindings (time, mass, etc.) - keeps only if column selected
/// 5. Clones parent ExperimentMeta but replaces ID with new_id
/// 6. Returns new TGAExperiment with fresh history
///
/// # Difference from TGASeries method:
/// This version doesn't require mutable access to TGASeries. Use when you want
/// to create a derived experiment without immediately adding it to a series.
///
/// # Arguments:
/// * `parent` - Reference to parent experiment
/// * `new_id` - Identifier for new experiment
/// * `columns` - Column names to extract
///
/// # Returns:
/// New TGAExperiment containing only selected columns with properly filtered schema.
///
/// # Errors:
/// Returns error if column list is empty or any column not found in parent.
///
/// # Example:
/// ```ignore
/// let parent_exp = series.experiments[0];
/// let new_exp = create_experiment_from_columns(
///     &parent_exp,
///     "smoothed_data".to_string(),
///     &["time", "mass_smooth", "temperature"]
/// )?;
/// ```
pub fn create_experiment_from_columns(
    parent: &TGAExperiment,
    new_id: String,
    columns: &[&str],
) -> Result<TGAExperiment, TGADomainError> {
    if columns.is_empty() {
        return Err(TGADomainError::InvalidOperation(
            "Column list cannot be empty".to_string(),
        ));
    }

    let select_exprs: Vec<Expr> = columns.iter().map(|&c| col(c)).collect();
    let new_frame = parent.dataset.frame.clone().select(select_exprs);

    let mut new_schema_columns = HashMap::new();
    for &col_name in columns {
        if let Some(meta) = parent.dataset.schema.columns.get(col_name) {
            new_schema_columns.insert(col_name.to_string(), meta.clone());
        } else {
            return Err(TGADomainError::ColumnNotFound(col_name.to_string()));
        }
    }

    let filter_field = |field: &Option<String>| -> Option<String> {
        field.as_ref().and_then(|name| {
            if columns.contains(&name.as_str()) {
                Some(name.clone())
            } else {
                None
            }
        })
    };

    let new_schema = TGASchema {
        columns: new_schema_columns,
        time: filter_field(&parent.dataset.schema.time),
        temperature: filter_field(&parent.dataset.schema.temperature),
        mass: filter_field(&parent.dataset.schema.mass),
        alpha: filter_field(&parent.dataset.schema.alpha),
        dm_dt: filter_field(&parent.dataset.schema.dm_dt),
        eta: filter_field(&parent.dataset.schema.eta),
        deta_dt: filter_field(&parent.dataset.schema.deta_dt),
        dalpha_dt: filter_field(&parent.dataset.schema.dalpha_dt),
        dT_dt: filter_field(&parent.dataset.schema.dT_dt),
    };

    let new_dataset = TGADataset {
        frame: new_frame,
        schema: new_schema,
        oneframeplot: None,
        history_of_operations:
            crate::Kinetics::experimental_kinetics::one_experiment_dataset::History::new(),
    };

    let mut new_meta = parent.meta.clone();
    new_meta.id = new_id;

    Ok(TGAExperiment {
        dataset: new_dataset,
        meta: new_meta,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Kinetics::experimental_kinetics::one_experiment_dataset::{
        ColumnMeta, ColumnOrigin, History, Unit,
    };
    use polars::prelude::*;

    fn make_mock_dataset() -> TGADataset {
        let df = df! {
            "time" => &[0.0, 1.0, 2.0, 3.0, 4.0],
            "temperature" => &[300.0, 310.0, 320.0, 330.0, 340.0],
            "mass" => &[10.0, 9.8, 9.5, 9.2, 9.0],
            "mass_smooth" => &[10.0, 9.85, 9.55, 9.25, 9.0],
        }
        .unwrap();

        let mut columns = HashMap::new();
        columns.insert(
            "time".to_string(),
            ColumnMeta {
                name: "time".to_string(),
                unit: Unit::Second,
                origin: ColumnOrigin::Raw,
                nature: ColumnNature::Time,
            },
        );
        columns.insert(
            "temperature".to_string(),
            ColumnMeta {
                name: "temperature".to_string(),
                unit: Unit::Kelvin,
                origin: ColumnOrigin::Raw,
                nature: ColumnNature::Temperature,
            },
        );
        columns.insert(
            "mass".to_string(),
            ColumnMeta {
                name: "mass".to_string(),
                unit: Unit::Milligram,
                origin: ColumnOrigin::Raw,
                nature: ColumnNature::Mass,
            },
        );
        columns.insert(
            "mass_smooth".to_string(),
            ColumnMeta {
                name: "mass_smooth".to_string(),
                unit: Unit::Milligram,
                origin: ColumnOrigin::PolarsDerived,
                nature: ColumnNature::Mass,
            },
        );

        let schema = TGASchema {
            columns,
            time: Some("time".to_string()),
            temperature: Some("temperature".to_string()),
            mass: Some("mass".to_string()),
            alpha: None,
            dm_dt: None,
            eta: None,
            deta_dt: None,
            dalpha_dt: None,
            dT_dt: None,
        };

        TGADataset {
            frame: df.lazy(),
            schema,
            oneframeplot: None,
            history_of_operations: History::new(),
        }
    }

    #[test]
    fn test_create_experiment_from_columns_basic() {
        let dataset = make_mock_dataset();
        let mut meta = ExperimentMeta::new();
        meta.id = "parent".to_string();
        meta.heating_rate = Some(10.0);
        let parent = TGAExperiment { dataset, meta };

        let new_exp =
            create_experiment_from_columns(&parent, "child".to_string(), &["time", "mass_smooth"])
                .unwrap();

        assert_eq!(new_exp.meta.id, "child");
        assert_eq!(new_exp.meta.heating_rate, Some(10.0));

        let df = new_exp.dataset.frame.collect().unwrap();
        assert_eq!(df.width(), 2);
        assert!(df.column("time").is_ok());
        assert!(df.column("mass_smooth").is_ok());
        assert!(df.column("temperature").is_err());
    }

    #[test]
    fn test_schema_filtered_correctly() {
        let dataset = make_mock_dataset();
        let parent = TGAExperiment {
            dataset,
            meta: ExperimentMeta::new(),
        };

        let new_exp =
            create_experiment_from_columns(&parent, "test".to_string(), &["time", "temperature"])
                .unwrap();

        assert_eq!(new_exp.dataset.schema.columns.len(), 2);
        assert!(new_exp.dataset.schema.columns.contains_key("time"));
        assert!(new_exp.dataset.schema.columns.contains_key("temperature"));
        assert!(!new_exp.dataset.schema.columns.contains_key("mass"));
    }

    #[test]
    fn test_semantic_fields_filtered() {
        let dataset = make_mock_dataset();
        let parent = TGAExperiment {
            dataset,
            meta: ExperimentMeta::new(),
        };

        let new_exp =
            create_experiment_from_columns(&parent, "test".to_string(), &["time", "mass_smooth"])
                .unwrap();

        assert_eq!(new_exp.dataset.schema.time, Some("time".to_string()));
        assert_eq!(new_exp.dataset.schema.temperature, None);
        assert_eq!(new_exp.dataset.schema.mass, None);
    }

    #[test]
    fn test_empty_columns_error() {
        let dataset = make_mock_dataset();
        let parent = TGAExperiment {
            dataset,
            meta: ExperimentMeta::new(),
        };

        let result = create_experiment_from_columns(&parent, "test".to_string(), &[]);

        assert!(result.is_err());
    }

    #[test]
    fn test_nonexistent_column_error() {
        let dataset = make_mock_dataset();
        let parent = TGAExperiment {
            dataset,
            meta: ExperimentMeta::new(),
        };

        let result =
            create_experiment_from_columns(&parent, "test".to_string(), &["time", "nonexistent"]);

        assert!(result.is_err());
    }

    #[test]
    fn test_series_method_adds_to_experiments() {
        let dataset = make_mock_dataset();
        let exp = TGAExperiment {
            dataset,
            meta: ExperimentMeta {
                id: "exp1".to_string(),
                heating_rate: Some(5.0),
                isothermal_temperature: None,
                comment: Some("test".to_string()),
            },
        };

        let mut series = TGASeries {
            experiments: vec![exp],
            exp_map: HashMap::new(),
        };

        series
            .create_experiment_from_columns(0, "exp2".to_string(), &["time", "temperature"])
            .unwrap();

        assert_eq!(series.experiments.len(), 2);
        assert_eq!(series.experiments[1].meta.id, "exp2");
        assert_eq!(series.experiments[1].meta.heating_rate, Some(5.0));
        assert_eq!(series.experiments[1].meta.comment, Some("test".to_string()));
    }

    #[test]
    fn test_series_method_invalid_index() {
        let mut series = TGASeries {
            experiments: vec![],
            exp_map: HashMap::new(),
        };

        let result = series.create_experiment_from_columns(0, "test".to_string(), &["time"]);

        assert!(result.is_err());
    }

    #[test]
    fn test_column_metadata_preserved() {
        let dataset = make_mock_dataset();
        let parent = TGAExperiment {
            dataset,
            meta: ExperimentMeta::new(),
        };

        let new_exp =
            create_experiment_from_columns(&parent, "test".to_string(), &["mass_smooth"]).unwrap();

        let meta = new_exp.dataset.schema.columns.get("mass_smooth").unwrap();
        assert_eq!(meta.unit, Unit::Milligram);
        assert_eq!(meta.origin, ColumnOrigin::PolarsDerived);
    }
}
