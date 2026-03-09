use crate::Kinetics::experimental_kinetics::LSQSplines::SolverKind;
use crate::Kinetics::experimental_kinetics::exp_engine_api::{
    GoldenPipelineConfig, PlotSeries, Ranges, ViewRange, XY,
};
use crate::Kinetics::experimental_kinetics::lowess_wrapper::LowessConfig;
use crate::Kinetics::experimental_kinetics::one_experiment_dataset::ColumnRole;
use crate::Kinetics::experimental_kinetics::one_experiment_dataset::History;
use crate::Kinetics::experimental_kinetics::one_experiment_dataset::{
    ColumnHistory, ColumnMeta, ColumnNature, ColumnOrigin, OperationRecord, TGADataset,
    TGADomainError, TGASchema, UnaryOp, Unit,
};
use crate::Kinetics::experimental_kinetics::splines::SplineKind;
use crate::Kinetics::experimental_kinetics::testing_mod::VirtualTGA;
use polars::prelude::{
    Column, CsvWriter, DataFrame, DataType, IntoLazy, LazyCsvReader, LazyFileListReader, NamedFrom,
    PlRefPath, SerWriter, Series,
};
use std::collections::HashMap;
use std::io::Write;
use std::path::Path;

#[derive(Clone, Debug, Default)]
/// ExperimentMeta data structure for the experimental kinetics API layer.
/// It groups domain data and exposes a stable facade over
/// dataset-level processing modules used throughout KiThe.
pub struct ExperimentMeta {
    /// experiment identifier (e.g. "β = 10 K/min")
    pub id: String,

    /// Heating rate for non-isothermal experiments (K/min)
    pub heating_rate: Option<f64>,

    /// Isothermal temperature (K)
    pub isothermal_temperature: Option<f64>,

    /// Free-form user comment
    pub comment: Option<String>,
}
impl ExperimentMeta {
    /// new API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn new() -> Self {
        Self {
            id: String::new(),
            heating_rate: None,
            isothermal_temperature: None,
            comment: None,
        }
    }
}

#[derive(Clone, Debug)]
/// TGAExperiment data structure for the experimental kinetics API layer.
/// It groups domain data and exposes a stable facade over
/// dataset-level processing modules used throughout KiThe.
pub struct TGAExperiment {
    pub dataset: TGADataset,
    pub meta: ExperimentMeta,
}

impl TGAExperiment {
    /// new API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn new(dataset: TGADataset) -> Self {
        Self {
            dataset,
            meta: ExperimentMeta::default(),
        }
    }

    /// with_id API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn with_id<S: Into<String>>(mut self, id: S) -> Self {
        self.meta.id = id.into();
        self
    }

    /// with_heating_rate API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn with_heating_rate(mut self, beta: f64) -> Self {
        self.meta.heating_rate = Some(beta);
        self
    }

    /// with_isothermal_temperature API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn with_isothermal_temperature(mut self, t: f64) -> Self {
        self.meta.isothermal_temperature = Some(t);
        self
    }

    /// with_comment API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn with_comment<S: Into<String>>(mut self, c: S) -> Self {
        self.meta.comment = Some(c.into());
        self
    }

    //==============================================================================
    // THIN WRAPPERS AROUND TGADataset FUNCTIONS
    //==============================================================================
    // INPUT/OUTPUT
    /// from_csv API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn from_csv(
        path: &str,
        time_col: &str,
        temp_col: &str,
        mass_col: &str,
    ) -> Result<Self, TGADomainError> {
        let dataset = TGADataset::from_csv(path, time_col, temp_col, mass_col)?;
        Ok(Self {
            dataset,
            meta: ExperimentMeta::new(),
        })
    }

    /// from_csv_raw API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn from_csv_raw(path: &Path) -> Result<Self, TGADomainError> {
        let dataset = TGADataset::from_csv_raw(path)?;
        Ok(Self {
            dataset,
            meta: ExperimentMeta::new(),
        })
    }

    /// from_csv_with_units API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn from_csv_with_units(path: &Path) -> Result<Self, TGADomainError> {
        let dataset = TGADataset::from_csv_with_units(path)?;
        Ok(Self {
            dataset,
            meta: ExperimentMeta::new(),
        })
    }

    /// from_csv_universal API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn from_csv_universal(path: &Path) -> Result<Self, TGADomainError> {
        let dataset = TGADataset::from_csv_universal(path)?;
        Ok(Self {
            dataset,
            meta: ExperimentMeta::new(),
        })
    }
    /// to_csv_with_units API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn to_csv_with_units(&self, path: &Path) -> Result<(), TGADomainError> {
        self.dataset.to_csv_with_units(path)?;
        Ok(())
    }

    /// to_csv_with_units_selected API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn to_csv_with_units_selected(
        &self,
        path: &Path,
        columns: &[&str],
    ) -> Result<(), TGADomainError> {
        self.dataset.to_csv_with_units_selected(path, columns)?;
        Ok(())
    }
    /// bind_time API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn bind_time(self, col: &str, unit: Unit) -> Result<Self, TGADomainError> {
        let dataset = self.dataset.bind_time(col, unit)?;
        Ok(Self {
            dataset,
            meta: self.meta,
        })
    }

    /// bind_temperature API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn bind_temperature(self, col: &str, unit: Unit) -> Result<Self, TGADomainError> {
        let dataset = self.dataset.bind_temperature(col, unit)?;
        Ok(Self {
            dataset,
            meta: self.meta,
        })
    }

    /// bind_mass API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn bind_mass(self, col: &str, unit: Unit) -> Result<Self, TGADomainError> {
        let dataset = self.dataset.bind_mass(col, unit)?;
        Ok(Self {
            dataset,
            meta: self.meta,
        })
    }

    /// bind_column API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn bind_column(
        self,
        role: ColumnRole,
        name: &str,
        unit: Unit,
    ) -> Result<Self, TGADomainError> {
        let dataset = self.dataset.bind_column(role, name, unit)?;
        Ok(Self {
            dataset,
            meta: self.meta,
        })
    }

    //================================================================
    // BASIC TRANSFORMATIONS: COLUMN CUT AND FILTER
    /// with_column_expr API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn with_column_expr(self, meta: ColumnMeta, expr: polars::prelude::Expr) -> Self {
        let dataset = self.dataset.with_column_expr(meta, expr);
        Self {
            dataset,
            meta: self.meta,
        }
    }

    /// filter_rows API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn filter_rows(self, predicate: polars::prelude::Expr) -> Self {
        let dataset = self.dataset.filter_rows(predicate);
        Self {
            dataset,
            meta: self.meta,
        }
    }

    /// filter_by_mask API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn filter_by_mask(self, mask: polars::prelude::Expr) -> Self {
        let dataset = self.dataset.filter_by_mask(mask);
        Self {
            dataset,
            meta: self.meta,
        }
    }

    /// cut_interval API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn cut_interval(self, colmn: &str, from: f64, to: f64) -> Self {
        let dataset = self.dataset.cut_interval(colmn, from, to);
        Self {
            dataset,
            meta: self.meta,
        }
    }

    /// cut_time_interval API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn cut_time_interval(self, from: f64, to: f64) -> Self {
        let dataset = self.dataset.cut_time_interval(from, to);
        Self {
            dataset,
            meta: self.meta,
        }
    }

    /// cut_temperature_interval API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn cut_temperature_interval(self, from: f64, to: f64) -> Self {
        let dataset = self.dataset.cut_temperature_interval(from, to);
        Self {
            dataset,
            meta: self.meta,
        }
    }

    /// cut_mass_interval API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn cut_mass_interval(self, from: f64, to: f64) -> Self {
        let dataset = self.dataset.cut_mass_interval(from, to);
        Self {
            dataset,
            meta: self.meta,
        }
    }

    /// cut_before_time API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn cut_before_time(self, t_start: f64) -> Self {
        let dataset = self.dataset.cut_before_time(t_start);
        Self {
            dataset,
            meta: self.meta,
        }
    }

    /// trim_range API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn trim_range(self, column: &str, from: f64, to: f64) -> Self {
        let dataset = self.dataset.trim_range(column, from, to);
        Self {
            dataset,
            meta: self.meta,
        }
    }

    /// trim_edges API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn trim_edges(self, left: usize, right: usize) -> Self {
        let dataset = self.dataset.trim_edges(left, right);
        Self {
            dataset,
            meta: self.meta,
        }
    }

    /// trim_null_edges API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn trim_null_edges(self) -> Self {
        let dataset = self.dataset.trim_null_edges();
        Self {
            dataset,
            meta: self.meta,
        }
    }

    /// trim_null_edges_for_columns API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn trim_null_edges_for_columns(self, cols_to_trim: &[&str]) -> Self {
        let dataset = self.dataset.trim_null_edges_for_columns(cols_to_trim);
        Self {
            dataset,
            meta: self.meta,
        }
    }

    /// rename_column API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn rename_column(self, old: &str, new: &str) -> Result<Self, TGADomainError> {
        let dataset = self.dataset.rename_column(old, new)?;
        Ok(Self {
            dataset,
            meta: self.meta,
        })
    }

    /// drop_column API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn drop_column(self, col_name: &str) -> Result<Self, TGADomainError> {
        let dataset = self.dataset.drop_column(col_name)?;
        Ok(Self {
            dataset,
            meta: self.meta,
        })
    }

    //======================================================================
    // UNIT TRANSFORMATIONS
    /// celsius_to_kelvin API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn celsius_to_kelvin(self) -> Self {
        let dataset = self.dataset.celsius_to_kelvin();
        Self {
            dataset,
            meta: self.meta,
        }
    }

    /// seconds_to_hours API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn seconds_to_hours(self) -> Self {
        let dataset = self.dataset.seconds_to_hours();
        Self {
            dataset,
            meta: self.meta,
        }
    }

    //=======================================================================
    // ARBITRARY ALGEBRAIC TRANSFORMATIONS ON COLUMN DATA
    /// scale_columns API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn scale_columns(self, cols: &[&str], factor: f64) -> Self {
        let dataset = self.dataset.scale_columns(cols, factor);
        Self {
            dataset,
            meta: self.meta,
        }
    }
    pub fn scale_column(self, col: &str, factor: f64) -> Self {
        let dataset = self.dataset.scale_column(col, factor);
        Self {
            dataset,
            meta: self.meta,
        }
    }
    pub fn scale_column_in_its_range(self, colmn: &str, s: f64, from: f64, to: f64) -> Self {
        let dataset = self.dataset.scale_column_in_its_range(colmn, s, from, to);
        Self {
            dataset,
            meta: self.meta,
        }
    }

    pub fn scale_column_in_range_by_reference(
        self,
        target_col: &str,
        reference_col: &str,
        s: f64,
        from: f64,
        to: f64,
    ) -> Self {
        let dataset =
            self.dataset
                .scale_column_in_range_by_reference(target_col, reference_col, s, from, to);
        Self {
            dataset,
            meta: self.meta,
        }
    }

    /// move_time_to_zero API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn move_time_to_zero(self) -> Result<Self, TGADomainError> {
        let dataset = self.dataset.move_time_to_zero()?;
        Ok(Self {
            dataset,
            meta: self.meta,
        })
    }

    /// offset_column API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn offset_column(self, colmn: &str, offset: f64) -> Self {
        let dataset = self.dataset.offset_column(colmn, offset);
        Self {
            dataset,
            meta: self.meta,
        }
    }

    pub fn offset_column_in_its_range(self, colmn: &str, offset: f64, from: f64, to: f64) -> Self {
        let dataset = self
            .dataset
            .offset_column_in_its_range(colmn, offset, from, to);
        Self {
            dataset,
            meta: self.meta,
        }
    }

    pub fn offset_column_in_range_by_reference(
        self,
        target_col: &str,
        reference_col: &str,
        offset: f64,
        from: f64,
        to: f64,
    ) -> Self {
        let dataset = self.dataset.offset_column_in_range_by_reference(
            target_col,
            reference_col,
            offset,
            from,
            to,
        );
        Self {
            dataset,
            meta: self.meta,
        }
    }
    /// calibrate_mass_from_voltage API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn calibrate_mass_from_voltage(self, k: f64, b: f64) -> Self {
        let dataset = self.dataset.calibrate_mass_from_voltage(k, b);
        Self {
            dataset,
            meta: self.meta,
        }
    }

    /// calibrate_mass API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn calibrate_mass(self, k: f64, b: f64, new_col: &str) -> Result<Self, TGADomainError> {
        let dataset = self.dataset.calibrate_mass(k, b, new_col)?;
        Ok(Self {
            dataset,
            meta: self.meta,
        })
    }

    /// unary_column_op API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn unary_column_op(
        self,
        src: &str,
        dst: &str,
        op: UnaryOp,
    ) -> Result<Self, TGADomainError> {
        let dataset = self.dataset.unary_column_op(src, dst, op)?;
        Ok(Self {
            dataset,
            meta: self.meta,
        })
    }

    /// exp_column API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn exp_column(self, col_name: &str) -> Result<Self, TGADomainError> {
        let dataset = self.dataset.exp_column(col_name)?;
        Ok(Self {
            dataset,
            meta: self.meta,
        })
    }

    /// ln_column API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn ln_column(self, col_name: &str) -> Result<Self, TGADomainError> {
        let dataset = self.dataset.ln_column(col_name)?;
        Ok(Self {
            dataset,
            meta: self.meta,
        })
    }

    //======================================================================
    // DIMENSIONLESS

    /// derive_dimensionless_mass API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn derive_dimensionless_mass(self, from: f64, to: f64) -> Result<Self, TGADomainError> {
        let dataset = self.dataset.derive_dimensionless_mass(from, to)?;
        Ok(Self {
            dataset,
            meta: self.meta,
        })
    }

    /// dimensionless_mass API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn dimensionless_mass(
        self,
        from: f64,
        to: f64,
        new_col: &str,
    ) -> Result<Self, TGADomainError> {
        let dataset = self.dataset.dimensionless_mass(from, to, new_col)?;
        Ok(Self {
            dataset,
            meta: self.meta,
        })
    }
    /// derive_conversion API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn derive_conversion(self) -> Self {
        let dataset = self.dataset.derive_conversion();
        Self {
            dataset,
            meta: self.meta,
        }
    }
    /// conversion API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn conversion(self, from: f64, to: f64, new_col: &str) -> Result<Self, TGADomainError> {
        let dataset = self.dataset.conversion(from, to, new_col)?;
        Ok(Self {
            dataset,
            meta: self.meta,
        })
    }

    //================================================================================================
    //      DIAGNOSTICS AND TESTING

    pub fn create_from_synthetic_data(
        virtga: &VirtualTGA,
        meta: ExperimentMeta,
    ) -> Result<Self, TGADomainError> {
        let dataset = TGADataset::create_from_synthetic_data(virtga)?;
        Ok(Self { dataset, meta })
    }
    /// check_nulls API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn check_nulls(self) -> Self {
        let dataset = self.dataset.check_nulls();
        Self {
            dataset,
            meta: self.meta,
        }
    }

    /// check_nulls_for_operation API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn check_nulls_for_operation(self, operation: &str) -> Self {
        let dataset = self.dataset.check_nulls_for_operation(operation);
        Self {
            dataset,
            meta: self.meta,
        }
    }

    pub fn check_nulls_for_operation_borrowed(&self, operation: &str) {
        self.dataset.check_nulls_for_operation_borrowed(operation);
    }
    /// assert_no_nulls API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn assert_no_nulls(self, operation: &str) -> Result<Self, TGADomainError> {
        let dataset = self.dataset.assert_no_nulls(operation)?;
        Ok(Self {
            dataset,
            meta: self.meta,
        })
    }

    pub fn sample_all_columns_even_layers_table(
        &self,
        n: usize,
        headers: Vec<Option<String>>,
    ) -> Result<(), TGADomainError> {
        self.dataset
            .sample_all_columns_even_layers_table(n, headers)
    }
    //================================================================================================
    // RATES
    /// derive_rate0 API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn derive_rate0(
        self,
        source_col: &str,
        new_col: &str,
        out_meta: ColumnMeta,
    ) -> Result<Self, TGADomainError> {
        let dataset = self.dataset.derive_rate0(source_col, new_col, out_meta)?;
        Ok(Self {
            dataset,
            meta: self.meta,
        })
    }

    /// derive_rate API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn derive_rate(
        self,
        source_col: &str,
        new_col: &str,
        new_col_nature: ColumnNature,
    ) -> Result<Self, TGADomainError> {
        let dataset = self
            .dataset
            .derive_rate(source_col, new_col, new_col_nature)?;
        Ok(Self {
            dataset,
            meta: self.meta,
        })
    }

    /// derive_mass_rate API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn derive_mass_rate(self, new_col: &str) -> Result<Self, TGADomainError> {
        let dataset = self.dataset.derive_mass_rate(new_col)?;
        Ok(Self {
            dataset,
            meta: self.meta,
        })
    }

    /// derive_temperature_rate API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn derive_temperature_rate(self, new_col: &str) -> Result<Self, TGADomainError> {
        let dataset = self.dataset.derive_temperature_rate(new_col)?;
        Ok(Self {
            dataset,
            meta: self.meta,
        })
    }

    /// derive_dimensionless_rate API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn derive_dimensionless_rate(
        self,
        col_name: &str,
        out_name: &str,
    ) -> Result<Self, TGADomainError> {
        let dataset = self.dataset.derive_dimensionless_rate(col_name, out_name)?;
        Ok(Self {
            dataset,
            meta: self.meta,
        })
    }

    /// derive_deta_dt API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn derive_deta_dt(self, out_name: &str) -> Result<Self, TGADomainError> {
        let dataset = self.dataset.derive_deta_dt(out_name)?;
        Ok(Self {
            dataset,
            meta: self.meta,
        })
    }

    /// derive_dalpha_dt API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn derive_dalpha_dt(self, out_name: &str) -> Result<Self, TGADomainError> {
        let dataset = self.dataset.derive_dalpha_dt(out_name)?;
        Ok(Self {
            dataset,
            meta: self.meta,
        })
    }
    //==================================================================
    //SMOOTHING AND FILTERING
    /// smooth_columns API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn smooth_columns(
        self,
        cols: &[&str],
        strategy: crate::Kinetics::experimental_kinetics::exp_kinetics_smooth_filter::SmoothStrategy,
    ) -> Result<Self, TGADomainError> {
        let dataset = self.dataset.smooth_columns(cols, strategy)?;
        Ok(Self {
            dataset,
            meta: self.meta,
        })
    }

    /// rolling_mean API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn rolling_mean(self, col_name: &str, window: usize) -> Self {
        let dataset = self.dataset.rolling_mean(col_name, window);
        Self {
            dataset,
            meta: self.meta,
        }
    }
    pub fn rolling_mean_as(self, col_name: &str, window: usize, out_col: Option<&str>) -> Self {
        let dataset = self.dataset.rolling_mean_as(col_name, window, out_col);
        Self {
            dataset,
            meta: self.meta,
        }
    }

    /// hampel_filter API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn hampel_filter(
        self,
        col: &str,
        window: usize,
        n_sigma: f64,
        strategy: crate::Kinetics::experimental_kinetics::exp_kinetics_smooth_filter::HampelStrategy,
    ) -> Result<Self, TGADomainError> {
        let dataset = self.dataset.hampel_filter(col, window, n_sigma, strategy)?;
        Ok(Self {
            dataset,
            meta: self.meta,
        })
    }
    pub fn hampel_filter_as(
        self,
        col: &str,
        window: usize,
        n_sigma: f64,
        strategy: crate::Kinetics::experimental_kinetics::exp_kinetics_smooth_filter::HampelStrategy,
        out_col: Option<&str>,
    ) -> Result<Self, TGADomainError> {
        let dataset = self
            .dataset
            .hampel_filter_as(col, window, n_sigma, strategy, out_col)?;
        Ok(Self {
            dataset,
            meta: self.meta,
        })
    }
    /// hampel_filter_null_safe API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn hampel_filter_null_safe(
        self,
        col: &str,
        window: usize,
        n_sigma: f64,
        strategy: crate::Kinetics::experimental_kinetics::exp_kinetics_smooth_filter::HampelStrategy,
    ) -> Result<Self, TGADomainError> {
        let dataset = self
            .dataset
            .hampel_filter_null_safe(col, window, n_sigma, strategy)?;
        Ok(Self {
            dataset,
            meta: self.meta,
        })
    }

    pub fn hampel_filter_null_safe_as(
        self,
        col: &str,
        window: usize,
        n_sigma: f64,
        strategy: crate::Kinetics::experimental_kinetics::exp_kinetics_smooth_filter::HampelStrategy,
        out_col: Option<&str>,
    ) -> Result<Self, TGADomainError> {
        let dataset = self
            .dataset
            .hampel_filter_null_safe_as(col, window, n_sigma, strategy, out_col)?;
        Ok(Self {
            dataset,
            meta: self.meta,
        })
    }

    /// sg_filter_column API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn sg_filter_column(
        self,
        col: &str,
        window: usize,
        poly_order: usize,
        deriv: usize,
        delta: f64,
    ) -> Result<Self, TGADomainError> {
        let dataset = self
            .dataset
            .sg_filter_column(col, window, poly_order, deriv, delta)?;
        Ok(Self {
            dataset,
            meta: self.meta,
        })
    }

    pub fn sg_filter_column_as(
        self,
        col: &str,
        window: usize,
        poly_order: usize,
        deriv: usize,
        delta: f64,
        out_col: Option<&str>,
    ) -> Result<Self, TGADomainError> {
        let dataset = self
            .dataset
            .sg_filter_column_as(col, window, poly_order, deriv, delta, out_col)?;
        Ok(Self {
            dataset,
            meta: self.meta,
        })
    }
    pub fn lowess_smooth_columns(
        self,
        time_col: &str,
        columns: &[&str],
        config: LowessConfig,
    ) -> Result<Self, TGADomainError> {
        let dataset = self
            .dataset
            .lowess_smooth_columns(time_col, columns, config)?;
        Ok(Self {
            dataset,
            meta: self.meta,
        })
    }

    pub fn lowess_smooth_column(
        self,
        time_col: &str,
        column: &str,
        config: LowessConfig,
    ) -> Result<Self, TGADomainError> {
        let dataset = self
            .dataset
            .lowess_smooth_column(time_col, column, config)?;
        Ok(Self {
            dataset,
            meta: self.meta,
        })
    }

    pub fn lowess_smooth_column_as(
        self,
        time_col: &str,
        column: &str,
        out_column: Option<&str>,
        config: LowessConfig,
    ) -> Result<Self, TGADomainError> {
        let dataset = self
            .dataset
            .lowess_smooth_column_as(time_col, column, out_column, config)?;
        Ok(Self {
            dataset,
            meta: self.meta,
        })
    }

    pub fn lowess_smooth_columns_as(
        self,
        time_col: &str,
        columns: &[&str],
        out_columns: &[Option<&str>],
        config: LowessConfig,
    ) -> Result<Self, TGADomainError> {
        let dataset =
            self.dataset
                .lowess_smooth_columns_as(time_col, columns, out_columns, config)?;
        Ok(Self {
            dataset,
            meta: self.meta,
        })
    }

    pub fn spline_resample_columns_as(
        self,
        time_col: &str,
        new_time_col: &str,
        columns: &[&str],
        out_columns: &[Option<&str>],
        n_points: usize,
        kind: SplineKind,
    ) -> Result<Self, TGADomainError> {
        let dataset = self.dataset.spline_resample_columns_as(
            time_col,
            new_time_col,
            columns,
            out_columns,
            n_points,
            kind,
        )?;
        Ok(Self {
            dataset,
            meta: self.meta,
        })
    }

    pub fn spline_resample_columns(
        self,
        time_col: &str,
        new_time_col: &str,
        columns: &[&str],

        n_points: usize,
        kind: SplineKind,
    ) -> Result<Self, TGADomainError> {
        let dataset = self.dataset.spline_resample_columns(
            time_col,
            new_time_col,
            columns,
            n_points,
            kind,
        )?;
        Ok(Self {
            dataset,
            meta: self.meta,
        })
    }

    pub fn lsq_spline_resample_columns(
        self,
        time_col: &str,
        new_time_col: &str,
        columns: &[&str],
        n_points: usize,
    ) -> Result<Self, TGADomainError> {
        let dataset =
            self.dataset
                .lsq_spline_resample_columns(time_col, new_time_col, columns, n_points)?;
        Ok(Self {
            dataset,
            meta: self.meta,
        })
    }

    pub fn lsq_spline_resample_columns_as(
        mut self,
        time_col: &str,
        new_time_col: &str,
        columns: &[&str],
        out_columns: &[Option<&str>],
        n_points: usize,
        degree: usize,
        n_internal_knots: usize,
        solver: SolverKind,
    ) -> Result<Self, TGADomainError> {
        let dataset = self.dataset.lsq_spline_resample_columns_as(
            time_col,
            new_time_col,
            columns,
            out_columns,
            n_points,
            degree,
            n_internal_knots,
            solver,
        )?;
        self.dataset = dataset;
        Ok(self)
    }
    //==========================================================
    // for GUI API features
    /// set_oneframeplot_x API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn set_oneframeplot_x(self, x_col: &str) -> Result<Self, TGADomainError> {
        let dataset = self.dataset.set_oneframeplot_x(x_col)?;
        Ok(Self {
            dataset,
            meta: self.meta,
        })
    }
    pub fn oneframeplot_axis_name(&self, axis: XY) -> Result<String, TGADomainError> {
        self.dataset.oneframeplot_axis_name(axis)
    }
    /// set_oneframeplot_y API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn set_oneframeplot_y(self, y_col: &str) -> Result<Self, TGADomainError> {
        let dataset = self.dataset.set_oneframeplot_y(y_col)?;
        Ok(Self {
            dataset,
            meta: self.meta,
        })
    }

    /// sample_oneframeplot API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn sample_oneframeplot(
        &self,
        range: Option<ViewRange>,
        max_points: usize,
    ) -> Result<PlotSeries, TGADomainError> {
        self.dataset.sample_oneframeplot(range, max_points)
    }
    pub fn sample_column(
        &self,
        col_name: &str,
        range: Option<(f64, f64)>,
        n_points: usize,
    ) -> Result<Vec<f64>, TGADomainError> {
        self.dataset.sample_column(col_name, range, n_points)
    }
    /// spline_resample_oneframeplot API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn spline_resample_oneframeplot(
        self,
        new_time_col: &str,
        n_points: usize,
        kind: SplineKind,
    ) -> Result<Self, TGADomainError> {
        let dataset = self
            .dataset
            .spline_resample_oneframeplot(new_time_col, n_points, kind)?;
        Ok(Self {
            dataset,
            meta: self.meta,
        })
    }
    pub fn spline_resample_oneframeplot_as(
        self,
        new_time_col: &str,
        n_points: usize,

        kind: SplineKind,
        out_col: Option<&str>,
    ) -> Result<Self, TGADomainError> {
        let dataset =
            self.dataset
                .spline_resample_oneframeplot_as(new_time_col, n_points, kind, out_col)?;
        Ok(Self {
            dataset,
            meta: self.meta,
        })
    }

    /// with_x_or_y API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn with_x_or_y<F>(self, axis: XY, op: F) -> Result<Self, TGADomainError>
    where
        F: FnOnce(TGADataset, &str) -> Result<TGADataset, TGADomainError>,
    {
        let dataset = self.dataset.with_x_or_y(axis, op)?;
        Ok(Self {
            dataset,
            meta: self.meta,
        })
    }

    /// with_x_and_y API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn with_x_and_y<F>(self, op: F) -> Result<Self, TGADomainError>
    where
        F: FnOnce(TGADataset, &str, &str) -> Result<TGADataset, TGADomainError>,
    {
        let dataset = self.dataset.with_x_and_y(op)?;
        Ok(Self {
            dataset,
            meta: self.meta,
        })
    }

    /// cut_before_x_or_y API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn cut_before_x_or_y(self, axis: XY, start_value: f64) -> Result<Self, TGADomainError> {
        let dataset = self.dataset.cut_before_x_or_y(axis, start_value)?;
        Ok(Self {
            dataset,
            meta: self.meta,
        })
    }

    /// cut_after_x_or_y API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn cut_after_x_or_y(self, axis: XY, end_value: f64) -> Result<Self, TGADomainError> {
        let dataset = self.dataset.cut_after_x_or_y(axis, end_value)?;
        Ok(Self {
            dataset,
            meta: self.meta,
        })
    }

    /// cut_range_x_or_y API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn cut_range_x_or_y(self, axis: XY, from: f64, to: f64) -> Result<Self, TGADomainError> {
        let dataset = self.dataset.cut_range_x_or_y(axis, from, to)?;
        Ok(Self {
            dataset,
            meta: self.meta,
        })
    }
    pub fn cut_range_inverse_x_or_y(
        self,
        axis: XY,
        from: f64,
        to: f64,
    ) -> Result<Self, TGADomainError> {
        let dataset = self.dataset.cut_range_inverse_x_or_y(axis, from, to)?;
        Ok(Self {
            dataset,
            meta: self.meta,
        })
    }
    /// min_distance_to_oneframeplot_point API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn min_distance_to_oneframeplot_point(
        &self,
        point: (f64, f64),
    ) -> Result<f64, TGADomainError> {
        self.dataset.min_distance_to_oneframeplot_point(point)
    }

    /// sample_columns API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn sample_columns(
        &self,
        time_col: &str,
        value_cols: &[&str],
        range: Option<ViewRange>,
        max_points: usize,
    ) -> Result<Vec<PlotSeries>, TGADomainError> {
        self.dataset
            .sample_columns(time_col, value_cols, range, max_points)
    }

    /// list_of_columns API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn list_of_columns(&self) -> Vec<String> {
        self.dataset.list_of_columns()
    }

    pub fn get_axis_as_vec(&self, axis: XY) -> Result<Vec<f64>, TGADomainError> {
        self.dataset.get_axis_as_vec(axis)
    }
    pub fn get_x_as_vec(&self) -> Result<Vec<f64>, TGADomainError> {
        self.dataset.get_x_as_vec()
    }
    pub fn get_y_as_vec(&self) -> Result<Vec<f64>, TGADomainError> {
        self.dataset.get_y_as_vec()
    }
    pub fn get_plotseries(&self) -> Result<PlotSeries, TGADomainError> {
        self.dataset.get_plotseries()
    }
    pub fn plot_xy_ranges(&self) -> Result<Ranges, TGADomainError> {
        self.dataset.plot_xy_ranges()
    }

    pub fn get_time_col(&self) -> Result<String, TGADomainError> {
        self.dataset.get_time_col()
    }
    pub fn get_mass_col(&self) -> Result<String, TGADomainError> {
        self.dataset.get_mass_col()
    }
    pub fn get_temperature_col(&self) -> Result<String, TGADomainError> {
        self.dataset.get_temperature_col()
    }

    pub fn offset_y_column_in_range_by_x_reference(
        self,
        offset: f64,
        from: f64,
        to: f64,
    ) -> Result<Self, TGADomainError> {
        let dataset = self
            .dataset
            .offset_y_column_in_range_by_x_reference(offset, from, to)?;
        Ok(Self {
            dataset,
            meta: self.meta,
        })
    }
    pub fn offset_y_column_in_its_range(
        self,
        offset: f64,
        from: f64,
        to: f64,
    ) -> Result<Self, TGADomainError> {
        let dataset = self
            .dataset
            .offset_y_column_in_its_range(offset, from, to)?;
        Ok(Self {
            dataset,
            meta: self.meta,
        })
    }

    pub fn offset_x_column_in_its_range(
        self,
        offset: f64,
        from: f64,
        to: f64,
    ) -> Result<Self, TGADomainError> {
        let dataset = self
            .dataset
            .offset_x_column_in_its_range(offset, from, to)?;
        Ok(Self {
            dataset,
            meta: self.meta,
        })
    }
    pub fn scale_y_column_in_range_by_x_reference(
        self,
        s: f64,
        from: f64,
        to: f64,
    ) -> Result<Self, TGADomainError> {
        let dataset = self
            .dataset
            .scale_y_column_in_range_by_x_reference(s, from, to)?;
        Ok(Self {
            dataset,
            meta: self.meta,
        })
    }
    pub fn scale_y_column_in_its_range(
        self,
        s: f64,
        from: f64,
        to: f64,
    ) -> Result<Self, TGADomainError> {
        let dataset = self.dataset.scale_y_column_in_its_range(s, from, to)?;
        Ok(Self {
            dataset,
            meta: self.meta,
        })
    }
    pub fn scale_x_column_in_its_range(
        self,
        s: f64,
        from: f64,
        to: f64,
    ) -> Result<Self, TGADomainError> {
        let dataset = self.dataset.scale_x_column_in_its_range(s, from, to)?;
        Ok(Self {
            dataset,
            meta: self.meta,
        })
    }

    pub fn history_of_operations(&self) -> Vec<OperationRecord> {
        self.dataset
            .history_of_operations
            .vector_of_operations
            .clone()
    }
    pub fn operations_on_column(&self, col: &str) -> Vec<OperationRecord> {
        self.dataset.operations_on_column(col)
    }
    pub fn get_column_history(&self, col: &str) -> ColumnHistory {
        self.dataset.get_column_history(col)
    }

    /// Compute mean of `value_col` in the time interval `[from, to]`.
    pub fn mean_on_interval(
        &self,
        value_col: &str,
        time_col: &str,
        from: f64,
        to: f64,
    ) -> Result<f64, TGADomainError> {
        Ok(self
            .dataset
            .mean_on_interval(value_col, time_col, from, to)?)
    }

    /// Compute mean on a column constrained by the same column's range.
    pub fn mean_on_interval_on_own_range(
        &self,
        column: &str,
        from: f64,
        to: f64,
    ) -> Result<f64, TGADomainError> {
        Ok(self
            .dataset
            .mean_on_interval_on_own_range(column, from, to)?)
    }

    /// Compute mean of all entries in a column.
    pub fn mean_on_column(&self, column: &str) -> Result<f64, TGADomainError> {
        Ok(self.dataset.mean_on_column(column)?)
    }

    pub fn take_column(&mut self, column_name: &str) -> Option<String> {
        self.dataset.take_column(column_name)
    }
    pub fn list_of_columns_to_recalc(&mut self) -> Vec<String> {
        self.dataset.list_of_columns_to_recalc()
    }
    pub fn drop_nulls(&mut self) -> Result<(), TGADomainError> {
        self.dataset.drop_nulls()
    }
    pub fn apply_golden_pipeline(
        self,
        config: GoldenPipelineConfig,
    ) -> Result<(Self, Vec<String>), TGADomainError> {
        let (dataset, vec_of_new) = self.dataset.apply_golden_pipeline(config)?;
        let exp = Self {
            dataset,
            meta: self.meta.clone(),
        };
        Ok((exp, vec_of_new))
    }

    pub fn get_column_by_nature(&self, nature: ColumnNature) -> Option<String> {
        self.dataset.get_column_by_nature(nature)
    }
    pub fn get_columns_by_nature(&self, nature: Vec<ColumnNature>) -> Vec<Option<String>> {
        self.dataset.get_columns_by_nature(nature)
    }

    //===============================================================================
    //  END OF THIN WRAPPERS
}

impl TGAExperiment {
    /// validate_for_kinetics API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn validate_for_kinetics(&self) -> Result<(), TGADomainError> {
        self.dataset
            .validate_required_columns(&["temperature", "alpha", "dalpha_dt"])
    }
}
//======================================================================================
//      TGA SERIES
//======================================================================================
#[derive(Clone, Debug)]
/// TGASeries data structure for the experimental kinetics API layer.
/// It groups domain data and exposes a stable facade over
/// dataset-level processing modules used throughout KiThe.
pub struct TGASeries {
    pub experiments: Vec<TGAExperiment>,
    pub exp_map: HashMap<String, usize>,
}

#[derive(Clone, Debug, Default)]
struct SeriesExperimentRecord {
    id: String,
    heating_rate: Option<f64>,
    isothermal_temperature: Option<f64>,
    comment: Option<String>,
    len: usize,
    columns: Vec<SeriesColumnRecord>,
    binds: HashMap<String, Option<String>>,
}

#[derive(Clone, Debug)]
struct SeriesColumnRecord {
    column_name: String,
    header_name: String,
    unit: Unit,
    origin: ColumnOrigin,
}

impl TGASeries {
    /// new API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn new() -> Self {
        Self {
            experiments: Vec::new(),
            exp_map: HashMap::new(),
        }
    }

    /// push_from_file API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn push_from_file(&mut self, path: &Path) -> Result<(), TGADomainError> {
        let exp = TGAExperiment::from_csv_universal(path)?;
        self.push(exp);
        Ok(())
    }

    pub fn rebuild_index(&mut self) {
        self.exp_map.clear();
        for (idx, exp) in self.experiments.iter().enumerate() {
            self.exp_map.insert(exp.meta.id.clone(), idx);
        }
    }

    pub fn drop_experiment(&mut self, idx: usize) {
        self.experiments.remove(idx);
        self.rebuild_index();
    }

    pub fn drop_experiment_by_id(&mut self, id: &str) -> Result<(), TGADomainError> {
        let idx = self.exp_map.get(id).ok_or(TGADomainError::NotImplemented)?;
        self.experiments.remove(*idx);
        self.rebuild_index();
        Ok(())
    }
    /// push API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn push(&mut self, exp: TGAExperiment) {
        let idx = self.experiments.len();
        self.exp_map.insert(exp.meta.id.clone(), idx);
        self.experiments.push(exp);
    }

    /// len API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn len(&self) -> usize {
        self.experiments.len()
    }

    /// is_empty API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn is_empty(&self) -> bool {
        self.experiments.is_empty()
    }

    /// experiments API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn experiments(&self) -> &[TGAExperiment] {
        &self.experiments
    }

    /// experiments_mut API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn experiments_mut(&mut self) -> &mut [TGAExperiment] {
        &mut self.experiments
    }

    /// get_experiment API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn get_experiment(&self, index: usize) -> Option<&TGAExperiment> {
        self.experiments.get(index)
    }

    /// get_experiment_mut API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn get_experiment_mut(&mut self, index: usize) -> Option<&mut TGAExperiment> {
        self.experiments.get_mut(index)
    }

    /// get_experiment_by_id API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn get_experiment_by_id(&self, id: &str) -> Result<&TGAExperiment, TGADomainError> {
        let idx = self.exp_map.get(id).ok_or_else(|| {
            TGADomainError::InvalidOperation(format!("Experiment id '{}' not found", id))
        })?;
        self.experiments.get(*idx).ok_or_else(|| {
            TGADomainError::InvalidOperation(format!(
                "Experiment index '{}' for id '{}' is out of bounds",
                idx, id
            ))
        })
    }

    /// get_experiment_by_id_mut API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn get_experiment_by_id_mut(
        &mut self,
        id: &str,
    ) -> Result<&mut TGAExperiment, TGADomainError> {
        let idx = *self.exp_map.get(id).ok_or_else(|| {
            TGADomainError::InvalidOperation(format!("Experiment id '{}' not found", id))
        })?;
        self.experiments.get_mut(idx).ok_or_else(|| {
            TGADomainError::InvalidOperation(format!(
                "Experiment index '{}' for id '{}' is out of bounds",
                idx, id
            ))
        })
    }

    pub fn index_by_id(&self, id: &str) -> Result<usize, TGADomainError> {
        let idx = *self.exp_map.get(id).ok_or_else(|| {
            TGADomainError::InvalidOperation(format!("Experiment id '{}' not found", id))
        })?;
        if idx >= self.experiments.len() {
            return Err(TGADomainError::InvalidOperation(format!(
                "Experiment index '{}' for id '{}' is out of bounds",
                idx, id
            )));
        }
        Ok(idx)
    }

    /// set_experiment_id API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn set_experiment_id(
        &mut self,
        index: usize,
        new_id: impl Into<String>,
    ) -> Result<(), TGADomainError> {
        let exp = self.experiments.get_mut(index).ok_or_else(|| {
            TGADomainError::InvalidOperation(format!(
                "Experiment index '{}' is out of bounds",
                index
            ))
        })?;
        exp.meta.id = new_id.into();
        self.rebuild_index();
        Ok(())
    }

    /// remove_experiment_by_id API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn remove_experiment_by_id(&mut self, id: &str) -> Result<TGAExperiment, TGADomainError> {
        let idx = *self.exp_map.get(id).ok_or_else(|| {
            TGADomainError::InvalidOperation(format!("Experiment id '{}' not found", id))
        })?;
        let removed = self.experiments.remove(idx);
        self.rebuild_index();
        Ok(removed)
    }

    /// ids API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn ids(&self) -> Vec<String> {
        self.experiments
            .iter()
            .map(|exp| exp.meta.id.clone())
            .collect()
    }

    /// get_list_of_complex_id API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn get_list_of_complex_id(&self) -> Vec<String> {
        let mut vec_of_complex_ids = Vec::new();
        for experiment in &self.experiments {
            let vec_of_columns_for_one_exp = experiment.list_of_columns();
            let experiment_id = experiment.meta.id.clone();
            for column_name in vec_of_columns_for_one_exp {
                let complex_name = format!("{}_{}", experiment_id, column_name);
                vec_of_complex_ids.push(complex_name);
            }
        }
        vec_of_complex_ids
    }
    //==========================================================================
    //                    I/O
    //==========================================================================
    /// Save the full series to a single CSV file with `#META` helper lines.
    /// The data columns are written as `<experiment_id>_<column_name>` (sanitized and uniquified).
    pub fn to_csv_series(&self, path: &Path) -> Result<(), TGADomainError> {
        let mut collected: Vec<DataFrame> = Vec::with_capacity(self.experiments.len());
        let mut max_height = 0usize;
        for exp in &self.experiments {
            let df = exp.dataset.frame.clone().collect()?;
            max_height = max_height.max(df.height());
            collected.push(df);
        }

        let mut meta_lines: Vec<String> = Vec::new();
        let mut export_columns: Vec<Column> = Vec::new();
        let mut header_counts: HashMap<String, usize> = HashMap::new();

        for (idx, exp) in self.experiments.iter().enumerate() {
            let df = &collected[idx];
            let exp_id = exp.meta.id.clone();
            let exp_id_token = Self::escape_meta_field(&exp_id);

            meta_lines.push(format!(
                "#META\tEXP\t{}\t{}\t{}\t{}",
                exp_id_token,
                Self::encode_opt_f64(exp.meta.heating_rate),
                Self::encode_opt_f64(exp.meta.isothermal_temperature),
                Self::encode_opt_string(exp.meta.comment.as_deref()),
            ));
            meta_lines.push(format!("#META\tLEN\t{}\t{}", exp_id_token, df.height()));

            for (bind_key, bind_value) in [
                ("time", exp.dataset.schema.time.as_ref()),
                ("temperature", exp.dataset.schema.temperature.as_ref()),
                ("mass", exp.dataset.schema.mass.as_ref()),
                ("alpha", exp.dataset.schema.alpha.as_ref()),
                ("dm_dt", exp.dataset.schema.dm_dt.as_ref()),
                ("eta", exp.dataset.schema.eta.as_ref()),
                ("deta_dt", exp.dataset.schema.deta_dt.as_ref()),
                ("dalpha_dt", exp.dataset.schema.dalpha_dt.as_ref()),
                ("dT_dt", exp.dataset.schema.dT_dt.as_ref()),
            ] {
                meta_lines.push(format!(
                    "#META\tBIND\t{}\t{}\t{}",
                    exp_id_token,
                    bind_key,
                    Self::encode_opt_string(bind_value.map(String::as_str)),
                ));
            }

            for col in df.columns() {
                let col_name = col.name().to_string();
                let header_name =
                    Self::build_unique_series_header(&mut header_counts, &exp_id, &col_name, idx);
                let meta = exp.dataset.schema.columns.get(&col_name);
                let unit = meta.map(|m| m.unit).unwrap_or(Unit::Unknown);
                let origin = meta.map(|m| m.origin).unwrap_or(ColumnOrigin::Imported);

                meta_lines.push(format!(
                    "#META\tCOL\t{}\t{}\t{}\t{}\t{}",
                    exp_id_token,
                    Self::escape_meta_field(&col_name),
                    Self::escape_meta_field(&header_name),
                    Self::unit_tag(unit),
                    Self::origin_tag(origin),
                ));

                let casted = col.cast(&DataType::Float64)?;
                let values = casted.f64()?;
                let mut padded: Vec<Option<f64>> = Vec::with_capacity(max_height);
                for row in 0..max_height {
                    if row < values.len() {
                        padded.push(values.get(row));
                    } else {
                        padded.push(None);
                    }
                }
                export_columns.push(Series::new(header_name.clone().into(), padded).into());
            }
        }

        let mut file = std::fs::File::create(path).map_err(|e| {
            TGADomainError::InvalidOperation(format!(
                "Failed to create series CSV '{}': {}",
                path.display(),
                e
            ))
        })?;
        writeln!(file, "# KiThe TGA Series").map_err(|e| {
            TGADomainError::InvalidOperation(format!(
                "Failed to write series CSV '{}': {}",
                path.display(),
                e
            ))
        })?;
        writeln!(file, "# format_version=1").map_err(|e| {
            TGADomainError::InvalidOperation(format!(
                "Failed to write series CSV '{}': {}",
                path.display(),
                e
            ))
        })?;
        for line in &meta_lines {
            writeln!(file, "{line}").map_err(|e| {
                TGADomainError::InvalidOperation(format!(
                    "Failed to write series CSV '{}': {}",
                    path.display(),
                    e
                ))
            })?;
        }

        if !export_columns.is_empty() {
            let mut df = DataFrame::new(max_height, export_columns)?;
            CsvWriter::new(&mut file)
                .include_header(true)
                .finish(&mut df)?;
        }

        Ok(())
    }

    /// Read the CSV format produced by `to_csv_series` and reconstruct `TGASeries`.
    pub fn from_csv_series(path: &Path) -> Result<Self, TGADomainError> {
        let text = std::fs::read_to_string(path).map_err(|e| {
            TGADomainError::InvalidOperation(format!(
                "Failed to read series CSV '{}': {}",
                path.display(),
                e
            ))
        })?;

        let mut records: HashMap<String, SeriesExperimentRecord> = HashMap::new();
        let mut order: Vec<String> = Vec::new();
        let mut csv_data_lines: Vec<String> = Vec::new();

        for line in text.lines() {
            if let Some(meta_payload) = line.strip_prefix("#META\t") {
                Self::parse_series_meta_line(meta_payload, &mut records, &mut order)?;
                continue;
            }
            if line.trim().is_empty() || line.starts_with('#') {
                continue;
            }
            csv_data_lines.push(line.to_string());
        }

        let data_df = if csv_data_lines.is_empty() {
            None
        } else {
            let tmp = tempfile::NamedTempFile::new().map_err(|e| {
                TGADomainError::InvalidOperation(format!(
                    "Failed to create temp file while reading '{}': {}",
                    path.display(),
                    e
                ))
            })?;
            std::fs::write(tmp.path(), csv_data_lines.join("\n")).map_err(|e| {
                TGADomainError::InvalidOperation(format!(
                    "Failed to write temp CSV while reading '{}': {}",
                    path.display(),
                    e
                ))
            })?;
            let plpath = PlRefPath::try_from_path(tmp.path())?;
            Some(
                LazyCsvReader::new(plpath)
                    .with_has_header(true)
                    .finish()?
                    .collect()?,
            )
        };

        let mut series = TGASeries::new();
        for id in order {
            let rec = records.get(&id).cloned().ok_or_else(|| {
                TGADomainError::InvalidOperation(format!(
                    "Series metadata is inconsistent for experiment id '{}'",
                    id
                ))
            })?;

            let mut exp_columns: Vec<Column> = Vec::new();
            if let Some(df) = data_df.as_ref() {
                for c in &rec.columns {
                    let column = df.column(&c.header_name).map_err(|_| {
                        TGADomainError::InvalidOperation(format!(
                            "Series CSV is missing expected header '{}'",
                            c.header_name
                        ))
                    })?;
                    let mut owned = column.clone();
                    owned.rename(c.column_name.clone().into());
                    exp_columns.push(owned);
                }
            }

            let mut exp_df = if exp_columns.is_empty() {
                DataFrame::default()
            } else {
                let height = exp_columns.iter().map(|c| c.len()).max().unwrap_or(0);
                DataFrame::new(height, exp_columns)?
            };
            if rec.len < exp_df.height() {
                exp_df = exp_df.slice(0, rec.len);
            }
            if rec.len > exp_df.height() {
                return Err(TGADomainError::InvalidOperation(format!(
                    "Series CSV for experiment '{}' has less data rows than declared length {}",
                    id, rec.len
                )));
            }

            let mut schema_columns: HashMap<String, ColumnMeta> = HashMap::new();
            for c in &rec.columns {
                schema_columns.insert(
                    c.column_name.clone(),
                    ColumnMeta {
                        name: c.column_name.clone(),
                        unit: c.unit,
                        origin: c.origin,
                        nature: ColumnNature::Unknown,
                    },
                );
            }

            let dataset = TGADataset {
                frame: exp_df.lazy(),
                schema: TGASchema {
                    columns: schema_columns,
                    time: rec.binds.get("time").cloned().unwrap_or(None),
                    temperature: rec.binds.get("temperature").cloned().unwrap_or(None),
                    mass: rec.binds.get("mass").cloned().unwrap_or(None),
                    alpha: rec.binds.get("alpha").cloned().unwrap_or(None),
                    dm_dt: rec.binds.get("dm_dt").cloned().unwrap_or(None),
                    eta: rec.binds.get("eta").cloned().unwrap_or(None),
                    deta_dt: rec.binds.get("deta_dt").cloned().unwrap_or(None),
                    dalpha_dt: rec.binds.get("dalpha_dt").cloned().unwrap_or(None),
                    dT_dt: rec.binds.get("dT_dt").cloned().unwrap_or(None),
                },
                oneframeplot: None,
                history_of_operations: History {
                    vector_of_operations: Vec::new(),
                    columns_has_changed: Vec::new(),
                },
            };
            let experiment = TGAExperiment {
                dataset,
                meta: ExperimentMeta {
                    id: rec.id,
                    heating_rate: rec.heating_rate,
                    isothermal_temperature: rec.isothermal_temperature,
                    comment: rec.comment,
                },
            };
            series.push(experiment);
        }

        Ok(series)
    }

    fn parse_series_meta_line(
        payload: &str,
        records: &mut HashMap<String, SeriesExperimentRecord>,
        order: &mut Vec<String>,
    ) -> Result<(), TGADomainError> {
        let parts: Vec<&str> = payload.split('\t').collect();
        if parts.is_empty() {
            return Ok(());
        }

        match parts[0] {
            "EXP" => {
                if parts.len() != 5 {
                    return Err(TGADomainError::InvalidOperation(format!(
                        "Invalid EXP metadata line: '{}'",
                        payload
                    )));
                }
                let id = Self::unescape_meta_field(parts[1])?;
                if !records.contains_key(&id) {
                    order.push(id.clone());
                }
                let rec = records
                    .entry(id.clone())
                    .or_insert_with(|| SeriesExperimentRecord {
                        id: id.clone(),
                        ..Default::default()
                    });
                rec.id = id;
                rec.heating_rate = Self::decode_opt_f64(parts[2])?;
                rec.isothermal_temperature = Self::decode_opt_f64(parts[3])?;
                rec.comment = Self::decode_opt_string(parts[4])?;
            }
            "LEN" => {
                if parts.len() != 3 {
                    return Err(TGADomainError::InvalidOperation(format!(
                        "Invalid LEN metadata line: '{}'",
                        payload
                    )));
                }
                let id = Self::unescape_meta_field(parts[1])?;
                if !records.contains_key(&id) {
                    order.push(id.clone());
                }
                let rec = records
                    .entry(id.clone())
                    .or_insert_with(|| SeriesExperimentRecord {
                        id,
                        ..Default::default()
                    });
                rec.len = parts[2].parse::<usize>().map_err(|e| {
                    TGADomainError::InvalidOperation(format!(
                        "Invalid LEN value '{}' in metadata: {}",
                        parts[2], e
                    ))
                })?;
            }
            "BIND" => {
                if parts.len() != 4 {
                    return Err(TGADomainError::InvalidOperation(format!(
                        "Invalid BIND metadata line: '{}'",
                        payload
                    )));
                }
                let id = Self::unescape_meta_field(parts[1])?;
                if !records.contains_key(&id) {
                    order.push(id.clone());
                }
                let rec = records
                    .entry(id.clone())
                    .or_insert_with(|| SeriesExperimentRecord {
                        id,
                        ..Default::default()
                    });
                rec.binds
                    .insert(parts[2].to_string(), Self::decode_opt_string(parts[3])?);
            }
            "COL" => {
                if parts.len() != 6 {
                    return Err(TGADomainError::InvalidOperation(format!(
                        "Invalid COL metadata line: '{}'",
                        payload
                    )));
                }
                let id = Self::unescape_meta_field(parts[1])?;
                if !records.contains_key(&id) {
                    order.push(id.clone());
                }
                let rec = records
                    .entry(id.clone())
                    .or_insert_with(|| SeriesExperimentRecord {
                        id,
                        ..Default::default()
                    });
                rec.columns.push(SeriesColumnRecord {
                    column_name: Self::unescape_meta_field(parts[2])?,
                    header_name: Self::unescape_meta_field(parts[3])?,
                    unit: Self::parse_unit_tag(parts[4])?,
                    origin: Self::parse_origin_tag(parts[5])?,
                });
            }
            _ => {}
        }

        Ok(())
    }

    fn build_unique_series_header(
        header_counts: &mut HashMap<String, usize>,
        exp_id: &str,
        col_name: &str,
        exp_idx: usize,
    ) -> String {
        let left = Self::sanitize_for_header(exp_id);
        let right = Self::sanitize_for_header(col_name);
        let mut base = format!("{left}_{right}");
        if left == "unnamed" {
            base = format!("exp{exp_idx}_{right}");
        }

        let count = header_counts.entry(base.clone()).or_insert(0);
        *count += 1;
        if *count == 1 {
            base
        } else {
            format!("{}_{}", base, count)
        }
    }

    fn sanitize_for_header(raw: &str) -> String {
        let mut out = String::with_capacity(raw.len());
        for ch in raw.chars() {
            if ch.is_ascii_alphanumeric() || ch == '_' {
                out.push(ch);
            } else {
                out.push('_');
            }
        }
        let out = out.trim_matches('_').to_string();
        if out.is_empty() {
            "unnamed".to_string()
        } else {
            out
        }
    }

    fn unit_tag(unit: Unit) -> &'static str {
        match unit {
            Unit::Second => "Second",
            Unit::Hour => "Hour",
            Unit::Kelvin => "Kelvin",
            Unit::Celsius => "Celsius",
            Unit::MilliVolt => "MilliVolt",
            Unit::Milligram => "Milligram",
            Unit::MilligramPerSecond => "MilligramPerSecond",
            Unit::KelvinPerSecond => "KelvinPerSecond",
            Unit::CelsiusPerSecond => "CelsiusPerSecond",
            Unit::PerSecond => "PerSecond",
            Unit::Gram => "Gram",
            Unit::Dimensionless => "Dimensionless",
            Unit::Unknown => "Unknown",
        }
    }

    fn parse_unit_tag(tag: &str) -> Result<Unit, TGADomainError> {
        match tag {
            "Second" => Ok(Unit::Second),
            "Hour" => Ok(Unit::Hour),
            "Kelvin" => Ok(Unit::Kelvin),
            "Celsius" => Ok(Unit::Celsius),
            "MilliVolt" => Ok(Unit::MilliVolt),
            "Milligram" => Ok(Unit::Milligram),
            "MilligramPerSecond" => Ok(Unit::MilligramPerSecond),
            "KelvinPerSecond" => Ok(Unit::KelvinPerSecond),
            "CelsiusPerSecond" => Ok(Unit::CelsiusPerSecond),
            "PerSecond" => Ok(Unit::PerSecond),
            "Gram" => Ok(Unit::Gram),
            "Dimensionless" => Ok(Unit::Dimensionless),
            "Unknown" => Ok(Unit::Unknown),
            _ => Err(TGADomainError::InvalidOperation(format!(
                "Unknown unit tag '{}'",
                tag
            ))),
        }
    }

    fn origin_tag(origin: ColumnOrigin) -> &'static str {
        match origin {
            ColumnOrigin::Raw => "Raw",
            ColumnOrigin::PolarsDerived => "PolarsDerived",
            ColumnOrigin::NumericDerived => "NumericDerived",
            ColumnOrigin::Imported => "Imported",
        }
    }

    fn parse_origin_tag(tag: &str) -> Result<ColumnOrigin, TGADomainError> {
        match tag {
            "Raw" => Ok(ColumnOrigin::Raw),
            "PolarsDerived" => Ok(ColumnOrigin::PolarsDerived),
            "NumericDerived" => Ok(ColumnOrigin::NumericDerived),
            "Imported" => Ok(ColumnOrigin::Imported),
            _ => Err(TGADomainError::InvalidOperation(format!(
                "Unknown column origin tag '{}'",
                tag
            ))),
        }
    }

    fn encode_opt_f64(v: Option<f64>) -> String {
        match v {
            Some(x) => format!("1:{x}"),
            None => "0".to_string(),
        }
    }

    fn decode_opt_f64(token: &str) -> Result<Option<f64>, TGADomainError> {
        if token == "0" {
            return Ok(None);
        }
        let raw = token.strip_prefix("1:").ok_or_else(|| {
            TGADomainError::InvalidOperation(format!("Invalid optional-f64 token '{}'", token))
        })?;
        let parsed = raw.parse::<f64>().map_err(|e| {
            TGADomainError::InvalidOperation(format!(
                "Invalid optional-f64 value '{}' in token '{}': {}",
                raw, token, e
            ))
        })?;
        Ok(Some(parsed))
    }

    fn encode_opt_string(v: Option<&str>) -> String {
        match v {
            Some(x) => format!("1:{}", Self::escape_meta_field(x)),
            None => "0".to_string(),
        }
    }

    fn decode_opt_string(token: &str) -> Result<Option<String>, TGADomainError> {
        if token == "0" {
            return Ok(None);
        }
        let raw = token.strip_prefix("1:").ok_or_else(|| {
            TGADomainError::InvalidOperation(format!("Invalid optional-string token '{}'", token))
        })?;
        Ok(Some(Self::unescape_meta_field(raw)?))
    }

    fn escape_meta_field(raw: &str) -> String {
        let mut out = String::with_capacity(raw.len());
        for ch in raw.chars() {
            match ch {
                '\\' => out.push_str("\\\\"),
                '\t' => out.push_str("\\t"),
                '\n' => out.push_str("\\n"),
                '\r' => out.push_str("\\r"),
                _ => out.push(ch),
            }
        }
        out
    }

    fn unescape_meta_field(raw: &str) -> Result<String, TGADomainError> {
        let mut out = String::with_capacity(raw.len());
        let mut chars = raw.chars();
        while let Some(ch) = chars.next() {
            if ch != '\\' {
                out.push(ch);
                continue;
            }

            let next = chars.next().ok_or_else(|| {
                TGADomainError::InvalidOperation(format!(
                    "Invalid escaped metadata token '{}'",
                    raw
                ))
            })?;
            match next {
                '\\' => out.push('\\'),
                't' => out.push('\t'),
                'n' => out.push('\n'),
                'r' => out.push('\r'),
                _ => {
                    return Err(TGADomainError::InvalidOperation(format!(
                        "Unsupported escape sequence '\\{}' in metadata token '{}'",
                        next, raw
                    )));
                }
            }
        }
        Ok(out)
    }
    //==========================================================
    // GENERICS
    /// examples
    ///
    /// series.try_transform_by_id("exp1", |exp| exp.move_time_to_zero())?;
    /// let cols = series.apply_by_id("exp1", |exp| exp.list_of_columns())?;
    /// let plot = series.try_apply_by_id("exp1", |exp| exp.get_plotseries())?
    ///
    /// apply_by_id API method in the experiment/series facade.
    /// Calls an infallible read-only operation on one experiment selected by id.
    /// generic argument:  F: FnOnce(&TGAExperiment) -> R,
    pub fn apply_by_id<R, F>(&self, id: &str, op: F) -> Result<R, TGADomainError>
    where
        F: FnOnce(&TGAExperiment) -> R,
    {
        let idx = self.index_by_id(id)?;
        Ok(op(&self.experiments[idx]))
    }

    /// try_apply_by_id API method in the experiment/series facade.
    /// Calls a fallible read-only operation on one experiment selected by id.
    /// generic argument:  FnOnce(&TGAExperiment) -> Result<R, TGADomainError>,
    pub fn try_apply_by_id<R, F>(&self, id: &str, op: F) -> Result<R, TGADomainError>
    where
        F: FnOnce(&TGAExperiment) -> Result<R, TGADomainError>,
    {
        let idx = self.index_by_id(id)?;
        op(&self.experiments[idx])
    }

    /// mutate_by_id API method in the experiment/series facade.
    /// Calls an infallible mutable operation on one experiment selected by id.
    ///  generic argument:  FnOnce(&mut TGAExperiment) -> R,
    pub fn mutate_by_id<R, F>(&mut self, id: &str, op: F) -> Result<R, TGADomainError>
    where
        F: FnOnce(&mut TGAExperiment) -> R,
    {
        let idx = self.index_by_id(id)?;
        Ok(op(&mut self.experiments[idx]))
    }

    /// try_mutate_by_id API method in the experiment/series facade.
    /// Calls a fallible mutable operation on one experiment selected by id.
    ///generic argument:  FnOnce(&mut TGAExperiment) -> Result<R, TGADomainError>,
    pub fn try_mutate_by_id<R, F>(&mut self, id: &str, op: F) -> Result<R, TGADomainError>
    where
        F: FnOnce(&mut TGAExperiment) -> Result<R, TGADomainError>,
    {
        let idx = self.index_by_id(id)?;
        op(&mut self.experiments[idx])
    }

    /// transform_by_id API method in the experiment/series facade.
    /// Applies an infallible consuming transformation (`TGAExperiment -> TGAExperiment`)
    /// to one experiment selected by id.
    /// generic argument: FnOnce(TGAExperiment) -> TGAExperiment,
    pub fn transform_by_id<F>(&mut self, id: &str, op: F) -> Result<(), TGADomainError>
    where
        F: FnOnce(TGAExperiment) -> TGAExperiment,
    {
        let idx = self.index_by_id(id)?;
        let prev_id = self.experiments[idx].meta.id.clone();
        let updated = op(self.experiments[idx].clone());
        self.experiments[idx] = updated;
        if self.experiments[idx].meta.id != prev_id {
            self.rebuild_index();
        }
        Ok(())
    }

    /// try_transform_by_id API method in the experiment/series facade.
    /// Applies a fallible consuming transformation
    /// (`TGAExperiment -> Result<TGAExperiment, TGADomainError>`)
    /// to one experiment selected by id.
    /// generic argument: FnOnce(TGAExperiment) -> Result<TGAExperiment, TGADomainError>,
    pub fn try_transform_by_id<F>(&mut self, id: &str, op: F) -> Result<(), TGADomainError>
    where
        F: FnOnce(TGAExperiment) -> Result<TGAExperiment, TGADomainError>,
    {
        let idx = self.index_by_id(id)?;
        let prev_id = self.experiments[idx].meta.id.clone();
        let updated = op(self.experiments[idx].clone())?;
        self.experiments[idx] = updated;
        if self.experiments[idx].meta.id != prev_id {
            self.rebuild_index();
        }
        Ok(())
    }
    //========================================================================

    //========================================================================
    // THIN WRAPPERS AROUND TGAExperiment FUNCTIONS
    //========================================================================

    //========================================================================
    // INPUT/OUTPUT
    /// bind_time API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn bind_time(&mut self, id: &str, col: &str, unit: Unit) -> Result<(), TGADomainError> {
        self.try_transform_by_id(id, |exp| exp.bind_time(col, unit))
    }

    /// bind_temperature API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn bind_temperature(
        &mut self,
        id: &str,
        col: &str,
        unit: Unit,
    ) -> Result<(), TGADomainError> {
        self.try_transform_by_id(id, |exp| exp.bind_temperature(col, unit))
    }

    /// bind_mass API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn bind_mass(&mut self, id: &str, col: &str, unit: Unit) -> Result<(), TGADomainError> {
        self.try_transform_by_id(id, |exp| exp.bind_mass(col, unit))
    }

    /// bind_column API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn bind_column(
        &mut self,
        id: &str,
        role: ColumnRole,
        name: &str,
        unit: Unit,
    ) -> Result<(), TGADomainError> {
        self.try_transform_by_id(id, |exp| exp.bind_column(role, name, unit))
    }

    /// to_csv_with_units API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn to_csv_with_units(&self, id: &str, path: &Path) -> Result<(), TGADomainError> {
        self.try_apply_by_id(id, |exp| exp.to_csv_with_units(path))
    }

    /// to_csv_with_units_selected API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn to_csv_with_units_selected(
        &self,
        id: &str,
        path: &Path,
        columns: &[&str],
    ) -> Result<(), TGADomainError> {
        self.try_apply_by_id(id, |exp| exp.to_csv_with_units_selected(path, columns))
    }

    //================================================================
    // BASIC TRANSFORMATIONS: COLUMN CUT AND FILTER
    /// with_column_expr API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn with_column_expr(
        &mut self,
        id: &str,
        meta: ColumnMeta,
        expr: polars::prelude::Expr,
    ) -> Result<(), TGADomainError> {
        self.transform_by_id(id, |exp| exp.with_column_expr(meta, expr))
    }

    /// filter_rows API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn filter_rows(
        &mut self,
        id: &str,
        predicate: polars::prelude::Expr,
    ) -> Result<(), TGADomainError> {
        self.transform_by_id(id, |exp| exp.filter_rows(predicate))
    }

    /// filter_by_mask API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn filter_by_mask(
        &mut self,
        id: &str,
        mask: polars::prelude::Expr,
    ) -> Result<(), TGADomainError> {
        self.transform_by_id(id, |exp| exp.filter_by_mask(mask))
    }

    /// cut_interval API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn cut_interval(
        &mut self,
        id: &str,
        colmn: &str,
        from: f64,
        to: f64,
    ) -> Result<(), TGADomainError> {
        self.transform_by_id(id, |exp| exp.cut_interval(colmn, from, to))
    }

    /// cut_time_interval API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn cut_time_interval(
        &mut self,
        id: &str,
        from: f64,
        to: f64,
    ) -> Result<(), TGADomainError> {
        self.transform_by_id(id, |exp| exp.cut_time_interval(from, to))
    }

    /// cut_temperature_interval API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn cut_temperature_interval(
        &mut self,
        id: &str,
        from: f64,
        to: f64,
    ) -> Result<(), TGADomainError> {
        self.transform_by_id(id, |exp| exp.cut_temperature_interval(from, to))
    }

    /// cut_mass_interval API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn cut_mass_interval(
        &mut self,
        id: &str,
        from: f64,
        to: f64,
    ) -> Result<(), TGADomainError> {
        self.transform_by_id(id, |exp| exp.cut_mass_interval(from, to))
    }

    /// cut_before_time API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn cut_before_time(&mut self, id: &str, t_start: f64) -> Result<(), TGADomainError> {
        self.transform_by_id(id, |exp| exp.cut_before_time(t_start))
    }

    /// trim_range API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn trim_range(
        &mut self,
        id: &str,
        column: &str,
        from: f64,
        to: f64,
    ) -> Result<(), TGADomainError> {
        self.transform_by_id(id, |exp| exp.trim_range(column, from, to))
    }

    /// trim_edges API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn trim_edges(
        &mut self,
        id: &str,
        left: usize,
        right: usize,
    ) -> Result<(), TGADomainError> {
        self.transform_by_id(id, |exp| exp.trim_edges(left, right))
    }

    /// trim_null_edges API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn trim_null_edges(&mut self, id: &str) -> Result<(), TGADomainError> {
        self.transform_by_id(id, |exp| exp.trim_null_edges())
    }

    /// trim_null_edges_for_columns API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn trim_null_edges_for_columns(
        &mut self,
        id: &str,
        cols_to_trim: &[&str],
    ) -> Result<(), TGADomainError> {
        self.transform_by_id(id, |exp| exp.trim_null_edges_for_columns(cols_to_trim))
    }

    /// rename_column API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn rename_column(&mut self, id: &str, old: &str, new: &str) -> Result<(), TGADomainError> {
        self.try_transform_by_id(id, |exp| exp.rename_column(old, new))
    }

    /// drop_column API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn drop_column(&mut self, id: &str, col_name: &str) -> Result<(), TGADomainError> {
        self.try_transform_by_id(id, |exp| exp.drop_column(col_name))
    }

    //======================================================================
    // UNIT TRANSFORMATIONS
    /// celsius_to_kelvin API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn celsius_to_kelvin(&mut self, id: &str) -> Result<(), TGADomainError> {
        self.transform_by_id(id, |exp| exp.celsius_to_kelvin())
    }

    /// seconds_to_hours API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn seconds_to_hours(&mut self, id: &str) -> Result<(), TGADomainError> {
        self.transform_by_id(id, |exp| exp.seconds_to_hours())
    }

    //=======================================================================
    // ARBITRARY ALGEBRAIC TRANSFORMATIONS ON COLUMN DATA
    /// scale_columns API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn scale_columns(
        &mut self,
        id: &str,
        cols: &[&str],
        factor: f64,
    ) -> Result<(), TGADomainError> {
        self.transform_by_id(id, |exp| exp.scale_columns(cols, factor))
    }

    pub fn scale_column(&mut self, id: &str, col: &str, factor: f64) -> Result<(), TGADomainError> {
        self.transform_by_id(id, |exp| exp.scale_column(col, factor))
    }

    pub fn scale_column_in_its_range(
        &mut self,
        id: &str,
        colmn: &str,
        scale: f64,
        from: f64,
        to: f64,
    ) -> Result<(), TGADomainError> {
        self.transform_by_id(id, |exp| {
            exp.scale_column_in_its_range(colmn, scale, from, to)
        })
    }

    pub fn scale_column_in_range_by_reference(
        &mut self,
        id: &str,
        target_col: &str,
        reference_col: &str,
        scale: f64,
        from: f64,
        to: f64,
    ) -> Result<(), TGADomainError> {
        self.transform_by_id(id, |exp| {
            exp.scale_column_in_range_by_reference(target_col, reference_col, scale, from, to)
        })
    }
    /// move_time_to_zero API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn move_time_to_zero(&mut self, id: &str) -> Result<(), TGADomainError> {
        self.try_transform_by_id(id, |exp| exp.move_time_to_zero())
    }

    /// offset_column API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn offset_column(
        &mut self,
        id: &str,
        colmn: &str,
        offset: f64,
    ) -> Result<(), TGADomainError> {
        self.transform_by_id(id, |exp| exp.offset_column(colmn, offset))
    }

    pub fn offset_column_in_its_range(
        &mut self,
        id: &str,
        colmn: &str,
        offset: f64,
        from: f64,
        to: f64,
    ) -> Result<(), TGADomainError> {
        self.transform_by_id(id, |exp| {
            exp.offset_column_in_its_range(colmn, offset, from, to)
        })
    }

    pub fn offset_column_in_range_by_reference(
        &mut self,
        id: &str,
        target_col: &str,
        reference_col: &str,
        offset: f64,
        from: f64,
        to: f64,
    ) -> Result<(), TGADomainError> {
        self.transform_by_id(id, |exp| {
            exp.offset_column_in_range_by_reference(target_col, reference_col, offset, from, to)
        })
    }

    /// calibrate_mass_from_voltage API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn calibrate_mass_from_voltage(
        &mut self,
        id: &str,
        k: f64,
        b: f64,
    ) -> Result<(), TGADomainError> {
        self.transform_by_id(id, |exp| exp.calibrate_mass_from_voltage(k, b))
    }

    /// calibrate_mass API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn calibrate_mass(
        &mut self,
        id: &str,
        k: f64,
        b: f64,
        new_col: &str,
    ) -> Result<(), TGADomainError> {
        self.try_transform_by_id(id, |exp| exp.calibrate_mass(k, b, new_col))
    }

    /// unary_column_op API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn unary_column_op(
        &mut self,
        id: &str,
        src: &str,
        dst: &str,
        op: UnaryOp,
    ) -> Result<(), TGADomainError> {
        self.try_transform_by_id(id, |exp| exp.unary_column_op(src, dst, op))
    }

    /// exp_column API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn exp_column(&mut self, id: &str, col_name: &str) -> Result<(), TGADomainError> {
        self.try_transform_by_id(id, |exp| exp.exp_column(col_name))
    }

    /// ln_column API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn ln_column(&mut self, id: &str, col_name: &str) -> Result<(), TGADomainError> {
        self.try_transform_by_id(id, |exp| exp.ln_column(col_name))
    }

    //======================================================================
    // DIMENSIONLESS

    /// derive_dimensionless_mass API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn derive_dimensionless_mass(
        &mut self,
        id: &str,
        from: f64,
        to: f64,
    ) -> Result<(), TGADomainError> {
        self.try_transform_by_id(id, |exp| exp.derive_dimensionless_mass(from, to))
    }

    /// dimensionless_mass API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn dimensionless_mass(
        &mut self,
        id: &str,
        from: f64,
        to: f64,
        new_col: &str,
    ) -> Result<(), TGADomainError> {
        self.try_transform_by_id(id, |exp| exp.dimensionless_mass(from, to, new_col))
    }

    /// derive_conversion API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn derive_conversion(&mut self, id: &str) -> Result<(), TGADomainError> {
        self.transform_by_id(id, |exp| exp.derive_conversion())
    }

    /// conversion API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn conversion(
        &mut self,
        id: &str,
        from: f64,
        to: f64,
        new_col: &str,
    ) -> Result<(), TGADomainError> {
        self.try_transform_by_id(id, |exp| exp.conversion(from, to, new_col))
    }

    //================================================================================================
    //      DIAGNOSTICS AND TESTING
    pub fn create_from_synthetic_data(
        &mut self,
        virtga: &VirtualTGA,
        meta: ExperimentMeta,
    ) -> Result<(), TGADomainError> {
        let exp = TGAExperiment::create_from_synthetic_data(virtga, meta)?;
        self.push(exp);
        self.rebuild_index();
        Ok(())
    }

    /// check_nulls API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn check_nulls(&mut self, id: &str) -> Result<(), TGADomainError> {
        self.transform_by_id(id, |exp| exp.check_nulls())
    }

    /// check_nulls_for_operation API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn check_nulls_for_operation(
        &mut self,
        id: &str,
        operation: &str,
    ) -> Result<(), TGADomainError> {
        self.transform_by_id(id, |exp| exp.check_nulls_for_operation(operation))
    }

    /// assert_no_nulls API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn assert_no_nulls(&mut self, id: &str, operation: &str) -> Result<(), TGADomainError> {
        self.try_transform_by_id(id, |exp| exp.assert_no_nulls(operation))
    }

    pub fn sample_all_columns_even_layers_table(
        &mut self,
        id: &str,
        n: usize,
        headers: Vec<Option<String>>,
    ) -> Result<(), TGADomainError> {
        self.try_apply_by_id(id, |exp| {
            exp.sample_all_columns_even_layers_table(n, headers)
        })
    }

    //================================================================================================
    // RATES
    /// derive_rate0 API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn derive_rate0(
        &mut self,
        id: &str,
        source_col: &str,
        new_col: &str,
        out_meta: ColumnMeta,
    ) -> Result<(), TGADomainError> {
        self.try_transform_by_id(id, |exp| exp.derive_rate0(source_col, new_col, out_meta))
    }

    /// derive_rate API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn derive_rate(
        &mut self,
        id: &str,
        source_col: &str,
        new_col: &str,
        new_col_nature: ColumnNature,
    ) -> Result<(), TGADomainError> {
        self.try_transform_by_id(id, |exp| {
            exp.derive_rate(source_col, new_col, new_col_nature)
        })
    }

    /// derive_mass_rate API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn derive_mass_rate(&mut self, id: &str, new_col: &str) -> Result<(), TGADomainError> {
        self.try_transform_by_id(id, |exp| exp.derive_mass_rate(new_col))
    }

    /// derive_temperature_rate API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn derive_temperature_rate(
        &mut self,
        id: &str,
        new_col: &str,
    ) -> Result<(), TGADomainError> {
        self.try_transform_by_id(id, |exp| exp.derive_temperature_rate(new_col))
    }

    /// derive_dimensionless_rate API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn derive_dimensionless_rate(
        &mut self,
        id: &str,
        col_name: &str,
        out_name: &str,
    ) -> Result<(), TGADomainError> {
        self.try_transform_by_id(id, |exp| exp.derive_dimensionless_rate(col_name, out_name))
    }

    /// derive_deta_dt API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn derive_deta_dt(&mut self, id: &str, out_name: &str) -> Result<(), TGADomainError> {
        self.try_transform_by_id(id, |exp| exp.derive_deta_dt(out_name))
    }

    /// derive_dalpha_dt API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn derive_dalpha_dt(&mut self, id: &str, out_name: &str) -> Result<(), TGADomainError> {
        self.try_transform_by_id(id, |exp| exp.derive_dalpha_dt(out_name))
    }

    //==================================================================
    //SMOOTHING AND FILTERING
    /// smooth_columns API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn smooth_columns(
        &mut self,
        id: &str,
        cols: &[&str],
        strategy: crate::Kinetics::experimental_kinetics::exp_kinetics_smooth_filter::SmoothStrategy,
    ) -> Result<(), TGADomainError> {
        self.try_transform_by_id(id, |exp| exp.smooth_columns(cols, strategy))
    }

    /// rolling_mean API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn rolling_mean(
        &mut self,
        id: &str,
        col_name: &str,
        window: usize,
    ) -> Result<(), TGADomainError> {
        self.transform_by_id(id, |exp| exp.rolling_mean(col_name, window))
    }
    pub fn rolling_mean_as(
        &mut self,
        id: &str,
        col_name: &str,
        window: usize,
        out_col: Option<&str>,
    ) -> Result<(), TGADomainError> {
        self.transform_by_id(id, |exp| exp.rolling_mean_as(col_name, window, out_col))
    }

    /// hampel_filter API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn hampel_filter(
        &mut self,
        id: &str,
        col: &str,
        window: usize,
        n_sigma: f64,
        strategy: crate::Kinetics::experimental_kinetics::exp_kinetics_smooth_filter::HampelStrategy,
    ) -> Result<(), TGADomainError> {
        self.try_transform_by_id(id, |exp| exp.hampel_filter(col, window, n_sigma, strategy))
    }
    pub fn hampel_filter_as(
        &mut self,
        id: &str,
        col: &str,
        window: usize,
        n_sigma: f64,
        strategy: crate::Kinetics::experimental_kinetics::exp_kinetics_smooth_filter::HampelStrategy,
        out_col: Option<&str>,
    ) -> Result<(), TGADomainError> {
        self.try_transform_by_id(id, |exp| {
            exp.hampel_filter_as(col, window, n_sigma, strategy, out_col)
        })
    }

    /// hampel_filter_null_safe API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn hampel_filter_null_safe(
        &mut self,
        id: &str,
        col: &str,
        window: usize,
        n_sigma: f64,
        strategy: crate::Kinetics::experimental_kinetics::exp_kinetics_smooth_filter::HampelStrategy,
    ) -> Result<(), TGADomainError> {
        self.try_transform_by_id(id, |exp| {
            exp.hampel_filter_null_safe(col, window, n_sigma, strategy)
        })
    }

    pub fn hampel_filter_null_safe_as(
        &mut self,
        id: &str,
        col: &str,
        window: usize,
        n_sigma: f64,
        strategy: crate::Kinetics::experimental_kinetics::exp_kinetics_smooth_filter::HampelStrategy,
        out_col: Option<&str>,
    ) -> Result<(), TGADomainError> {
        self.try_transform_by_id(id, |exp| {
            exp.hampel_filter_null_safe_as(col, window, n_sigma, strategy, out_col)
        })
    }
    /// sg_filter_column API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn sg_filter_column(
        &mut self,
        id: &str,
        col: &str,
        window: usize,
        poly_order: usize,
        deriv: usize,
        delta: f64,
    ) -> Result<(), TGADomainError> {
        self.try_transform_by_id(id, |exp| {
            exp.sg_filter_column(col, window, poly_order, deriv, delta)
        })
    }

    pub fn sg_filter_column_as(
        &mut self,
        id: &str,
        col: &str,
        window: usize,
        poly_order: usize,
        deriv: usize,
        delta: f64,
        out_col: Option<&str>,
    ) -> Result<(), TGADomainError> {
        self.try_transform_by_id(id, |exp| {
            exp.sg_filter_column_as(col, window, poly_order, deriv, delta, out_col)
        })
    }
    pub fn lowess_smooth_columns(
        &mut self,
        id: &str,
        time_col: &str,
        columns: &[&str],
        config: LowessConfig,
    ) -> Result<(), TGADomainError> {
        self.try_transform_by_id(id, |exp| {
            exp.lowess_smooth_columns(time_col, columns, config)
        })
    }

    pub fn lowess_smooth_columns_as(
        &mut self,
        id: &str,
        time_col: &str,
        columns: &[&str],
        out_columns: &[Option<&str>],
        config: LowessConfig,
    ) -> Result<(), TGADomainError> {
        self.try_transform_by_id(id, |exp| {
            exp.lowess_smooth_columns_as(time_col, columns, out_columns, config)
        })
    }

    pub fn lowess_smooth_column(
        &mut self,
        id: &str,
        time_col: &str,
        column: &str,
        config: LowessConfig,
    ) -> Result<(), TGADomainError> {
        self.try_transform_by_id(id, |exp| exp.lowess_smooth_column(time_col, column, config))
    }

    pub fn lowess_smooth_column_as(
        &mut self,
        id: &str,
        time_col: &str,
        column: &str,
        out_column: Option<&str>,
        config: LowessConfig,
    ) -> Result<(), TGADomainError> {
        self.try_transform_by_id(id, |exp| {
            exp.lowess_smooth_column_as(time_col, column, out_column, config)
        })
    }
    pub fn spline_resample_columns(
        &mut self,
        id: &str,
        time_col: &str,
        new_time_col: &str,
        columns: &[&str],

        n_points: usize,
        kind: SplineKind,
    ) -> Result<(), TGADomainError> {
        self.try_transform_by_id(id, |exp| {
            exp.spline_resample_columns(time_col, new_time_col, columns, n_points, kind)
        })
    }
    pub fn spline_resample_columns_as(
        &mut self,
        id: &str,
        time_col: &str,
        new_time_col: &str,
        columns: &[&str],
        out_columns: &[Option<&str>],
        n_points: usize,
        kind: SplineKind,
    ) -> Result<(), TGADomainError> {
        self.try_transform_by_id(id, |exp| {
            exp.spline_resample_columns_as(
                time_col,
                new_time_col,
                columns,
                out_columns,
                n_points,
                kind,
            )
        })
    }

    pub fn lsq_spline_resample_columns(
        &mut self,
        id: &str,
        time_col: &str,
        new_time_col: &str,
        columns: &[&str],
        n_points: usize,
    ) -> Result<(), TGADomainError> {
        self.try_transform_by_id(id, |exp| {
            exp.lsq_spline_resample_columns(time_col, new_time_col, columns, n_points)
        })
    }

    pub fn lsq_spline_resample_columns_as(
        &mut self,
        id: &str,
        time_col: &str,
        new_time_col: &str,
        columns: &[&str],
        out_columns: &[Option<&str>],
        n_points: usize,
        degree: usize,
        n_internal_knots: usize,
        solver: SolverKind,
    ) -> Result<(), TGADomainError> {
        self.try_transform_by_id(id, |exp| {
            exp.lsq_spline_resample_columns_as(
                time_col,
                new_time_col,
                columns,
                out_columns,
                n_points,
                degree,
                n_internal_knots,
                solver,
            )
        })
    }
    //==========================================================
    // for GUI API features
    /// set_oneframeplot_x API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn set_oneframeplot_x(&mut self, id: &str, x_col: &str) -> Result<(), TGADomainError> {
        self.try_transform_by_id(id, |exp| exp.set_oneframeplot_x(x_col))
    }

    /// set_oneframeplot_y API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn set_oneframeplot_y(&mut self, id: &str, y_col: &str) -> Result<(), TGADomainError> {
        self.try_transform_by_id(id, |exp| exp.set_oneframeplot_y(y_col))
    }

    /// spline_resample_oneframeplot API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn spline_resample_oneframeplot(
        &mut self,
        id: &str,
        new_time_col: &str,
        n_points: usize,
        kind: SplineKind,
    ) -> Result<(), TGADomainError> {
        self.try_transform_by_id(id, |exp| {
            exp.spline_resample_oneframeplot(new_time_col, n_points, kind)
        })
    }
    pub fn spline_resample_oneframeplot_as(
        &mut self,
        id: &str,
        new_time_col: &str,
        n_points: usize,
        kind: SplineKind,
        out_col: Option<&str>,
    ) -> Result<(), TGADomainError> {
        self.try_transform_by_id(id, |exp| {
            exp.spline_resample_oneframeplot_as(new_time_col, n_points, kind, out_col)
        })
    }
    /// with_x_or_y API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn with_x_or_y<F>(&mut self, id: &str, axis: XY, op: F) -> Result<(), TGADomainError>
    where
        F: FnOnce(TGADataset, &str) -> Result<TGADataset, TGADomainError>,
    {
        self.try_transform_by_id(id, |exp| exp.with_x_or_y(axis, op))
    }

    /// with_x_and_y API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn with_x_and_y<F>(&mut self, id: &str, op: F) -> Result<(), TGADomainError>
    where
        F: FnOnce(TGADataset, &str, &str) -> Result<TGADataset, TGADomainError>,
    {
        self.try_transform_by_id(id, |exp| exp.with_x_and_y(op))
    }

    /// cut_before_x_or_y API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn cut_before_x_or_y(
        &mut self,
        id: &str,
        axis: XY,
        start_value: f64,
    ) -> Result<(), TGADomainError> {
        self.try_transform_by_id(id, |exp| exp.cut_before_x_or_y(axis, start_value))
    }

    /// cut_after_x_or_y API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn cut_after_x_or_y(
        &mut self,
        id: &str,
        axis: XY,
        end_value: f64,
    ) -> Result<(), TGADomainError> {
        self.try_transform_by_id(id, |exp| exp.cut_after_x_or_y(axis, end_value))
    }

    /// cut_range_x_or_y API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn cut_range_x_or_y(
        &mut self,
        id: &str,
        axis: XY,
        from: f64,
        to: f64,
    ) -> Result<(), TGADomainError> {
        self.try_transform_by_id(id, |exp| exp.cut_range_x_or_y(axis, from, to))
    }
    pub fn cut_range_inverse_x_or_y(
        &mut self,
        id: &str,
        axis: XY,
        from: f64,
        to: f64,
    ) -> Result<(), TGADomainError> {
        self.try_transform_by_id(id, |exp| exp.cut_range_inverse_x_or_y(axis, from, to))
    }
    pub fn sample_oneframeplot(
        &mut self,
        id: &str,
        range: Option<ViewRange>,
        max_points: usize,
    ) -> Result<PlotSeries, TGADomainError> {
        self.try_apply_by_id(id, |exp| exp.sample_oneframeplot(range, max_points))
    }
    pub fn sample_column(
        &self,
        id: &str,
        col_name: &str,
        range: Option<(f64, f64)>,
        n_points: usize,
    ) -> Result<Vec<f64>, TGADomainError> {
        self.try_apply_by_id(id, |exp| exp.sample_column(col_name, range, n_points))
    }
    pub fn list_of_columns(&self, id: &str) -> Result<Vec<String>, TGADomainError> {
        Ok(self.apply_by_id(id, |exp| exp.list_of_columns())?)
    }

    pub fn plot_xy_ranges(&self, id: &str) -> Result<Ranges, TGADomainError> {
        let r = self.apply_by_id(id, |exp| exp.plot_xy_ranges())?;
        r
    }

    pub fn get_temperature_col(&self, id: &str) -> Result<String, TGADomainError> {
        let r = self.apply_by_id(id, |exp| exp.get_temperature_col())?;
        r
    }

    pub fn get_mass_col(&self, id: &str) -> Result<String, TGADomainError> {
        let r = self.apply_by_id(id, |exp| exp.get_mass_col())?;
        r
    }
    pub fn get_time_col(&self, id: &str) -> Result<String, TGADomainError> {
        let r = self.apply_by_id(id, |exp| exp.get_time_col())?;
        r
    }
    pub fn oneframeplot_axis_name(&self, id: &str, axis: XY) -> Result<String, TGADomainError> {
        let r = self.apply_by_id(id, |exp| exp.oneframeplot_axis_name(axis))?;
        r
    }

    pub fn offset_y_column_in_range_by_x_reference(
        &mut self,
        id: &str,
        offset: f64,
        from: f64,
        to: f64,
    ) -> Result<(), TGADomainError> {
        self.try_transform_by_id(id, |exp| {
            exp.offset_y_column_in_range_by_x_reference(offset, from, to)
        })
    }
    pub fn offset_y_column_in_its_range(
        &mut self,
        id: &str,
        offset: f64,
        from: f64,
        to: f64,
    ) -> Result<(), TGADomainError> {
        self.try_transform_by_id(id, |exp| exp.offset_y_column_in_its_range(offset, from, to))
    }

    pub fn offset_x_column_in_its_range(
        &mut self,
        id: &str,
        offset: f64,
        from: f64,
        to: f64,
    ) -> Result<(), TGADomainError> {
        self.try_transform_by_id(id, |exp| exp.offset_x_column_in_its_range(offset, from, to))
    }

    pub fn scale_y_column_in_range_by_x_reference(
        &mut self,
        id: &str,
        scale: f64,
        from: f64,
        to: f64,
    ) -> Result<(), TGADomainError> {
        self.try_transform_by_id(id, |exp| {
            exp.scale_y_column_in_range_by_x_reference(scale, from, to)
        })
    }
    pub fn scale_y_column_in_its_range(
        &mut self,
        id: &str,
        scale: f64,
        from: f64,
        to: f64,
    ) -> Result<(), TGADomainError> {
        self.try_transform_by_id(id, |exp| exp.scale_y_column_in_its_range(scale, from, to))
    }

    pub fn scale_x_column_in_its_range(
        &mut self,
        id: &str,
        scale: f64,
        from: f64,
        to: f64,
    ) -> Result<(), TGADomainError> {
        self.try_transform_by_id(id, |exp| exp.scale_x_column_in_its_range(scale, from, to))
    }

    pub fn history_of_operations(&self, id: &str) -> Result<Vec<OperationRecord>, TGADomainError> {
        self.apply_by_id(id, |exp| exp.history_of_operations())
    }

    pub fn operations_on_column(
        &mut self,
        id: &str,
        col: &str,
    ) -> Result<Vec<OperationRecord>, TGADomainError> {
        self.apply_by_id(id, |exp| exp.operations_on_column(col).clone())
    }
    pub fn get_column_history(
        &mut self,
        id: &str,
        col: &str,
    ) -> Result<ColumnHistory, TGADomainError> {
        self.apply_by_id(id, |exp| exp.get_column_history(col))
    }

    /// series-level wrapper for `mean_on_interval`.
    pub fn mean_on_interval(
        &self,
        id: &str,
        value_col: &str,
        time_col: &str,
        from: f64,
        to: f64,
    ) -> Result<f64, TGADomainError> {
        self.try_apply_by_id(id, |exp| {
            exp.mean_on_interval(value_col, time_col, from, to)
        })
    }

    /// series-level wrapper for `mean_on_interval_on_own_range`.
    pub fn mean_on_interval_on_own_range(
        &self,
        id: &str,
        col: &str,
        from: f64,
        to: f64,
    ) -> Result<f64, TGADomainError> {
        self.try_apply_by_id(id, |exp| exp.mean_on_interval_on_own_range(col, from, to))
    }

    /// series-level wrapper for `mean_on_column`.
    pub fn mean_on_column(&self, id: &str, col: &str) -> Result<f64, TGADomainError> {
        self.try_apply_by_id(id, |exp| exp.mean_on_column(col))
    }

    pub fn take_column(
        &mut self,
        id: &str,
        column_name: &str,
    ) -> Result<Option<String>, TGADomainError> {
        self.mutate_by_id(id, |exp| exp.take_column(column_name))
    }

    pub fn list_of_columns_to_recalc(&mut self, id: &str) -> Result<Vec<String>, TGADomainError> {
        self.mutate_by_id(id, |exp| exp.list_of_columns_to_recalc())
    }
    pub fn set_heating_rate(&mut self, id: &str, rate: f64) -> Result<(), TGADomainError> {
        self.transform_by_id(id, |exp| exp.with_heating_rate(rate))
    }
    pub fn set_comment(&mut self, id: &str, comment: &str) -> Result<(), TGADomainError> {
        self.transform_by_id(id, |exp| exp.with_comment(comment))
    }

    pub fn set_experiment_temperature(&mut self, id: &str, T: f64) -> Result<(), TGADomainError> {
        self.transform_by_id(id, |exp| exp.with_isothermal_temperature(T))
    }

    pub fn get_column_by_nature(
        &self,
        id: &str,
        nature: ColumnNature,
    ) -> Result<Option<String>, TGADomainError> {
        self.apply_by_id(id, |exp| exp.get_column_by_nature(nature))
    }
    pub fn get_columns_by_nature(
        &self,
        id: &str,
        nature: Vec<ColumnNature>,
    ) -> Result<Vec<Option<String>>, TGADomainError> {
        self.apply_by_id(id, |exp| exp.get_columns_by_nature(nature))
    }
}
