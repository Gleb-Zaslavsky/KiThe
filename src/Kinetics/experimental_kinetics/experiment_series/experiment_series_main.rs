//! Публичный фасад и основная точка сборки для `TGASeries` и связанных оберток.
//!
//! Здесь остаются только общие экспорты, макросы-генераторы и high-level wrappers.
//! Более узкие блоки вынесены в `experiment_series_core.rs`, `experiment_series_io.rs`,
//! `experiment_series_ops.rs` и `experiment_series_analysis.rs`.

use crate::Kinetics::experimental_kinetics::exp_engine_api::{
    GoldenPipelineConfig, PlotSeries, Ranges, ViewRange, XY,
};
use crate::Kinetics::experimental_kinetics::one_experiment_dataset::ColumnRole;
use crate::Kinetics::experimental_kinetics::one_experiment_dataset::History;
use crate::Kinetics::experimental_kinetics::one_experiment_dataset::{
    ColumnHistory, ColumnMeta, ColumnNature, ColumnOrigin, OperationRecord, TGADataset,
    TGADomainError, TGASchema, UnaryOp, Unit,
};
use crate::Kinetics::experimental_kinetics::testing_mod::VirtualTGA;
use RustedSciThe::numerical::data_processing::LSQSplines::SolverKind;
use RustedSciThe::numerical::data_processing::lowess_wrapper::LowessConfig;
use RustedSciThe::numerical::data_processing::splines::SplineKind;
use polars::prelude::{
    Column, CsvWriter, DataFrame, DataType, IntoLazy, LazyCsvReader, LazyFileListReader, NamedFrom,
    PlRefPath, SerWriter, Series,
};
use std::collections::HashMap;
use std::io::Write;
use std::path::Path;

#[path = "experiment_series_core.rs"]
mod experiment_series_core;

#[path = "experiment_series_io.rs"]
mod experiment_series_io;

#[path = "experiment_series_ops.rs"]
mod experiment_series_ops;

#[path = "experiment_series_analysis.rs"]
mod experiment_series_analysis;

// Forwarding macro for `TGAExperiment -> Result<TGAExperiment, TGADomainError>` wrappers.
// It keeps the "call dataset method, preserve meta" boilerplate in one place.
macro_rules! tgaexperiment_forward_result {
    ($(pub fn $name:ident(self $(, $arg:ident : $ty:ty)* $(,)?) -> Result<Self, TGADomainError> => $target:ident;)+) => {
        $(
            pub fn $name(self $(, $arg: $ty)*) -> Result<Self, TGADomainError> {
                let dataset = self.dataset.$target($($arg),*)?;
                Ok(Self {
                    dataset,
                    meta: self.meta,
                })
            }
        )+
    };
}

// Forwarding macro for `TGAExperiment -> TGAExperiment` wrappers.
// These are infallible dataset transforms that only need to carry `meta` through.
macro_rules! tgaexperiment_forward_self {
    ($(pub fn $name:ident(self $(, $arg:ident : $ty:ty)* $(,)?) -> Self => $target:ident;)+) => {
        $(
            pub fn $name(self $(, $arg: $ty)*) -> Self {
                let dataset = self.dataset.$target($($arg),*);
                Self {
                    dataset,
                    meta: self.meta,
                }
            }
        )+
    };
    ($(pub fn $name:ident(self $(, $arg:ident : $ty:ty)* $(,)?) -> Result<Self, TGADomainError> => $target:ident;)+) => {
        $(
            pub fn $name(self $(, $arg: $ty)*) -> Result<Self, TGADomainError> {
                let dataset = self.dataset.$target($($arg),*)?;
                Ok(Self {
                    dataset,
                    meta: self.meta,
                })
            }
        )+
    };
}

// Series wrapper for fallible in-place mutable operations selected by experiment id.
// It routes to `try_mutate_by_id`, which expects `FnOnce(&mut TGAExperiment) -> Result<_, TGADomainError>`.
macro_rules! tgaseries_forward_try_mutate {
    ($(pub fn $name:ident(&mut self, id: &str $(, $arg:ident : $ty:ty)* $(,)?) -> Result<(), TGADomainError> => $target:ident;)+) => {
        $(
            pub fn $name(&mut self, id: &str $(, $arg: $ty)*) -> Result<(), TGADomainError> {
                self.try_mutate_by_id(id, |exp| exp.$target($($arg),*))
            }
        )+
    };
}

// Series wrapper for infallible mutable transforms selected by experiment id.
// It routes to `transform_by_id`, which expects `TGAExperiment -> TGAExperiment`.
macro_rules! tgaseries_forward_transform {
    ($(pub fn $name:ident(&mut self, id: &str $(, $arg:ident : $ty:ty)* $(,)?) -> Result<(), TGADomainError> => $target:ident;)+) => {
        $(
            pub fn $name(&mut self, id: &str $(, $arg: $ty)*) -> Result<(), TGADomainError> {
                self.transform_by_id(id, |exp| exp.$target($($arg),*))
            }
        )+
    };
}

// Series wrapper for fallible mutable transforms selected by experiment id.
// It routes to `try_transform_by_id`, which expects `TGAExperiment -> Result<TGAExperiment, TGADomainError>`.
macro_rules! tgaseries_forward_try_transform {
    ($(pub fn $name:ident(&mut self, id: &str $(, $arg:ident : $ty:ty)* $(,)?) -> Result<(), TGADomainError> => $target:ident;)+) => {
        $(
            pub fn $name(&mut self, id: &str $(, $arg: $ty)*) -> Result<(), TGADomainError> {
                self.try_transform_by_id(id, |exp| exp.$target($($arg),*))
            }
        )+
    };
}

// Series wrapper for infallible read-only operations selected by experiment id.
// It routes to `apply_by_id`, which expects `&TGAExperiment -> R`.
macro_rules! tgaseries_forward_apply {
    ($(pub fn $name:ident(&self, id: &str $(, $arg:ident : $ty:ty)* $(,)?) -> Result<$ret:ty, TGADomainError> => $target:ident;)+) => {
        $(
            pub fn $name(&self, id: &str $(, $arg: $ty)*) -> Result<$ret, TGADomainError> {
                self.apply_by_id(id, |exp| exp.$target($($arg),*))
            }
        )+
    };
}

// Series wrapper for fallible read-only operations selected by experiment id.
// It routes to `try_apply_by_id`, which expects `&TGAExperiment -> Result<R, TGADomainError>`.
macro_rules! tgaseries_forward_try_apply {
    ($(pub fn $name:ident(&self, id: &str $(, $arg:ident : $ty:ty)* $(,)?) -> Result<$ret:ty, TGADomainError> => $target:ident;)+) => {
        $(
            pub fn $name(&self, id: &str $(, $arg: $ty)*) -> Result<$ret, TGADomainError> {
                self.try_apply_by_id(id, |exp| exp.$target($($arg),*))
            }
        )+
    };
}

// Series wrapper for mutating read/write operations selected by experiment id.
// It routes to `mutate_by_id`, which expects `TGAExperiment -> R`.
macro_rules! tgaseries_forward_mutate {
    ($(pub fn $name:ident(&mut self, id: &str $(, $arg:ident : $ty:ty)* $(,)?) -> Result<$ret:ty, TGADomainError> => $target:ident;)+) => {
        $(
            pub fn $name(&mut self, id: &str $(, $arg: $ty)*) -> Result<$ret, TGADomainError> {
                self.mutate_by_id(id, |exp| exp.$target($($arg),*))
            }
        )+
    };
}

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

#[derive(Clone, Debug, Default)]
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

    /// In-place binding variant used by the series facade to avoid cloning experiments.
    pub fn bind_time_inplace(&mut self, col: &str, unit: Unit) -> Result<(), TGADomainError> {
        self.dataset.bind_time_inplace(col, unit)
    }

    /// In-place binding variant used by the series facade to avoid cloning experiments.
    pub fn bind_temperature_inplace(
        &mut self,
        col: &str,
        unit: Unit,
    ) -> Result<(), TGADomainError> {
        self.dataset.bind_temperature_inplace(col, unit)
    }

    /// In-place binding variant used by the series facade to avoid cloning experiments.
    pub fn bind_mass_inplace(&mut self, col: &str, unit: Unit) -> Result<(), TGADomainError> {
        self.dataset.bind_mass_inplace(col, unit)
    }

    /// In-place binding variant used by the series facade to avoid cloning experiments.
    pub fn bind_column_inplace(
        &mut self,
        role: ColumnRole,
        name: &str,
        unit: Unit,
    ) -> Result<(), TGADomainError> {
        self.dataset.bind_column_inplace(role, name, unit)
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
    tgaexperiment_forward_result! {
        pub fn bind_time(self, col: &str, unit: Unit) -> Result<Self, TGADomainError> => bind_time;
        pub fn bind_temperature(self, col: &str, unit: Unit) -> Result<Self, TGADomainError> => bind_temperature;
        pub fn bind_mass(self, col: &str, unit: Unit) -> Result<Self, TGADomainError> => bind_mass;
        pub fn bind_column(self, role: ColumnRole, name: &str, unit: Unit) -> Result<Self, TGADomainError> => bind_column;
        pub fn rename_column(self, old: &str, new: &str) -> Result<Self, TGADomainError> => rename_column;
        pub fn drop_column(self, col_name: &str) -> Result<Self, TGADomainError> => drop_column;
        pub fn move_time_to_zero(self) -> Result<Self, TGADomainError> => move_time_to_zero;
        pub fn sort_by_bound_time(self) -> Result<Self, TGADomainError> => sort_by_bound_time;
        pub fn calibrate_mass(self, k: f64, b: f64, new_col: &str) -> Result<Self, TGADomainError> => calibrate_mass;
    }

    tgaexperiment_forward_result! {
        pub fn with_column_expr(self, meta: ColumnMeta, expr: polars::prelude::Expr) -> Result<Self, TGADomainError> => with_column_expr;
    }

    tgaexperiment_forward_self! {
        pub fn filter_rows(self, predicate: polars::prelude::Expr) -> Self => filter_rows;
        pub fn filter_by_mask(self, mask: polars::prelude::Expr) -> Self => filter_by_mask;
        pub fn cut_interval(self, colmn: &str, from: f64, to: f64) -> Self => cut_interval;
        pub fn cut_time_interval(self, from: f64, to: f64) -> Self => cut_time_interval;
        pub fn cut_temperature_interval(self, from: f64, to: f64) -> Self => cut_temperature_interval;
        pub fn cut_mass_interval(self, from: f64, to: f64) -> Self => cut_mass_interval;
        pub fn cut_before_time(self, t_start: f64) -> Self => cut_before_time;
        pub fn trim_range(self, column: &str, from: f64, to: f64) -> Self => trim_range;
        pub fn trim_edges(self, left: usize, right: usize) -> Self => trim_edges;
        pub fn trim_null_edges(self) -> Self => trim_null_edges;
        pub fn trim_null_edges_for_columns(self, cols_to_trim: &[&str]) -> Self => trim_null_edges_for_columns;
        pub fn celsius_to_kelvin(self) -> Self => celsius_to_kelvin;
        pub fn seconds_to_hours(self) -> Self => seconds_to_hours;
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
    // check_nulls API method in the experiment/series facade.
    // It keeps call-sites ergonomic while delegating core logic
    // to dataset and engine modules in experimental kinetics.
    tgaexperiment_forward_self! {
        pub fn check_nulls(self) -> Self => check_nulls;
        pub fn check_nulls_for_operation(self, operation: &str) -> Self => check_nulls_for_operation;
    }

    pub fn check_nulls_for_operation_borrowed(&self, operation: &str) {
        self.dataset.check_nulls_for_operation_borrowed(operation);
    }
    // assert_no_nulls API method in the experiment/series facade.
    // It keeps call-sites ergonomic while delegating core logic
    // to dataset and engine modules in experimental kinetics.
    tgaexperiment_forward_result! {
        pub fn assert_no_nulls(self, operation: &str) -> Result<Self, TGADomainError> => assert_no_nulls;
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
    // derive_rate0 API method in the experiment/series facade.
    // It keeps call-sites ergonomic while delegating core logic
    // to dataset and engine modules in experimental kinetics.
    tgaexperiment_forward_result! {
        pub fn derive_rate0(self, source_col: &str, new_col: &str, out_meta: ColumnMeta) -> Result<Self, TGADomainError> => derive_rate0;
        pub fn derive_rate(self, source_col: &str, new_col: &str, new_col_nature: ColumnNature) -> Result<Self, TGADomainError> => derive_rate;
        pub fn derive_mass_rate(self, new_col: &str) -> Result<Self, TGADomainError> => derive_mass_rate;
        pub fn derive_temperature_rate(self, new_col: &str) -> Result<Self, TGADomainError> => derive_temperature_rate;
        pub fn derive_dimensionless_rate(self, col_name: &str, out_name: &str) -> Result<Self, TGADomainError> => derive_dimensionless_rate;
        pub fn derive_deta_dt(self, out_name: &str) -> Result<Self, TGADomainError> => derive_deta_dt;
        pub fn derive_dalpha_dt(self, out_name: &str) -> Result<Self, TGADomainError> => derive_dalpha_dt;
    }
    // analysis helpers moved to experiment_series_analysis.rs
    //===============================================================================
    //  END OF THIN WRAPPERS
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
impl TGASeries {
    // core/index helpers live in experiment_series_core.rs

    //==========================================================================
    //                    I/O
    //==========================================================================
    /// Save the full series to a single CSV file with `#META` helper lines.
    /// The data columns are written as `<experiment_id>_<column_name>` (sanitized and uniquified).
    pub fn to_csv_series(&self, path: &Path) -> Result<(), TGADomainError> {
        experiment_series_io::to_csv_series_impl(self, path)
    }

    /// Read the CSV format produced by `to_csv_series` and reconstruct `TGASeries`.
    pub fn from_csv_series(path: &Path) -> Result<Self, TGADomainError> {
        experiment_series_io::from_csv_series_impl(path)
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
        let experiment = std::mem::take(&mut self.experiments[idx]);
        let updated = op(experiment);
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
    // bind_time API method in the experiment/series facade.
    // It keeps call-sites ergonomic while delegating core logic
    // to dataset and engine modules in experimental kinetics.
    tgaseries_forward_try_mutate! {
        pub fn bind_time(&mut self, id: &str, col: &str, unit: Unit) -> Result<(), TGADomainError> => bind_time_inplace;
    }

    // bind_temperature API method in the experiment/series facade.
    // It keeps call-sites ergonomic while delegating core logic
    // to dataset and engine modules in experimental kinetics.
    tgaseries_forward_try_mutate! {
        pub fn bind_temperature(&mut self, id: &str, col: &str, unit: Unit) -> Result<(), TGADomainError> => bind_temperature_inplace;
    }

    // bind_mass API method in the experiment/series facade.
    // It keeps call-sites ergonomic while delegating core logic
    // to dataset and engine modules in experimental kinetics.
    tgaseries_forward_try_mutate! {
        pub fn bind_mass(&mut self, id: &str, col: &str, unit: Unit) -> Result<(), TGADomainError> => bind_mass_inplace;
    }

    // bind_column API method in the experiment/series facade.
    // It keeps call-sites ergonomic while delegating core logic
    // to dataset and engine modules in experimental kinetics.
    tgaseries_forward_try_mutate! {
        pub fn bind_column(&mut self, id: &str, role: ColumnRole, name: &str, unit: Unit) -> Result<(), TGADomainError> => bind_column_inplace;
    }

    // to_csv_with_units API method in the experiment/series facade.
    // It keeps call-sites ergonomic while delegating core logic
    // to dataset and engine modules in experimental kinetics.
    tgaseries_forward_try_apply! {
        pub fn to_csv_with_units(&self, id: &str, path: &Path) -> Result<(), TGADomainError> => to_csv_with_units;
        pub fn to_csv_with_units_selected(&self, id: &str, path: &Path, columns: &[&str]) -> Result<(), TGADomainError> => to_csv_with_units_selected;
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
        self.try_transform_by_id(id, |exp| exp.with_column_expr(meta, expr))?;
        Ok(())
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

    // cut_before_time API method in the experiment/series facade.
    // It keeps call-sites ergonomic while delegating core logic
    // to dataset and engine modules in experimental kinetics.
    tgaseries_forward_transform! {
        pub fn cut_before_time(&mut self, id: &str, t_start: f64) -> Result<(), TGADomainError> => cut_before_time;
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

    // trim_null_edges API method in the experiment/series facade.
    // It keeps call-sites ergonomic while delegating core logic
    // to dataset and engine modules in experimental kinetics.
    tgaseries_forward_transform! {
        pub fn trim_null_edges(&mut self, id: &str) -> Result<(), TGADomainError> => trim_null_edges;
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

    // rename_column API method in the experiment/series facade.
    // It keeps call-sites ergonomic while delegating core logic
    // to dataset and engine modules in experimental kinetics.
    tgaseries_forward_try_transform! {
        pub fn rename_column(&mut self, id: &str, old: &str, new: &str) -> Result<(), TGADomainError> => rename_column;
    }

    // drop_column API method in the experiment/series facade.
    // It keeps call-sites ergonomic while delegating core logic
    // to dataset and engine modules in experimental kinetics.
    tgaseries_forward_try_transform! {
        pub fn drop_column(&mut self, id: &str, col_name: &str) -> Result<(), TGADomainError> => drop_column;
    }

    //======================================================================
    // UNIT TRANSFORMATIONS
    // celsius_to_kelvin API method in the experiment/series facade.
    // It keeps call-sites ergonomic while delegating core logic
    // to dataset and engine modules in experimental kinetics.
    tgaseries_forward_transform! {
        pub fn celsius_to_kelvin(&mut self, id: &str) -> Result<(), TGADomainError> => celsius_to_kelvin;
    }

    // seconds_to_hours API method in the experiment/series facade.
    // It keeps call-sites ergonomic while delegating core logic
    // to dataset and engine modules in experimental kinetics.
    tgaseries_forward_transform! {
        pub fn seconds_to_hours(&mut self, id: &str) -> Result<(), TGADomainError> => seconds_to_hours;
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
    // move_time_to_zero API method in the experiment/series facade.
    // It keeps call-sites ergonomic while delegating core logic
    // to dataset and engine modules in experimental kinetics.
    tgaseries_forward_try_transform! {
        pub fn move_time_to_zero(&mut self, id: &str) -> Result<(), TGADomainError> => move_time_to_zero;
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

    // calibrate_mass API method in the experiment/series facade.
    // It keeps call-sites ergonomic while delegating core logic
    // to dataset and engine modules in experimental kinetics.
    tgaseries_forward_try_transform! {
        pub fn calibrate_mass(&mut self, id: &str, k: f64, b: f64, new_col: &str) -> Result<(), TGADomainError> => calibrate_mass;
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

    // exp_column API method in the experiment/series facade.
    // It keeps call-sites ergonomic while delegating core logic
    // to dataset and engine modules in experimental kinetics.
    tgaseries_forward_try_transform! {
        pub fn exp_column(&mut self, id: &str, col_name: &str) -> Result<(), TGADomainError> => exp_column;
    }

    // ln_column API method in the experiment/series facade.
    // It keeps call-sites ergonomic while delegating core logic
    // to dataset and engine modules in experimental kinetics.
    tgaseries_forward_try_transform! {
        pub fn ln_column(&mut self, id: &str, col_name: &str) -> Result<(), TGADomainError> => ln_column;
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

    /// dimensionless_mass_from_column API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while allowing an explicit mass source.
    pub fn dimensionless_mass_from_column(
        &mut self,
        id: &str,
        source_col: &str,
        from: f64,
        to: f64,
        new_col: &str,
    ) -> Result<(), TGADomainError> {
        self.try_transform_by_id(id, |exp| {
            exp.dimensionless_mass_from_column(source_col, from, to, new_col)
        })
    }

    /// derive_conversion API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while delegating core logic
    /// to dataset and engine modules in experimental kinetics.
    pub fn derive_conversion(&mut self, id: &str) -> Result<(), TGADomainError> {
        self.try_transform_by_id(id, |exp| exp.derive_conversion())
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

    /// conversion_from_column API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while allowing an explicit mass source.
    pub fn conversion_from_column(
        &mut self,
        id: &str,
        source_col: &str,
        from: f64,
        to: f64,
        new_col: &str,
    ) -> Result<(), TGADomainError> {
        self.try_transform_by_id(id, |exp| {
            exp.conversion_from_column(source_col, from, to, new_col)
        })
    }

    pub fn conversion_with_m0(
        &mut self,
        id: &str,
        m0: f64,
        new_col: &str,
    ) -> Result<(), TGADomainError> {
        self.try_transform_by_id(id, |exp| exp.conversion_with_m0(m0, new_col))
    }
    //================================================================================================
    //      DIAGNOSTICS AND TESTING
    pub fn create_from_synthetic_data(
        &mut self,
        virtga: &VirtualTGA,
        meta: ExperimentMeta,
    ) -> Result<(), TGADomainError> {
        let exp = TGAExperiment::create_from_synthetic_data(virtga, meta)?;
        self.push(exp)?;
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

    /// derive_mass_rate_from_column API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while allowing an explicit mass source.
    pub fn derive_mass_rate_from_column(
        &mut self,
        id: &str,
        source_col: &str,
        new_col: &str,
    ) -> Result<(), TGADomainError> {
        self.try_transform_by_id(id, |exp| {
            exp.derive_mass_rate_from_column(source_col, new_col)
        })
    }

    /// derive_temperature_rate_from_column API method in the experiment/series facade.
    /// It keeps call-sites ergonomic while allowing an explicit temperature source.
    pub fn derive_temperature_rate_from_column(
        &mut self,
        id: &str,
        source_col: &str,
        new_col: &str,
    ) -> Result<(), TGADomainError> {
        self.try_transform_by_id(id, |exp| {
            exp.derive_temperature_rate_from_column(source_col, new_col)
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
    tgaseries_forward_try_apply! {
        pub fn sample_oneframeplot(&self, id: &str, range: Option<ViewRange>, max_points: usize) -> Result<PlotSeries, TGADomainError> => sample_oneframeplot;
        pub fn sample_column(&self, id: &str, col_name: &str, range: Option<(f64, f64)>, n_points: usize) -> Result<Vec<f64>, TGADomainError> => sample_column;
    }
    tgaseries_forward_apply! {
        pub fn list_of_columns(&self, id: &str) -> Result<Vec<String>, TGADomainError> => list_of_columns;
    }
    tgaseries_forward_try_apply! {
        pub fn plot_xy_ranges(&self, id: &str) -> Result<Ranges, TGADomainError> => plot_xy_ranges;
        pub fn get_temperature_col(&self, id: &str) -> Result<String, TGADomainError> => get_temperature_col;
        pub fn get_mass_col(&self, id: &str) -> Result<String, TGADomainError> => get_mass_col;
        pub fn get_time_col(&self, id: &str) -> Result<String, TGADomainError> => get_time_col;
        pub fn oneframeplot_axis_name(&self, id: &str, axis: XY) -> Result<String, TGADomainError> => oneframeplot_axis_name;
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

    tgaseries_forward_apply! {
        pub fn history_of_operations(&self, id: &str) -> Result<Vec<OperationRecord>, TGADomainError> => history_of_operations;
    }

    tgaseries_forward_apply! {
        pub fn operations_on_column(&self, id: &str, col: &str) -> Result<Vec<OperationRecord>, TGADomainError> => operations_on_column;
        pub fn get_column_history(&self, id: &str, col: &str) -> Result<ColumnHistory, TGADomainError> => get_column_history;
    }

    tgaseries_forward_try_apply! {
        pub fn mean_on_interval(&self, id: &str, value_col: &str, time_col: &str, from: f64, to: f64) -> Result<f64, TGADomainError> => mean_on_interval;
        pub fn mean_on_interval_on_own_range(&self, id: &str, col: &str, from: f64, to: f64) -> Result<f64, TGADomainError> => mean_on_interval_on_own_range;
        pub fn mean_on_column(&self, id: &str, col: &str) -> Result<f64, TGADomainError> => mean_on_column;
    }

    tgaseries_forward_mutate! {
        pub fn take_column(&mut self, id: &str, column_name: &str) -> Result<Option<String>, TGADomainError> => take_column;
        pub fn list_of_columns_to_recalc(&mut self, id: &str) -> Result<Vec<String>, TGADomainError> => list_of_columns_to_recalc;
    }

    tgaseries_forward_transform! {
        pub fn set_heating_rate(&mut self, id: &str, rate: f64) -> Result<(), TGADomainError> => with_heating_rate;
        pub fn set_comment(&mut self, id: &str, comment: &str) -> Result<(), TGADomainError> => with_comment;
        pub fn set_experiment_temperature(&mut self, id: &str, T: f64) -> Result<(), TGADomainError> => with_isothermal_temperature;
    }

    /// Revert the last dataset mutation for one experiment selected by id.
    pub fn undo_last(&mut self, id: &str) -> Result<(), TGADomainError> {
        self.try_transform_by_id(id, |exp| exp.undo_last())
    }

    tgaseries_forward_apply! {
        pub fn get_column_by_nature(&self, id: &str, nature: ColumnNature) -> Result<Option<String>, TGADomainError> => get_column_by_nature;
        pub fn get_columns_by_nature(&self, id: &str, nature: Vec<ColumnNature>) -> Result<Vec<Option<String>>, TGADomainError> => get_columns_by_nature;
    }
}
