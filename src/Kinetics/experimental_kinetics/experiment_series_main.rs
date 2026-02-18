use crate::Kinetics::experimental_kinetics::exp_engine_api::{PlotSeries, Ranges, ViewRange, XY};
use crate::Kinetics::experimental_kinetics::one_experiment_dataset::ColumnRole;
use crate::Kinetics::experimental_kinetics::one_experiment_dataset::{
    ColumnMeta, TGADataset, TGADomainError, UnaryOp, Unit,
};
use crate::Kinetics::experimental_kinetics::splines::SplineKind;
use std::collections::HashMap;
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
    pub fn conversion(self, src: &str, new_col: &str) -> Result<Self, TGADomainError> {
        let dataset = self.dataset.conversion(src, new_col)?;
        Ok(Self {
            dataset,
            meta: self.meta,
        })
    }

    //================================================================================================
    //      DIAGNOSTICS AND TESTING
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
    pub fn derive_rate(self, source_col: &str, new_col: &str) -> Result<Self, TGADomainError> {
        let dataset = self.dataset.derive_rate(source_col, new_col)?;
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
#[derive(Clone, Debug)]
/// TGASeries data structure for the experimental kinetics API layer.
/// It groups domain data and exposes a stable facade over
/// dataset-level processing modules used throughout KiThe.
pub struct TGASeries {
    pub experiments: Vec<TGAExperiment>,
    pub exp_map: HashMap<String, usize>,
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

    fn rebuild_index(&mut self) {
        self.exp_map.clear();
        for (idx, exp) in self.experiments.iter().enumerate() {
            self.exp_map.insert(exp.meta.id.clone(), idx);
        }
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
    pub fn conversion(&mut self, id: &str, src: &str, new_col: &str) -> Result<(), TGADomainError> {
        self.try_transform_by_id(id, |exp| exp.conversion(src, new_col))
    }

    //================================================================================================
    //      DIAGNOSTICS AND TESTING
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
    ) -> Result<(), TGADomainError> {
        self.try_transform_by_id(id, |exp| exp.derive_rate(source_col, new_col))
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

    pub fn sample_oneframeplot(
        &mut self,
        id: &str,
        range: Option<ViewRange>,
        max_points: usize,
    ) -> Result<PlotSeries, TGADomainError> {
        self.try_apply_by_id(id, |exp| exp.sample_oneframeplot(range, max_points))
    }

    pub fn list_of_columns(&self, id: &str) -> Result<Vec<String>, TGADomainError> {
        Ok(self.apply_by_id(id, |exp| exp.list_of_columns())?)
    }

    pub fn plot_xy_ranges(&self, id: &str) -> Result<Ranges, TGADomainError> {
        let r = self.apply_by_id(id, |exp| exp.plot_xy_ranges())?;
        r
    }
}
