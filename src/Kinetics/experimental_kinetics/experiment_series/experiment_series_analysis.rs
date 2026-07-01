//! Аналитический слой `TGAExperiment`.
//!
//! Здесь живут сглаживание, фильтры, ресемплинг, oneframeplot helpers,
//! статистика по окнам и вспомогательные wrappers для работы с осями.
//! Модуль оставлен плоским намеренно: здесь важнее читабельность API,
//! чем ещё одна прослойка абстракции.

use super::*;

impl TGAExperiment {
    //======================================================================
    // SMOOTHING AND FILTERING
    //======================================================================
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

    //======================================================================
    // GUI / ONEFRAMEPLOT HELPERS
    //======================================================================
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

    pub fn set_oneframeplot_y(self, y_col: &str) -> Result<Self, TGADomainError> {
        let dataset = self.dataset.set_oneframeplot_y(y_col)?;
        Ok(Self {
            dataset,
            meta: self.meta,
        })
    }

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

    pub fn cut_before_x_or_y(self, axis: XY, start_value: f64) -> Result<Self, TGADomainError> {
        let dataset = self.dataset.cut_before_x_or_y(axis, start_value)?;
        Ok(Self {
            dataset,
            meta: self.meta,
        })
    }

    pub fn cut_after_x_or_y(self, axis: XY, end_value: f64) -> Result<Self, TGADomainError> {
        let dataset = self.dataset.cut_after_x_or_y(axis, end_value)?;
        Ok(Self {
            dataset,
            meta: self.meta,
        })
    }

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

    pub fn min_distance_to_oneframeplot_point(
        &self,
        point: (f64, f64),
    ) -> Result<f64, TGADomainError> {
        self.dataset.min_distance_to_oneframeplot_point(point)
    }

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

    //======================================================================
    // STATS / COLUMN HELPERS
    //======================================================================
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

    /// Returns `true` when the experiment has at least one reversible snapshot.
    pub fn can_undo(&self) -> bool {
        self.dataset.can_undo()
    }

    /// Revert the last dataset mutation for this experiment.
    pub fn undo_last(mut self) -> Result<Self, TGADomainError> {
        self.dataset.undo_last()?;
        Ok(self)
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
        scale: f64,
        from: f64,
        to: f64,
    ) -> Result<Self, TGADomainError> {
        let dataset = self
            .dataset
            .scale_y_column_in_range_by_x_reference(scale, from, to)?;
        Ok(Self {
            dataset,
            meta: self.meta,
        })
    }

    pub fn scale_y_column_in_its_range(
        self,
        scale: f64,
        from: f64,
        to: f64,
    ) -> Result<Self, TGADomainError> {
        let dataset = self.dataset.scale_y_column_in_its_range(scale, from, to)?;
        Ok(Self {
            dataset,
            meta: self.meta,
        })
    }

    pub fn scale_x_column_in_its_range(
        self,
        scale: f64,
        from: f64,
        to: f64,
    ) -> Result<Self, TGADomainError> {
        let dataset = self.dataset.scale_x_column_in_its_range(scale, from, to)?;
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

    /// Human-readable provenance chain for one column.
    pub fn column_provenance_text(&self, col: &str) -> Option<String> {
        self.dataset.column_provenance_text(col)
    }

    pub fn validate_for_kinetics(&self) -> Result<(), TGADomainError> {
        self.dataset
            .validate_required_columns(&["temperature", "alpha", "dalpha_dt"])
    }
}
