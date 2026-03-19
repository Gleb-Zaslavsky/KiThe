use crate::Kinetics::experimental_kinetics::LSQSplines::SolverKind;
use crate::Kinetics::experimental_kinetics::exp_engine_api::GoldenPipelineConfig;
use crate::Kinetics::experimental_kinetics::exp_engine_api::{PlotSeries, Ranges, ViewRange, XY};
use crate::Kinetics::experimental_kinetics::exp_kinetics_smooth_filter::HampelStrategy;
use crate::Kinetics::experimental_kinetics::experiment_series_main::{
    ExperimentMeta, TGAExperiment, TGASeries,
};
use crate::Kinetics::experimental_kinetics::experiment_series2::SampledColumns;
use crate::Kinetics::experimental_kinetics::kinetic_methods::KineticDataView;
use crate::Kinetics::experimental_kinetics::kinetic_methods::integral_isoconversion::IsoconversionalResult;
use crate::Kinetics::experimental_kinetics::kinetic_methods::is_this_a_sublimation::{
    SublimationMethod, SublimationResult,
};
use crate::Kinetics::experimental_kinetics::kinetic_methods::isoconversion::IsoconversionalMethod;
use crate::Kinetics::experimental_kinetics::lowess_wrapper::LowessConfig;
use crate::Kinetics::experimental_kinetics::one_experiment_dataset::{
    ColumnHistory, ColumnNature, OperationRecord, TGADomainError, Unit,
};
use crate::Kinetics::experimental_kinetics::splines::SplineKind;
use crate::Kinetics::experimental_kinetics::testing_mod::{AdvancedTGAConfig, VirtualTGA};
use crate::gui::experimental_kinetics_gui::model::{Colours, PlotCurve, PlotModel, TGAGUIError};
use log::{debug, info};
use std::collections::HashMap;
use std::path::Path;

//=========================================================================================
#[derive(Debug, Clone)]
pub struct SampledColumn {
    pub col: Vec<f64>,
    pub name: String,
}
//===========================================================================================
/// Builder для пошагового создания PlotCurve / Builder for step-by-step PlotCurve construction

pub struct PlotCurveBuilder {
    curve: PlotCurve,
}

impl PlotCurveBuilder {
    /// Создает builder из сырых данных графика / Creates builder from raw plot data
    pub fn new_from_plot(
        plot: PlotSeries,
        id: &str,
        index: usize,
        ranges: Ranges,
        short_name: String,
    ) -> Self {
        let x_name = plot.name_x;
        let y_name = plot.name_y;
        let points: Vec<[f64; 2]> = plot
            .x
            .into_iter()
            .zip(plot.y)
            .map(|(xi, yi)| [xi, yi])
            .collect();

        let mut curve = PlotCurve::default();
        curve.experiment_id = id.to_string();
        curve.experiment_index = index;
        curve.plot_short_name = short_name;
        curve.x_name = x_name;
        curve.y_name = y_name;
        curve.ranges = ranges;
        curve.points = points;

        Self { curve }
    }

    /// Устанавливает цвет через enum / Sets colour via enum
    pub fn with_colour(mut self, colour: Colours) -> Self {
        self.curve.set_colour(colour);
        self
    }

    /// Устанавливает цвет напрямую RGB / Sets colour directly as RGB
    pub fn with_colour_rgb(mut self, rgb: [u8; 3]) -> Self {
        self.curve.set_colour_rgb(rgb);
        self
    }

    /// Устанавливает флаг видимости / Sets visibility flag
    pub fn with_shown(mut self, shown: bool) -> Self {
        self.curve.shown = shown;
        self
    }

    /// Устанавливает флаг выбора / Sets selection flag
    pub fn with_selected(mut self, selected: bool) -> Self {
        self.curve.selected = selected;
        self
    }

    /// Устанавливает флаг подсветки / Sets highlight flag
    pub fn with_highlighted(mut self, highlighted: bool) -> Self {
        self.curve.highlighted = highlighted;
        self
    }

    /// Завершает сборку и возвращает кривую / Finalizes and returns the curve
    pub fn finish(self) -> PlotCurve {
        self.curve
    }
}
///===============================================================================================
impl PlotModel {
    pub fn push_message(&mut self, message: &str) -> Result<(), TGAGUIError> {
        self.message = message.to_string();
        Ok(())
    }

    pub fn clear_message(&mut self) -> Result<(), TGAGUIError> {
        self.message = String::new();
        Ok(())
    }
    //////////////////////////////////////////////////////////////////////////////////
    //============================================================================
    // INPUT/OUTPUT
    //============================================================================
    //////////////////////////////////////////////////////////////////////////////////
    pub fn push_from_file(&mut self, path: &Path) -> Result<(), TGAGUIError> {
        self.series.push_from_file(path)?;
        Ok(())
    }
    pub fn to_csv_series(&self, path: &Path) -> Result<(), TGAGUIError> {
        self.series.to_csv_series(path)?;
        Ok(())
    }

    pub fn from_csv_series(&mut self, path: &Path) -> Result<(), TGAGUIError> {
        self.series = TGASeries::from_csv_series(path)?;
        Ok(())
    }
    //////////////////////////////////////////////////////////////////////////////////
    //=================================================================
    // SETTERS
    //=======================================================================
    //////////////////////////////////////////////////////////////////////////////////
    pub fn drop_experiment(&mut self, indx: usize) -> Result<(), TGAGUIError> {
        self.series.drop_experiment(indx);
        Ok(())
    }

    pub fn drop_experiment_by_id(&mut self, id: &str) -> Result<(), TGAGUIError> {
        self.series.drop_experiment_by_id(id)?;
        Ok(())
    }
    /// list of experiments
    pub fn list_of_experiments(&self) -> Vec<String> {
        self.series.ids()
    }

    pub fn get_list_of_complex_id(&self) -> Vec<String> {
        self.series.get_list_of_complex_id()
    }

    pub fn get_experiment_by_id_mut(
        &mut self,
        id: &str,
    ) -> Result<&mut TGAExperiment, TGAGUIError> {
        let exp = self.series.get_experiment_by_id_mut(id)?;
        Ok(exp)
    }

    pub fn bind_time_of_experiment(
        &mut self,
        id: &str,
        col: &str,
        unit: Unit,
    ) -> Result<(), TGAGUIError> {
        self.series.bind_time(id, col, unit)?;
        Ok(())
    }

    pub fn bind_temperature_of_experiment(
        &mut self,
        id: &str,
        col: &str,
        unit: Unit,
    ) -> Result<(), TGAGUIError> {
        self.series.bind_temperature(id, col, unit)?;
        Ok(())
    }

    pub fn bind_mass_of_experiment(
        &mut self,
        id: &str,
        col: &str,
        unit: Unit,
    ) -> Result<(), TGAGUIError> {
        self.series.bind_mass(id, col, unit)?;
        Ok(())
    }

    pub fn bind_time_of_selected(&mut self, col: &str, unit: Unit) -> Result<(), TGAGUIError> {
        let id = self.get_experiment_by_selected_curve()?;
        debug!("selected experiment {}", id);
        self.bind_time_of_experiment(&id, col, unit)?;
        Ok(())
    }

    pub fn bind_temperature_of_selected(
        &mut self,
        col: &str,
        unit: Unit,
    ) -> Result<(), TGAGUIError> {
        let id = self.get_experiment_by_selected_curve()?;
        debug!("selected experiment {}", id);
        self.bind_temperature_of_experiment(&id, col, unit)?;
        Ok(())
    }

    pub fn bind_mass_of_selected(&mut self, col: &str, unit: Unit) -> Result<(), TGAGUIError> {
        let id = self.get_experiment_by_selected_curve()?;
        debug!("selected experiment {}", id);
        self.bind_mass_of_experiment(&id, col, unit)?;
        Ok(())
    }

    pub fn create_from_synthetic_data(
        &mut self,
        virtga: &VirtualTGA,
        meta: ExperimentMeta,
    ) -> Result<(), TGAGUIError> {
        self.series.create_from_synthetic_data(virtga, meta)?;
        Ok(())
    }
    pub fn generate_synthetic_data_from_config(
        &mut self,
        config: &AdvancedTGAConfig,
    ) -> Result<(), TGAGUIError> {
        let simulated = VirtualTGA::generate_series(config);
        for exp in simulated {
            let meta = exp.meta;
            let vtga = exp.virtual_tga;
            self.series.create_from_synthetic_data(&vtga, meta)?;
        }
        Ok(())
    }

    pub fn compare_new_and_old_list_of_columns(
        &self,
        id: &str,
        old_list: Vec<String>,
    ) -> Result<Vec<String>, TGAGUIError> {
        info!("old list of columns: {:?}", old_list);
        let new_list = self.list_of_columns(id)?;
        let mut diff_list = Vec::new();

        for col in new_list {
            if !old_list.contains(&col) {
                diff_list.push(col);
            }
        }
        Ok(diff_list)
    }

    pub fn create_plots_with_columns(
        &mut self,
        id: &str,
        new_list: Vec<String>,
        preferred_x: Option<&str>,
    ) -> Result<(), TGAGUIError> {
        // Default: keep current x of selected curve.
        // But when a new time-like column is produced (e.g. splining),
        // force it as x to avoid plotting new columns against old x.
        let mut x_col = self.this_is_x()?;
        if let Some(px) = preferred_x {
            if new_list.iter().any(|c| c == px) {
                x_col = px.to_string();
            }
        }
        self.set_x(id, &x_col)?;
        for col in new_list {
            if col == x_col {
                continue;
            }
            self.set_y(id, &col)?;
            self.create_points_for_curve_with_builder(id, Colours::Red)?;
            self.plots.last_mut().unwrap().shown = true;
            info!("created plot for column {}", col);
        }
        Ok(())
    }

    pub fn create_plots_with_new_columns(
        &mut self,
        id: &str,
        old_list: Vec<String>,
        preferred_x: Option<&str>,
    ) -> Result<(), TGAGUIError> {
        let new_list = self.compare_new_and_old_list_of_columns(id, old_list)?;
        info!("new columns: {:?}", new_list);
        self.create_plots_with_columns(id, new_list, preferred_x)?;
        self.reset_view();
        Ok(())
    }

    pub fn create_experiment_from_columns(
        &mut self,
        parent_idx: usize,
        new_id: String,
        columns: &[&str],
    ) -> Result<(), TGAGUIError> {
        self.series
            .create_experiment_from_columns(parent_idx, new_id, columns)?;
        Ok(())
    }
    pub fn create_experiment_from_columns_for_experiment(
        &mut self,
        id: &str,
        new_id: String,
        columns: &[&str],
    ) -> Result<(), TGAGUIError> {
        let idx = self.index_by_id(id)?;
        self.create_experiment_from_columns(idx, new_id, columns)
    }

    //========================================================================
    //////////////////////////////////////////////////////////////////////////////////
    //========================================================================
    // GENERIC METHODS
    //=================================================================================
    //////////////////////////////////////////////////////////////////////////////////
    pub fn apply_by_id<R, F>(&self, id: &str, op: F) -> Result<R, TGAGUIError>
    where
        F: FnOnce(&TGAExperiment) -> R,
    {
        let res = self.series.apply_by_id(id, op)?;
        Ok(res)
    }

    /// try_apply_by_id API method in the experiment/series facade.
    /// Calls a fallible read-only operation on one experiment selected by id.
    /// generic argument:  FnOnce(&TGAExperiment) -> Result<R, TGADomainError>,
    pub fn try_apply_by_id<R, F>(&self, id: &str, op: F) -> Result<R, TGAGUIError>
    where
        F: FnOnce(&TGAExperiment) -> Result<R, TGADomainError>,
    {
        let res = self.series.try_apply_by_id(id, op)?;
        Ok(res)
    }

    /// mutate_by_id API method in the experiment/series facade.
    /// Calls an infallible mutable operation on one experiment selected by id.
    ///  generic argument:  FnOnce(&mut TGAExperiment) -> R,
    pub fn mutate_by_id<R, F>(&mut self, id: &str, op: F) -> Result<R, TGAGUIError>
    where
        F: FnOnce(&mut TGAExperiment) -> R,
    {
        let res = self.series.mutate_by_id(id, op)?;
        Ok(res)
    }

    /// try_mutate_by_id API method in the experiment/series facade.
    /// Calls a fallible mutable operation on one experiment selected by id.
    ///generic argument:  FnOnce(&mut TGAExperiment) -> Result<R, TGADomainError>,
    pub fn try_mutate_by_id<R, F>(&mut self, id: &str, op: F) -> Result<R, TGAGUIError>
    where
        F: FnOnce(&mut TGAExperiment) -> Result<R, TGADomainError>,
    {
        let res = self.series.try_mutate_by_id(id, op)?;
        Ok(res)
    }

    /// transform_by_id API method in the experiment/series facade.
    /// Applies an infallible consuming transformation (`TGAExperiment -> TGAExperiment`)
    /// to one experiment selected by id.
    /// generic argument: FnOnce(TGAExperiment) -> TGAExperiment,
    pub fn transform_by_id<F>(&mut self, id: &str, op: F) -> Result<(), TGAGUIError>
    where
        F: FnOnce(TGAExperiment) -> TGAExperiment,
    {
        self.series.transform_by_id(id, op)?;
        Ok(())
    }

    /// try_transform_by_id API method in the experiment/series facade.
    /// Applies a fallible consuming transformation
    /// (`TGAExperiment -> Result<TGAExperiment, TGADomainError>`)
    /// to one experiment selected by id.
    /// generic argument: FnOnce(TGAExperiment) -> Result<TGAExperiment, TGADomainError>,
    pub fn try_transform_by_id<F>(&mut self, id: &str, op: F) -> Result<(), TGAGUIError>
    where
        F: FnOnce(TGAExperiment) -> Result<TGAExperiment, TGADomainError>,
    {
        self.series.try_transform_by_id(id, op)?;

        Ok(())
    }
    //////////////////////////////////////////////////////////////////////////////////
    //=========================================================================
    //              CREATING TABLE
    //========================================================================
    //////////////////////////////////////////////////////////////////////////////////
    pub fn history_of_operation_for_selected(&self) -> Result<Vec<OperationRecord>, TGAGUIError> {
        let id = self.get_experiment_by_selected_curve()?;
        let h = self.series.history_of_operations(&id)?;
        Ok(h)
    }

    pub fn operations_on_column_for_experiment(
        &mut self,
        id: &str,
        col: &str,
    ) -> Result<Vec<OperationRecord>, TGAGUIError> {
        let r = self.series.operations_on_column(id, col)?;
        Ok(r)
    }
    /// Returns the history of operations on a column for an experiment.
    /// The `id` parameter should be the identifier of the experiment,
    /// and the `col` parameter should be the name of the column.
    pub fn get_column_history_for_experiment(
        &mut self,
        id: &str,
        col: &str,
    ) -> Result<ColumnHistory, TGAGUIError> {
        let r = self.series.get_column_history(id, col)?;
        Ok(r)
    }
    /// takes from vector of changed column names Vec<Option<String>> a name leaving None
    pub fn take_column(&mut self, id: &str, col: &str) -> Result<Option<String>, TGAGUIError> {
        let name = self.series.take_column(id, col)?;
        Ok(name)
    }

    pub fn list_of_columns_to_recalc(&mut self, id: &str) -> Result<Vec<String>, TGAGUIError> {
        let list_of_columns_to_recalc = self.series.list_of_columns_to_recalc(id)?;
        Ok(list_of_columns_to_recalc)
    }

    pub fn list_of_columns_to_recalc_for_all_experiments(
        &mut self,
    ) -> Result<HashMap<String, Vec<String>>, TGAGUIError> {
        let mut map_of_columns_to_recalc = HashMap::new();
        for id in self.list_of_experiments() {
            let list_of_columns_to_recalc = self.list_of_columns_to_recalc(&id)?;
            map_of_columns_to_recalc.insert(id, list_of_columns_to_recalc);
        }

        Ok(map_of_columns_to_recalc)
    }

    pub fn sample_column(
        &self,
        id: &str,
        col_name: &str,
        range: Option<(f64, f64)>,
        n_points: usize,
    ) -> Result<Vec<f64>, TGAGUIError> {
        let v = self.series.sample_column(id, col_name, range, n_points)?;
        Ok(v)
    }
    /// returns HashMap<experiment_id, HashMap<column_name, SampledColumn>>
    pub fn column_samples_for_all_experiments(
        &mut self,
        n_points: usize,
    ) -> Result<HashMap<String, HashMap<String, SampledColumn>>, TGAGUIError> {
        let mut data_map = HashMap::new();
        for id in self.list_of_experiments() {
            let list_of_columns_to_recalc = &self.list_of_columns_to_recalc(&id)?;
            info!(
                "this columns will be resampled {:?}",
                list_of_columns_to_recalc
            );
            let mut exp_columns = HashMap::new();
            for col in list_of_columns_to_recalc {
                let sampled_col = self.sample_column(&id, &col, None, n_points)?;
                info!(
                    "for column {} created sample of length {}",
                    col,
                    sampled_col.len()
                );
                let sampled_col_struct = SampledColumn {
                    name: col.clone(),
                    col: sampled_col,
                };
                exp_columns.insert(col.clone(), sampled_col_struct);
            }
            data_map.insert(id, exp_columns);
        }
        Ok(data_map)
    }

    pub fn column_samples_for_all_experiment_for_plotting(
        &mut self,
        n_points: usize,
    ) -> Result<(HashMap<String, Vec<f64>>, SampledColumns), TGAGUIError> {
        let res = self
            .series
            .column_samples_for_all_experiment_for_plotting(n_points)?;
        Ok(res)
    }

    pub fn set_heating_rate(&mut self, id: &str, rate: f64) -> Result<(), TGAGUIError> {
        self.series.set_heating_rate(id, rate)?;
        Ok(())
    }
    pub fn set_comment(&mut self, id: &str, comment: &str) -> Result<(), TGAGUIError> {
        self.series.set_comment(id, comment)?;
        Ok(())
    }

    pub fn set_experiment_temperature(&mut self, id: &str, T: f64) -> Result<(), TGAGUIError> {
        self.series.set_experiment_temperature(id, T)?;
        Ok(())
    }

    pub fn monotony_of_time_check_for_experiment(&self, id: &str) -> Result<Vec<f64>, TGAGUIError> {
        Ok(self.series.monotony_of_time_check(id)?)
    }
    //////////////////////////////////////////////////////////////////////////////////
    //=================================================================================
    //=================================================================================
    // CREATING PLOT
    //=========================================================================
    //=================================================================================
    //
    pub fn set_x(&mut self, id: &str, x_col: &str) -> Result<(), TGAGUIError> {
        self.series.set_oneframeplot_x(id, x_col)?;
        Ok(())
    }
    pub fn set_y(&mut self, id: &str, y_col: &str) -> Result<(), TGAGUIError> {
        self.series.set_oneframeplot_y(id, y_col)?;
        Ok(())
    }
    pub fn list_of_columns(&self, id: &str) -> Result<Vec<String>, TGAGUIError> {
        let list = self.series.list_of_columns(id)?;
        Ok(list)
    }
    pub fn list_of_columns_for_selected(&self) -> Result<Vec<String>, TGAGUIError> {
        let id = self.get_experiment_by_selected_curve()?;
        self.list_of_columns(&id)
    }
    /// max and min of x and y
    pub fn plot_range(&mut self, id: &str) -> Result<Ranges, TGAGUIError> {
        debug!("setting plot named {} the range", id);
        let ranges = self.series.plot_xy_ranges(id)?;
        Ok(ranges)
    }
    // initial ranges - now combines with existing plots
    pub fn give_ranges_to_plot(&mut self, id: &str) -> Result<Ranges, TGAGUIError> {
        debug!("giving ranges to plot named {}", id);

        // Range of the curve currently being (re)built.
        let ranges = self.plot_range(id)?;

        // Build this curve short-name so we can ignore stale previous version of it
        // while computing global view bounds.
        let x_name = self.series.oneframeplot_axis_name(id, XY::X)?;
        let y_name = self.series.oneframeplot_axis_name(id, XY::Y)?;
        let short_name = self.create_plot_short_name(id, &x_name, &y_name);

        let mut x_min = ranges.x_min;
        let mut x_max = ranges.x_max;
        let mut y_min = ranges.y_min;
        let mut y_max = ranges.y_max;

        for curve in self.plots.iter().filter(|curve| curve.is_shown()) {
            if curve.plot_short_name == short_name {
                continue;
            }
            x_min = x_min.min(curve.ranges.x_min);
            x_max = x_max.max(curve.ranges.x_max);
            y_min = y_min.min(curve.ranges.y_min);
            y_max = y_max.max(curve.ranges.y_max);
        }

        self.interaction.view_range = (x_min, x_max);
        self.interaction.view_y_range = (y_min, y_max);

        Ok(ranges)
    }
    pub fn calculate_samplepoints_for_plot(
        &mut self,
        id: &str,
        x_range: [f64; 2],
    ) -> Result<PlotSeries, TGAGUIError> {
        let range = Some({
            ViewRange {
                t_min: x_range[0],
                t_max: x_range[1],
            }
        });
        let plot = self
            .series
            .sample_oneframeplot(id, range, self.settings.n_points().unwrap())?;
        Ok(plot)
    }
    /// get number of plot by itd String id
    pub fn index_by_id(&self, id: &str) -> Result<usize, TGAGUIError> {
        let res = self.series.index_by_id(id)?;
        Ok(res)
    }
    pub fn create_points_for_curve(&mut self, id: &str) -> Result<(), TGAGUIError> {
        debug!("calculating points for curve");
        let ranges = self.give_ranges_to_plot(id)?;
        // Always sample full current curve extent after transformations (cut, scale, etc.),
        // not just the current viewport.
        let x_range = [ranges.x_min, ranges.x_max];
        let plot = self.calculate_samplepoints_for_plot(id, x_range)?;
        let index = self.index_by_id(id)?;
        self.push_to_vec_of_curves(plot, id, index, ranges)?;
        Ok(())
    }
    pub fn create_points_for_selected_curve(&mut self) -> Result<(), TGAGUIError> {
        let id = self.get_experiment_by_selected_curve()?;
        self.create_points_for_curve(&id)?;
        Ok(())
    }

    /// Альтернативный путь: создание кривой через builder с выбором цвета / Alternative path: create curve via builder with colour selection
    pub fn create_points_for_curve_with_builder(
        &mut self,
        id: &str,
        colour: Colours,
    ) -> Result<(), TGAGUIError> {
        debug!("calculating points for curve via builder");
        let ranges = self.give_ranges_to_plot(id)?;
        // Всегда берём полный диапазон кривой / Always use full curve range.
        let x_range = [ranges.x_min, ranges.x_max];
        let plot = self.calculate_samplepoints_for_plot(id, x_range)?;
        let index = self.index_by_id(id)?;
        self.push_to_vec_of_curves_with_builder(plot, id, index, ranges, colour)?;
        Ok(())
    }

    /// Альтернативный путь для выбранной кривой / Alternative builder path for selected curve
    pub fn create_points_for_selected_curve_with_builder(
        &mut self,
        colour: Colours,
    ) -> Result<(), TGAGUIError> {
        let id = self.get_experiment_by_selected_curve()?;
        self.create_points_for_curve_with_builder(&id, colour)?;
        Ok(())
    }

    pub fn create_plot_short_name(&self, id: &str, x_name: &String, y_name: &String) -> String {
        format!("{}_{}_{}", id.to_string(), x_name, y_name)
    }

    pub fn push_to_vec_of_curves(
        &mut self,
        plot: PlotSeries,
        id: &str,
        index: usize,
        ranges: Ranges,
    ) -> Result<(), TGAGUIError> {
        let mut new_curve = PlotCurve::default();
        new_curve.experiment_id = id.to_string().clone();
        new_curve.experiment_index = index;
        let x_name = plot.name_x;
        let y_name = plot.name_y;
        let short_name = self.create_plot_short_name(id, &x_name, &y_name);
        new_curve.x_name = x_name;
        new_curve.y_name = y_name;
        let x_vec = plot.x;
        let y_vec = plot.y;
        let mut points = Vec::new();
        for (xi, yi) in x_vec.iter().zip(y_vec) {
            points.push([*xi, yi]);
        }
        new_curve.ranges = ranges;
        new_curve.points = points;
        new_curve.plot_short_name = short_name;
        // Replace an existing curve for the same experiment/x/y instead of duplicating
        // stale pre-recalculation data.
        if let Some(existing) = self
            .plots
            .iter()
            .position(|curve| curve.plot_short_name == new_curve.plot_short_name)
        {
            let prev = &self.plots[existing];
            new_curve.color = prev.color;
            new_curve.shown = prev.shown;
            new_curve.selected = prev.selected;
            new_curve.highlighted = prev.highlighted;
            self.plots[existing] = new_curve;
        } else {
            self.plots.push(new_curve);
        }
        Ok(())
    }

    /// Альтернативный push через builder / Alternative push path via builder
    pub fn push_to_vec_of_curves_with_builder(
        &mut self,
        plot: PlotSeries,
        id: &str,
        index: usize,
        ranges: Ranges,
        colour: Colours,
    ) -> Result<(), TGAGUIError> {
        let short_name = self.create_plot_short_name(id, &plot.name_x, &plot.name_y);
        let base_builder = PlotCurveBuilder::new_from_plot(plot, id, index, ranges, short_name);
        let configured_builder = base_builder.with_colour(colour);
        let mut new_curve = configured_builder.finish();

        // При замене сохраняем UI-состояние старой кривой / Preserve old UI state on replacement.
        if let Some(existing) = self
            .plots
            .iter()
            .position(|curve| curve.plot_short_name == new_curve.plot_short_name)
        {
            let prev = &self.plots[existing];
            new_curve = PlotCurveBuilder { curve: new_curve }
                .with_shown(prev.shown)
                .with_selected(prev.selected)
                .with_highlighted(prev.highlighted)
                .finish();
            self.plots[existing] = new_curve;
        } else {
            self.plots.push(new_curve);
        }
        Ok(())
    }
    //////////////////////////////////////////////////////////////////////////////////
    //=================================================================================
    //COLUMNS TRASFORMATION
    //================================================================================
    //////////////////////////////////////////////////////////////////////////////////
    pub fn rename_column_for_experiment(
        &mut self,
        id: &str,
        col_name: &str,
        new_name: &str,
    ) -> Result<(), TGAGUIError> {
        self.series.rename_column(id, col_name, new_name)?;
        self.reset_view();
        Ok(())
    }
    pub fn drop_column_for_experiment(
        &mut self,
        id: &str,
        col_name: &str,
    ) -> Result<(), TGAGUIError> {
        self.series.drop_column(id, col_name)?;
        self.reset_view();
        Ok(())
    }

    pub fn drop_column_for_selected(&mut self, col_name: &str) -> Result<(), TGAGUIError> {
        let id = self.get_experiment_by_selected_curve()?;
        self.drop_column_for_experiment(&id, col_name)?;

        Ok(())
    }

    pub fn move_time_to_zero_for_experiment(&mut self, id: &str) -> Result<(), TGAGUIError> {
        self.series.move_time_to_zero(id)?;
        self.create_points_for_curve(&id)?;
        self.reset_view();
        Ok(())
    }

    pub fn move_time_to_zero_for_selected(&mut self) -> Result<(), TGAGUIError> {
        let id = self.get_experiment_by_selected_curve()?;
        debug!("selected experiment {}", id);
        info!("moving time to zero!");
        self.move_time_to_zero_for_experiment(&id)?;
        Ok(())
    }

    // UNITS TRANSFORMATIONS
    //=====================================
    #[allow(non_snake_case)]
    pub fn from_C_to_K_for_experiment(&mut self, id: &str) -> Result<(), TGAGUIError> {
        self.series.celsius_to_kelvin(id)?;
        self.create_points_for_curve(&id)?;
        self.reset_view();
        Ok(())
    }
    #[allow(non_snake_case)]
    pub fn from_C_to_K_of_selected(&mut self) -> Result<(), TGAGUIError> {
        let id = self.get_experiment_by_selected_curve()?;
        debug!("selected experiment {}", id);
        info!("T converted from C to K");
        self.from_C_to_K_for_experiment(&id)?;

        Ok(())
    }

    pub fn from_s_to_h_for_experiment(&mut self, id: &str) -> Result<(), TGAGUIError> {
        self.series.seconds_to_hours(id)?;
        self.create_points_for_curve(&id)?;
        self.reset_view();
        Ok(())
    }

    pub fn from_s_to_h_of_selected(&mut self) -> Result<(), TGAGUIError> {
        let id = self.get_experiment_by_selected_curve()?;
        debug!("selected experiment {}", id);
        info!("time turned from seconds to hours");
        self.from_s_to_h_for_experiment(&id)?;
        Ok(())
    }
    // CUTTING AND DELETING
    //===============================================================================
    pub fn cut_before_time_for_selected(&mut self, time: f64) -> Result<(), TGAGUIError> {
        let id = self.get_experiment_by_selected_curve()?;
        debug!("selected experiment {}", id);
        info!("cut all data before time {}", time);
        self.series.cut_before_time(&id, time)?;
        self.create_points_for_curve(&id)?;
        self.reset_view();
        Ok(())
    }
    /// эта функция должна удалять выбранный график
    /// внимание она не дропает колонку датафрейма - только график!
    pub fn delete_selected_curve(&mut self) -> Result<(), TGAGUIError> {
        let maybe_index = self.get_selected_curve_index().or_else(|| {
            // Fallback: when only one curve is shown, allow deleting it without explicit selection.
            if self.is_only_one_shown() {
                self.plots.iter().position(|curve| curve.is_shown())
            } else {
                None
            }
        });

        if let Some(index) = maybe_index {
            self.plots.remove(index);
            self.clear_selection_rect();
            self.reset_view();
            Ok(())
        } else {
            Err(TGAGUIError::BindingError(
                "No selected curve to delete".to_string(),
            ))
        }
    }

    pub fn delete_all_plots(&mut self) -> Result<(), TGAGUIError> {
        self.plots.clear();
        Ok(())
    }
    pub fn cut_range_x_or_y_for_experiment(
        &mut self,
        id: &str,
        axis: XY,
    ) -> Result<(), TGAGUIError> {
        if let Some(range) = self.interaction.selection_rect {
            let bounds = range.bounds();
            let x_min = bounds.0;
            let x_max = bounds.1;
            let y_min = bounds.2;
            let y_max = bounds.3;
            match axis {
                XY::X => {
                    info!("cutting range {:?} on X axe", (x_min, x_max));
                    self.series
                        .cut_range_inverse_x_or_y(id, XY::X, x_min, x_max)?;
                }
                XY::Y => {
                    info!("cutting range {:?} on Y axe", (y_min, y_max));
                    self.series
                        .cut_range_inverse_x_or_y(id, XY::Y, y_min, y_max)?;
                }
            }
            self.create_points_for_curve(&id)?;
            self.reset_view();
        } else {
            return Err(TGAGUIError::BindingError(
                "No selection rectangle. Drag-select a region on the plot first.".to_string(),
            ));
        }

        Ok(())
    }

    pub fn cut_range_x_or_y_for_selected(&mut self, axis: XY) -> Result<(), TGAGUIError> {
        let id = self.get_experiment_by_selected_curve()?;
        self.cut_range_x_or_y_for_experiment(&id, axis)
    }

    // RELATIVE MASS AND CONVERSION
    //=============================================================
    pub fn dimensionless_mass_for_experiment(
        &mut self,
        id: String,
        t_end: f64,
        new_col: &str,
    ) -> Result<(), TGAGUIError> {
        self.series.dimensionless_mass(&id, 0.0, t_end, new_col)?;
        Ok(())
    }

    pub fn relative_mass_for_selected(
        &mut self,
        t_end: f64,
        new_col: Option<&str>,
    ) -> Result<(), TGAGUIError> {
        let id = self.get_experiment_by_selected_curve()?;
        match new_col {
            Some(new_col) => {
                self.dimensionless_mass_for_experiment(id, t_end, new_col)?;
            }
            None => {
                self.dimensionless_mass_for_experiment(id, t_end, "relative_mass")?;
            }
        }
        Ok(())
    }

    pub fn conversion_for_experiment(
        &mut self,
        id: String,
        t_end: f64,
        new_col: &str,
    ) -> Result<(), TGAGUIError> {
        self.series.conversion(&id, 0.0, t_end, new_col)?;
        Ok(())
    }

    pub fn conversion_for_selected(
        &mut self,
        t_end: f64,
        new_col: Option<&str>,
    ) -> Result<(), TGAGUIError> {
        let id = self.get_experiment_by_selected_curve()?;
        match new_col {
            Some(new_col) => {
                self.conversion_for_experiment(id, t_end, new_col)?;
            }
            None => {
                self.conversion_for_experiment(id, t_end, "conversion")?;
            }
        }
        Ok(())
    }
    // CALIBRATION: MASS FROM VOLTAGE
    //===========================================================================
    pub fn calibrate_mass_from_voltage(&mut self, id: &str) -> Result<(), TGAGUIError> {
        let set = self.settings.calibration_line().unwrap();
        let k = set.k();
        let b = set.b();
        info!(
            "calibrating mass from voltage with coefficients: k= {}, b= {} ",
            k, b
        );
        self.series.calibrate_mass_from_voltage(id, k, b)?;
        self.create_points_for_curve(&id)?;
        self.reset_view();
        Ok(())
    }

    pub fn calibrate_mass_from_voltage_for_selected(&mut self) -> Result<(), TGAGUIError> {
        let id = self.get_experiment_by_selected_curve()?;
        debug!("selected experiment {}", id);
        self.calibrate_mass_from_voltage(&id)?;
        Ok(())
    }

    pub fn calibrate_mass_from_voltage_new_column(
        &mut self,
        id: &str,
        col: &str,
    ) -> Result<(), TGAGUIError> {
        let set = self.settings.calibration_line().unwrap();
        let k = set.k();
        let b = set.b();
        self.series.calibrate_mass(id, k, b, col)?;
        self.create_points_for_curve(&id)?;
        self.reset_view();
        Ok(())
    }

    pub fn calibrate_mass_from_voltage_new_colum_for_selected(
        &mut self,
        col: &str,
    ) -> Result<(), TGAGUIError> {
        let id = self.get_experiment_by_selected_curve()?;
        debug!("selected experiment {}", id);
        self.calibrate_mass_from_voltage_new_column(&id, col)?;
        Ok(())
    }

    pub fn calibrate_mass_from_voltage_with_new_optional_column_for_selected(
        &mut self,
        col: Option<&str>,
    ) -> Result<(), TGAGUIError> {
        match col {
            Some(col) => self.calibrate_mass_from_voltage_new_colum_for_selected(col),
            None => self.calibrate_mass_from_voltage_for_selected(),
        }
    }

    //===================================================================================
    //  X or Y column operations
    //===================================================================================
    pub fn this_is_x_or_y(&self, xy: XY) -> Result<String, TGAGUIError> {
        let id = self.get_experiment_by_selected_curve()?;
        let x_name = self.series.oneframeplot_axis_name(&id, xy)?;
        Ok(x_name)
    }
    pub fn this_is_x(&self) -> Result<String, TGAGUIError> {
        let id = self.get_experiment_by_selected_curve()?;
        let x_name = self.series.oneframeplot_axis_name(&id, XY::X)?;
        Ok(x_name)
    }
    pub fn this_is_y(&self) -> Result<String, TGAGUIError> {
        let id = self.get_experiment_by_selected_curve()?;
        let x_name = self.series.oneframeplot_axis_name(&id, XY::Y)?;
        Ok(x_name)
    }

    // + AND -
    //=====================================================================
    // OFFSET ALL COLUMN
    pub fn offset_column(
        &mut self,
        id: &str,
        column: &str,
        offset: f64,
    ) -> Result<(), TGAGUIError> {
        self.series.offset_column(id, column, offset)?;

        Ok(())
    }

    pub fn offset_column_for_selected(
        &mut self,
        column: &str,
        offset: f64,
    ) -> Result<(), TGAGUIError> {
        let id = self.get_experiment_by_selected_curve()?;
        debug!("for selected experiment {}", id);
        self.series.offset_column(&id, column, offset)?;
        self.create_points_for_curve(&id)?;
        self.reset_view();
        Ok(())
    }

    pub fn add_column_for_selected(
        &mut self,
        column: &str,
        offset: f64,
    ) -> Result<(), TGAGUIError> {
        info!("add for column {} by {}", column, offset);
        self.offset_column_for_selected(column, offset)
    }
    pub fn sub_column_for_selected(
        &mut self,
        column: &str,
        offset: f64,
    ) -> Result<(), TGAGUIError> {
        info!("subtraction for column {} by {}", column, offset);
        self.offset_column_for_selected(column, -offset)
    }
    // OFFSET ONLY SELECTED
    pub fn offset_column_in_its_range_for_experiment(
        &mut self,
        id: &str,
        colmn: &str,
        offset: f64,
        from: f64,
        to: f64,
    ) -> Result<(), TGAGUIError> {
        self.series
            .offset_column_in_its_range(id, colmn, offset, from, to)?;

        Ok(())
    }
    pub fn offset_column_in_its_range_for_selected(
        &mut self,
        colmn: &str,
        offset: f64,
    ) -> Result<(), TGAGUIError> {
        let id = self.get_experiment_by_selected_curve()?;
        if let Some(range) = self.interaction.selection_rect {
            let bounds = range.bounds();

            let x_min = bounds.0;
            let x_max = bounds.1;
            let x_col = self.series.oneframeplot_axis_name(&id, XY::X)?;
            self.series
                .offset_column_in_range_by_reference(&id, colmn, &x_col, offset, x_min, x_max)?;
            self.create_points_for_curve(&id)?;
            self.reset_view();
        } else {
            return Err(TGAGUIError::BindingError(
                "No selection rectangle. Drag-select a region on the plot first.".to_string(),
            ));
        }

        Ok(())
    }
    pub fn add_column_in_its_range_for_selected(
        &mut self,
        colmn: &str,
        offset: f64,
    ) -> Result<(), TGAGUIError> {
        self.offset_column_in_its_range_for_selected(colmn, offset)
    }

    pub fn sub_column_in_its_range_for_selected(
        &mut self,
        colmn: &str,
        offset: f64,
    ) -> Result<(), TGAGUIError> {
        self.offset_column_in_its_range_for_selected(colmn, -offset)
    }
    /// нужно чтобы график был выбран
    pub fn offset_y_column_in_its_range_for_experiment(
        &mut self,
        id: &str,
        offset: f64,
        from: f64,
        to: f64,
    ) -> Result<(), TGAGUIError> {
        self.series
            .offset_y_column_in_its_range(id, offset, from, to)?;

        Ok(())
    }

    pub fn offset_y_column_in_its_range_for_selected(
        &mut self,
        offset: f64,
    ) -> Result<(), TGAGUIError> {
        let id = self.get_experiment_by_selected_curve()?;
        if let Some(range) = self.interaction.selection_rect {
            let bounds = range.bounds();

            let y_min = bounds.2;
            let y_max = bounds.3;
            self.offset_y_column_in_its_range_for_experiment(&id, offset, y_min, y_max)?;
            self.create_points_for_curve(&id)?;
            self.reset_view();
        } else {
            return Err(TGAGUIError::BindingError(
                "No selection rectangle. Drag-select a region on the plot first.".to_string(),
            ));
        }

        Ok(())
    }
    pub fn add_y_column_in_its_range_for_selected(
        &mut self,
        offset: f64,
    ) -> Result<(), TGAGUIError> {
        self.offset_y_column_in_its_range_for_selected(offset)
    }

    pub fn sub_y_column_in_its_range_for_selected(
        &mut self,
        offset: f64,
    ) -> Result<(), TGAGUIError> {
        self.offset_y_column_in_its_range_for_selected(-offset)
    }

    // * AND /
    //=========================================================
    pub fn scale_columns(
        &mut self,
        id: &str,
        cols: &[&str],
        factor: f64,
    ) -> Result<(), TGAGUIError> {
        self.series.scale_columns(id, cols, factor)?;
        self.create_points_for_curve(&id)?;
        self.reset_view();
        Ok(())
    }

    pub fn scale_column(&mut self, id: &str, col: &str, factor: f64) -> Result<(), TGAGUIError> {
        self.scale_columns(id, &[col], factor)?;
        Ok(())
    }

    pub fn scale_column_of_selected(&mut self, col: &str, factor: f64) -> Result<(), TGAGUIError> {
        let id = self.get_experiment_by_selected_curve()?;
        debug!("for selected experiment {}", id);
        self.scale_column(&id, col, factor)?;
        Ok(())
    }

    pub fn mul_column_of_selected(&mut self, col: &str, factor: f64) -> Result<(), TGAGUIError> {
        info!("multiplication for column {} by factor {}", col, factor);
        self.scale_column_of_selected(col, factor)
    }

    pub fn div_column_of_selected(&mut self, col: &str, factor: f64) -> Result<(), TGAGUIError> {
        info!("division for column {} by factor {}", col, factor);
        self.scale_column_of_selected(col, 1.0 / factor)
    }

    // SCALE ONLY SELECTED
    pub fn scale_column_in_its_range_for_experiment(
        &mut self,
        id: &str,
        colmn: &str,
        scale: f64,
        from: f64,
        to: f64,
    ) -> Result<(), TGAGUIError> {
        self.series
            .scale_column_in_its_range(id, colmn, scale, from, to)?;

        Ok(())
    }
    pub fn scale_column_in_its_range_for_selected(
        &mut self,
        colmn: &str,
        scale: f64,
    ) -> Result<(), TGAGUIError> {
        let id = self.get_experiment_by_selected_curve()?;
        if let Some(range) = self.interaction.selection_rect {
            let bounds = range.bounds();

            let x_min = bounds.0;
            let x_max = bounds.1;
            let x_col = self.series.oneframeplot_axis_name(&id, XY::X)?;
            self.series
                .scale_column_in_range_by_reference(&id, colmn, &x_col, scale, x_min, x_max)?;
            self.create_points_for_curve(&id)?;
            self.reset_view();
        } else {
            return Err(TGAGUIError::BindingError(
                "No selection rectangle. Drag-select a region on the plot first.".to_string(),
            ));
        }

        Ok(())
    }
    pub fn mul_column_in_its_range_for_selected(
        &mut self,
        colmn: &str,
        scale: f64,
    ) -> Result<(), TGAGUIError> {
        self.scale_column_in_its_range_for_selected(colmn, scale)
    }

    pub fn div_column_in_its_range_for_selected(
        &mut self,
        colmn: &str,
        scale: f64,
    ) -> Result<(), TGAGUIError> {
        self.scale_column_in_its_range_for_selected(colmn, 1.0 / scale)
    }
    /// нужно чтобы график был выбран
    pub fn scale_y_column_in_its_range_for_experiment(
        &mut self,
        id: &str,
        scale: f64,
        from: f64,
        to: f64,
    ) -> Result<(), TGAGUIError> {
        self.series
            .scale_y_column_in_its_range(id, scale, from, to)?;

        Ok(())
    }

    pub fn scale_y_column_in_its_range_for_selected(
        &mut self,
        scale: f64,
    ) -> Result<(), TGAGUIError> {
        let id = self.get_experiment_by_selected_curve()?;
        if let Some(range) = self.interaction.selection_rect {
            let bounds = range.bounds();

            let y_min = bounds.2;
            let y_max = bounds.3;
            self.scale_y_column_in_its_range_for_experiment(&id, scale, y_min, y_max)?;
            self.create_points_for_curve(&id)?;
            self.reset_view();
        } else {
            return Err(TGAGUIError::BindingError(
                "No selection rectangle. Drag-select a region on the plot first.".to_string(),
            ));
        }

        Ok(())
    }
    pub fn mul_y_column_in_its_range_for_selected(
        &mut self,
        scale: f64,
    ) -> Result<(), TGAGUIError> {
        self.scale_y_column_in_its_range_for_selected(scale)
    }

    pub fn div_y_column_in_its_range_for_selected(
        &mut self,
        scale: f64,
    ) -> Result<(), TGAGUIError> {
        self.scale_y_column_in_its_range_for_selected(1.0 / scale)
    }
    // exp and ln
    //===========================================================================================
    pub fn exp_column(&mut self, id: &str, col_name: &str) -> Result<(), TGAGUIError> {
        self.series.exp_column(id, col_name)?;
        self.create_points_for_curve(&id)?;
        self.reset_view();
        Ok(())
    }

    pub fn exp_column_for_selected(&mut self, col_name: &str) -> Result<(), TGAGUIError> {
        let id = self.get_experiment_by_selected_curve()?;
        debug!("selected experiment {}", id);
        self.exp_column(&id, col_name)?;

        Ok(())
    }

    pub fn ln_column(&mut self, id: &str, col_name: &str) -> Result<(), TGAGUIError> {
        self.series.ln_column(id, col_name)?;
        self.create_points_for_curve(&id)?;
        self.reset_view();
        Ok(())
    }

    pub fn ln_column_for_selected(&mut self, col_name: &str) -> Result<(), TGAGUIError> {
        let id = self.get_experiment_by_selected_curve()?;
        debug!("selected experiment {}", id);
        self.ln_column(&id, col_name)?;
        Ok(())
    }

    //===========================================================================
    // SMOOTHING, SPLINING, FILTERING
    //==========================================================================
    // Hampel filter
    pub fn hampel_filter_for_selected(
        &mut self,
        col_name: &str,
        window: usize,
        sigma: f64,
        strategy: HampelStrategy,
    ) -> Result<(), TGAGUIError> {
        let id = self.get_experiment_by_selected_curve()?;
        let old_cols = self.list_of_columns(&id)?;
        self.series
            .hampel_filter(&id, col_name, window, sigma, strategy)?;
        self.create_plots_with_new_columns(&id, old_cols, None)?;
        Ok(())
    }
    pub fn hampel_filter_for_selected_as(
        &mut self,
        col_name: &str,
        window: usize,
        sigma: f64,
        strategy: HampelStrategy,
        out_col: Option<&str>,
    ) -> Result<(), TGAGUIError> {
        let id = self.get_experiment_by_selected_curve()?;
        let old_cols = self.list_of_columns(&id)?;
        self.series
            .hampel_filter_as(&id, col_name, window, sigma, strategy, out_col)?;
        self.create_plots_with_new_columns(&id, old_cols, None)?;
        Ok(())
    }
    // Savitzky Golay
    pub fn sg_filter_column_for_selected(
        &mut self,
        col: &str,
        window: usize,
        poly_order: usize,
        deriv: usize,
        delta: f64,
    ) -> Result<(), TGAGUIError> {
        let id = self.get_experiment_by_selected_curve()?;
        let old_cols = self.list_of_columns(&id)?;
        self.series
            .sg_filter_column(&id, col, window, poly_order, deriv, delta)?;
        self.create_plots_with_new_columns(&id, old_cols, None)?;

        Ok(())
    }
    pub fn sg_filter_column_for_selected_as(
        &mut self,
        col: &str,
        window: usize,
        poly_order: usize,
        deriv: usize,
        delta: f64,
        out_col: Option<&str>,
    ) -> Result<(), TGAGUIError> {
        let id = self.get_experiment_by_selected_curve()?;
        let old_cols = self.list_of_columns(&id)?;
        self.series
            .sg_filter_column_as(&id, col, window, poly_order, deriv, delta, out_col)?;
        self.create_plots_with_new_columns(&id, old_cols, None)?;
        Ok(())
    }
    // rolling mean
    pub fn rolling_mean_for_selected(
        &mut self,
        col_name: &str,
        window: usize,
    ) -> Result<(), TGAGUIError> {
        let id = self.get_experiment_by_selected_curve()?;
        let old_cols = self.list_of_columns(&id)?;
        self.series.rolling_mean(&id, col_name, window)?;
        self.series.trim_null_edges(&id)?;
        self.create_plots_with_new_columns(&id, old_cols, None)?;
        Ok(())
    }

    pub fn rolling_mean_for_selected_as(
        &mut self,
        col_name: &str,
        window: usize,
        out_col: Option<&str>,
    ) -> Result<(), TGAGUIError> {
        let id = self.get_experiment_by_selected_curve()?;
        let old_cols = self.list_of_columns(&id)?;
        self.series
            .rolling_mean_as(&id, col_name, window, out_col)?;
        self.series.trim_null_edges(&id)?;
        self.create_plots_with_new_columns(&id, old_cols, None)?;
        Ok(())
    }
    pub fn lowess_smooth_columns_for_selected(
        &mut self,
        time_col: &str,
        columns: &[&str],
        config: LowessConfig,
    ) -> Result<(), TGAGUIError> {
        let id = self.get_experiment_by_selected_curve()?;
        let old_cols = self.list_of_columns(&id)?;
        self.series
            .lowess_smooth_columns(&id, time_col, columns, config)?;
        self.create_plots_with_new_columns(&id, old_cols, Some(time_col))?;
        Ok(())
    }

    pub fn lowess_smooth_columns_for_selected_as(
        &mut self,
        time_col: &str,
        columns: &[&str],
        out_columns: &[Option<&str>],
        config: LowessConfig,
    ) -> Result<(), TGAGUIError> {
        let id = self.get_experiment_by_selected_curve()?;
        let old_cols = self.list_of_columns(&id)?;
        self.series
            .lowess_smooth_columns_as(&id, time_col, columns, out_columns, config)?;
        self.create_plots_with_new_columns(&id, old_cols, Some(time_col))?;
        Ok(())
    }
    // SPLINES

    pub fn splines_for_selected(
        &mut self,
        new_time_col: &str,
        n_points: usize,
        kind: SplineKind,
    ) -> Result<(), TGAGUIError> {
        let id = self.get_experiment_by_selected_curve()?;
        let old_cols = self.list_of_columns(&id)?;
        self.series
            .spline_resample_oneframeplot(&id, new_time_col, n_points, kind)?;
        self.create_plots_with_new_columns(&id, old_cols, Some(new_time_col))?;
        Ok(())
    }

    pub fn splines_for_selected_as(
        &mut self,
        new_time_col: &str,
        n_points: usize,
        kind: SplineKind,
        new_col: Option<&str>,
    ) -> Result<(), TGAGUIError> {
        let id = self.get_experiment_by_selected_curve()?;
        let old_cols = self.list_of_columns(&id)?;
        self.series
            .spline_resample_oneframeplot_as(&id, new_time_col, n_points, kind, new_col)?;
        self.create_plots_with_new_columns(&id, old_cols, Some(new_time_col))?;
        Ok(())
    }

    pub fn lsq_spline_resample_columns(
        &mut self,

        time_col: &str,
        new_time_col: &str,
        columns: &[&str],
        n_points: usize,
    ) -> Result<(), TGAGUIError> {
        let id = self.get_experiment_by_selected_curve()?;
        let old_cols = self.list_of_columns(&id)?;
        self.series
            .lsq_spline_resample_columns(&id, time_col, new_time_col, columns, n_points)?;
        self.create_plots_with_new_columns(&id, old_cols, Some(new_time_col))?;
        Ok(())
    }

    pub fn lsq_spline_resample_columns_as(
        &mut self,

        time_col: &str,
        new_time_col: &str,
        columns: &[&str],
        out_columns: &[Option<&str>],
        n_points: usize,
        degree: usize,
        n_internal_knots: usize,
    ) -> Result<(), TGAGUIError> {
        let id = self.get_experiment_by_selected_curve()?;
        let old_cols = self.list_of_columns(&id)?;
        self.series.lsq_spline_resample_columns_as(
            &id,
            time_col,
            new_time_col,
            columns,
            out_columns,
            n_points,
            degree,
            n_internal_knots,
            SolverKind::Banded,
        )?;
        self.create_plots_with_new_columns(&id, old_cols, Some(new_time_col))?;
        Ok(())
    }
    //================================================================================================
    //          RATES
    //==================================================================================================

    pub fn derive_rate(
        &mut self,
        id: &str,
        source_col: &str,
        new_col: &str,
        new_col_nature: ColumnNature,
    ) -> Result<(), TGAGUIError> {
        self.series
            .derive_rate(id, source_col, new_col, new_col_nature)?;
        Ok(())
    }

    pub fn derive_mass_rate(&mut self, id: &str, new_col: Option<&str>) -> Result<(), TGAGUIError> {
        let new_col = new_col.unwrap_or("dm_dt");
        self.series.derive_mass_rate(id, new_col)?;
        Ok(())
    }

    pub fn derive_temperature_rate(
        &mut self,
        id: &str,
        new_col: Option<&str>,
    ) -> Result<(), TGAGUIError> {
        let new_col = new_col.unwrap_or("dT_dt");
        self.series.derive_temperature_rate(id, new_col)?;
        Ok(())
    }

    pub fn derive_dimensionless_rate(
        &mut self,
        id: &str,
        col_name: &str,
        out_name: &str,
    ) -> Result<(), TGAGUIError> {
        self.series
            .derive_dimensionless_rate(id, col_name, out_name)?;
        Ok(())
    }

    pub fn derive_deta_dt(&mut self, id: &str, out_name: Option<&str>) -> Result<(), TGAGUIError> {
        let out_name = out_name.unwrap_or("deta_dt");
        self.series.derive_deta_dt(id, out_name)?;
        Ok(())
    }

    pub fn derive_dalpha_dt(
        &mut self,
        id: &str,
        out_name: Option<&str>,
    ) -> Result<(), TGAGUIError> {
        let out_name = out_name.unwrap_or("dalpha_dt");
        self.series.derive_dalpha_dt(id, out_name)?;
        Ok(())
    }

    pub fn derive_rate_for_selected(
        &mut self,

        source_col: &str,
        new_col: &str,
        new_col_nature: ColumnNature,
    ) -> Result<(), TGAGUIError> {
        let id = self.get_experiment_by_selected_curve()?;
        self.derive_rate(&id, source_col, new_col, new_col_nature)?;

        Ok(())
    }

    pub fn derive_mass_rate_for_selected(
        &mut self,
        new_col: Option<&str>,
    ) -> Result<(), TGAGUIError> {
        let id = self.get_experiment_by_selected_curve()?;
        self.derive_mass_rate(&id, new_col)?;
        Ok(())
    }

    pub fn derive_temperature_rate_for_selected(
        &mut self,
        new_col: Option<&str>,
    ) -> Result<(), TGAGUIError> {
        let id = self.get_experiment_by_selected_curve()?;
        self.derive_temperature_rate(&id, new_col)?;
        Ok(())
    }

    pub fn derive_dimensionless_rate_for_selected(
        &mut self,

        col_name: &str,
        out_name: &str,
    ) -> Result<(), TGAGUIError> {
        let id = self.get_experiment_by_selected_curve()?;
        self.derive_dimensionless_rate(&id, col_name, out_name)?;
        Ok(())
    }
    // convesion rate
    pub fn derive_deta_dt_for_selected(
        &mut self,
        out_name: Option<&str>,
    ) -> Result<(), TGAGUIError> {
        let id = self.get_experiment_by_selected_curve()?;
        self.derive_deta_dt(&id, out_name)?;
        Ok(())
    }
    // dimesionless mass rate
    pub fn derive_dalpha_dt_for_selected(
        &mut self,
        out_name: Option<&str>,
    ) -> Result<(), TGAGUIError> {
        let id = self.get_experiment_by_selected_curve()?;
        self.derive_dalpha_dt(&id, out_name)?;
        Ok(())
    }
    // AVERAGES
    /// Compute mean of `value_col` in the time interval `[from, to]`.
    pub fn mean_on_interval(
        &self,
        id: &str,
        value_col: &str,
        time_col: &str,
        from: f64,
        to: f64,
    ) -> Result<f64, TGAGUIError> {
        let r = self
            .series
            .mean_on_interval(id, value_col, time_col, from, to)?;
        Ok(r)
    }

    /// Compute mean on a column constrained by the same column's range.
    pub fn mean_on_interval_on_own_range(
        &self,
        id: &str,
        col: &str,
        from: f64,
        to: f64,
    ) -> Result<f64, TGAGUIError> {
        let r = self
            .series
            .mean_on_interval_on_own_range(id, col, from, to)?;
        Ok(r)
    }

    pub fn mean_on_interval_on_own_range_for_selected(&self) -> Result<f64, TGAGUIError> {
        let id = self.get_experiment_by_selected_curve()?;
        let r = if let Some(range) = self.interaction.selection_rect {
            let bounds = range.bounds();

            let y_min = bounds.2;
            let y_max = bounds.3;
            let y_col = self.series.oneframeplot_axis_name(&id, XY::Y)?;
            let r = self.mean_on_interval_on_own_range(&id, &y_col, y_min, y_max)?;
            r
        } else {
            return Err(TGAGUIError::BindingError(
                "No selection rectangle. Drag-select a region on the plot first.".to_string(),
            ));
        };
        Ok(r)
    }
    /// Compute mean of all entries in a column.
    pub fn mean_on_column(&self, id: &str, col: &str) -> Result<f64, TGAGUIError> {
        let r = self.series.mean_on_column(id, col)?;
        Ok(r)
    }
    pub fn mean_on_column_for_selected(&self) -> Result<f64, TGAGUIError> {
        let id = self.get_experiment_by_selected_curve()?;
        let y_col = self.series.oneframeplot_axis_name(&id, XY::Y)?;
        let r = self.mean_on_column(&id, &y_col)?;
        Ok(r)
    }
    // GOLDEN PIPELINE
    pub fn apply_golden_pipeline(
        &mut self,
        id: &str,
        config: GoldenPipelineConfig,
    ) -> Result<(), TGAGUIError> {
        let mut config = config;
        let do_we_need_new_exp = config.save_to_new_experiment;
        let line = self.settings.calibration_line().unwrap();
        let k = line.k();
        let b = line.b();
        config.b = b;
        config.k = k;
        self.series.apply_golden_pipeline(id, config)?;

        if do_we_need_new_exp {
            self.delete_all_plots()?;

            let last_exp = &self.series.experiments.last().unwrap();
            let new_id = last_exp.meta.id.clone();
            let deta_dt = &last_exp.dataset.schema.deta_dt.clone().unwrap();
            let eta = &last_exp.dataset.schema.eta.clone().unwrap();
            self.set_x(&new_id, &eta)?;
            self.set_y(&new_id, &deta_dt)?;
            self.create_points_for_curve(&new_id)?;
            info!("new experiment created");
        }
        Ok(())
    }
    //===========================================================================================
    //  KINETIC METHODS
    //===========================================================================================
    /// Build a materialized KineticDataView from a subset of experiments and selected column natures.
    /// This is a GUI-facing wrapper around the series-layer API.
    pub fn create_kinetic_data_view(
        &self,
        what_exp_to_take: Option<&[&str]>,
        what_cols_take: Vec<ColumnNature>,
    ) -> Result<KineticDataView, TGAGUIError> {
        let view = self
            .series
            .create_kinetic_data_view(what_exp_to_take, what_cols_take)?;
        Ok(view)
    }

    /// Build a KineticDataView using the column requirements of a selected isoconversional method.
    /// This is a GUI-facing wrapper around the series-layer API.
    pub fn create_kinetic_data_view_for_method(
        &self,
        what_exp_to_take: Option<&[&str]>,
        method: &IsoconversionalMethod,
    ) -> Result<KineticDataView, TGAGUIError> {
        let view = self
            .series
            .create_kinetic_data_view_for_method(what_exp_to_take, method)?;
        Ok(view)
    }

    /// Push an isoconversional result into TGASeries and build a fresh plot for it.
    ///
    /// Behavior:
    /// 1. Writes a new experiment into the series under the provided `id`.
    /// 2. Sets X = "eta" and Y = "Ea".
    /// 3. Hides all currently shown plots.
    /// 4. Creates a new curve in red and makes it visible.
    pub fn push_isoconversional_result(
        &mut self,
        result: &IsoconversionalResult,
        id: &str,
    ) -> Result<(), TGAGUIError> {
        self.series.push_isoconversional_result(result, id)?;

        // Make sure the new result is the only visible curve on the screen.
        for curve in &mut self.plots {
            curve.shown = false;
        }

        let x_col = "eta";
        let y_col = "Ea";

        self.set_x(id, x_col)?;
        self.set_y(id, y_col)?;

        self.create_points_for_curve_with_builder(id, Colours::Red)?;

        // Ensure the newly created curve is visible even if we replaced an older one.
        let short_name = self.create_plot_short_name(id, &x_col.to_string(), &y_col.to_string());
        if let Some(curve) = self
            .plots
            .iter_mut()
            .find(|curve| curve.plot_short_name == short_name)
        {
            curve.shown = true;
        }

        self.reset_view();

        Ok(())
    }

    /// Add a numeric column from Vec<f64> to an experiment by id.
    ///
    /// This is a GUI wrapper around the series-level method, which checks
    /// length, marks the column as NumericDerived, and logs the operation.
    pub fn add_column_from_vec(
        &mut self,
        id: &str,
        name: &str,
        unit: Unit,
        nature: ColumnNature,
        data: Vec<f64>,
    ) -> Result<(), TGAGUIError> {
        self.series
            .add_column_from_vec(id, name, unit, nature, data)?;
        Ok(())
    }

    /// Push fitted sublimation rates into the series and build plots for each experiment.
    ///
    /// For every `(conversion, fitted_conversion_rate)` pair we create a new curve:
    /// X = conversion column, Y = fitted rate column.
    pub fn push_sublimation_fitted_rates(
        &mut self,
        results: &[SublimationResult],
    ) -> Result<(), TGAGUIError> {
        SublimationMethod::push_fitted_rates_to_series(results, &mut self.series)?;

        for res in results {
            let conv_col = self
                .series
                .get_column_by_nature(&res.experiment_id, ColumnNature::Conversion)?
                .ok_or_else(|| {
                    TGAGUIError::BindingError(format!(
                        "Conversion column not found for experiment '{}'",
                        res.experiment_id
                    ))
                })?;

            self.set_x(&res.experiment_id, &conv_col)?;
            self.set_y(&res.experiment_id, "fitted_sublim_rate")?;
            self.create_points_for_curve_with_builder(&res.experiment_id, Colours::Red)?;
        }

        Ok(())
    }
}
