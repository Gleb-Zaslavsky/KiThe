use crate::Kinetics::experimental_kinetics::kinetic_methods::is_this_a_sublimation::{
    SublimationMethod, SublimationResult,
};
use crate::Kinetics::experimental_kinetics::math::fitting::{
    FitColumnRequest, FitColumnResult, FitColumnSession, FittedEvaluationBatch,
};
use crate::Kinetics::experimental_kinetics::one_experiment_dataset::ColumnNature;
use crate::gui::experimental_kinetics_gui::model::{Colours, PlotModel, TGAGUIError};

impl PlotModel {
    /// Fit one `(x, y)` column pair, publish fitted y, and create its plot.
    ///
    /// The domain layer owns fitting semantics and provenance. The GUI layer only
    /// asks for the operation, then wires `x` vs `fitted_y` into the normal plot
    /// machinery so the result is immediately visible.
    pub fn fit_column_and_plot(
        &mut self,
        request: FitColumnRequest,
    ) -> Result<FitColumnResult, TGAGUIError> {
        self.settings
            .n_points()
            .filter(|count| *count > 0)
            .ok_or_else(|| {
                TGAGUIError::SettingsErrors(
                    "Plot sample count must be a positive integer.".to_string(),
                )
            })?;

        let result = self.series.fit_column(request)?;
        let short_name =
            self.create_plot_short_name(&result.experiment_id, &result.x_col, &result.output_col);
        let previous_curve = self
            .plots
            .iter()
            .find(|curve| curve.plot_short_name == short_name)
            .cloned();
        let previous_interaction = self.interaction.clone();
        let previous_x = self
            .series
            .oneframeplot_axis_name(
                &result.experiment_id,
                crate::Kinetics::experimental_kinetics::exp_engine_api::XY::X,
            )
            .ok();
        let previous_y = self
            .series
            .oneframeplot_axis_name(
                &result.experiment_id,
                crate::Kinetics::experimental_kinetics::exp_engine_api::XY::Y,
            )
            .ok();
        let plot_result = (|| -> Result<(), TGAGUIError> {
            self.set_x(&result.experiment_id, &result.x_col)?;
            self.set_y(&result.experiment_id, &result.output_col)?;
            self.create_points_for_curve_with_builder(&result.experiment_id, Colours::Red)?;
            if let Some(curve) = self
                .plots
                .iter_mut()
                .find(|curve| curve.plot_short_name == short_name)
            {
                curve.shown = true;
            }
            self.reset_view();
            Ok(())
        })();
        if let Err(error) = plot_result {
            self.plots
                .retain(|curve| curve.plot_short_name != short_name);
            if let Some(curve) = previous_curve {
                let index = self
                    .plots
                    .iter()
                    .position(|existing| existing.experiment_index > curve.experiment_index)
                    .unwrap_or(self.plots.len());
                self.plots.insert(index, curve);
            }
            if let Some(x) = previous_x {
                let _ = self.set_x(&result.experiment_id, &x);
            }
            if let Some(y) = previous_y {
                let _ = self.set_y(&result.experiment_id, &y);
            }
            self.interaction = previous_interaction;
            return Err(TGAGUIError::post_commit_warning(format!(
                "fit result '{}' was committed in experiment '{}', but its curve could not be plotted ({error})",
                result.output_col, result.experiment_id
            )));
        }
        Ok(result)
    }

    /// Fit without changing the source experiment; evaluation is a separate user action.
    pub fn fit_column_for_evaluation(
        &self,
        request: FitColumnRequest,
    ) -> Result<FitColumnSession, TGAGUIError> {
        Ok(self.series.fit_column_for_evaluation(request)?)
    }

    /// Commit an evaluated fit batch as a new experiment and display each fitted column.
    pub fn push_fitted_evaluation_batch_and_plot(
        &mut self,
        batch: &FittedEvaluationBatch,
        new_id: impl Into<String>,
    ) -> Result<(), TGAGUIError> {
        let new_id = new_id.into();
        let x_col = batch.x_col.clone();
        let y_columns = batch
            .columns
            .iter()
            .map(|column| column.name.clone())
            .collect::<Vec<_>>();
        let previous_interaction = self.interaction.clone();
        let previous_reset_view_requested = self.reset_view_requested;
        let previous_plot_recreation_required = self.plot_recreation_required;
        self.series
            .push_fitted_evaluation_batch(batch, new_id.clone())?;
        let plot_result = (|| -> Result<(), TGAGUIError> {
            for y_col in &y_columns {
                self.set_x(&new_id, &x_col)?;
                self.set_y(&new_id, y_col)?;
                self.create_points_for_curve_with_builder(&new_id, Colours::Red)?;
                let short_name = self.create_plot_short_name(&new_id, &x_col, y_col);
                if let Some(curve) = self
                    .plots
                    .iter_mut()
                    .find(|plot| plot.plot_short_name == short_name)
                {
                    curve.shown = true;
                }
            }
            self.reset_view();
            Ok(())
        })();
        if let Err(error) = plot_result {
            self.plots.retain(|plot| plot.experiment_id != new_id);
            self.interaction = previous_interaction;
            self.reset_view_requested = previous_reset_view_requested;
            self.plot_recreation_required = previous_plot_recreation_required;
            return Err(TGAGUIError::post_commit_warning(format!(
                "reconstructed experiment '{new_id}' was added, but its curves could not be plotted; stale curves were removed ({error})"
            )));
        }
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

        let plot_result = (|| -> Result<(), TGAGUIError> {
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
        })();

        if let Err(error) = plot_result {
            for result in results {
                self.plots.retain(|plot| {
                    plot.experiment_id != result.experiment_id
                        || plot.y_name != "fitted_sublim_rate"
                });
            }
            self.reset_view();
            return Err(TGAGUIError::post_commit_warning(format!(
                "fitted sublimation-rate data were added to the series, but not all result curves could be plotted; the data remain available ({error})"
            )));
        }

        Ok(())
    }
}
