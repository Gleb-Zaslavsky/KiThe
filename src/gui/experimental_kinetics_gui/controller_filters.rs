use crate::Kinetics::experimental_kinetics::math::exp_kinetics_smooth_filter::HampelStrategy;
use crate::Kinetics::experimental_kinetics::one_experiment_dataset::Unit;
use crate::gui::experimental_kinetics_gui::model::{PlotModel, TGAGUIError};
use RustedSciThe::numerical::data_processing::lowess_wrapper::LowessConfig;
use RustedSciThe::numerical::data_processing::splines::SplineKind;
use std::collections::HashMap;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum SplineKindChoice {
    Linear,
    Cosine,
    Bezier,
}

impl Default for SplineKindChoice {
    fn default() -> Self {
        Self::Linear
    }
}

// Group: Mathematics
pub struct Mathematics {
    show_hampel_window: bool,
    show_sg_window: bool,
    show_rolling_window: bool,
    show_splines_window: bool,
    show_sql_splines_window: bool,
    show_lowess_window: bool,

    hampel_col: Option<String>,
    hampel_window: usize,
    hampel_sigma: f64,
    hampel_strategy: HampelStrategy,
    hampel_out_col: String,

    sg_col: Option<String>,
    sg_window: usize,
    sg_poly_order: usize,
    sg_deriv: usize,
    sg_delta: f64,
    sg_out_col: String,

    rolling_col: Option<String>,
    rolling_window: usize,
    rolling_out_col: String,

    splines_new_time_col: String,
    splines_n_points: usize,
    splines_kind: SplineKindChoice,
    splines_bezier_tension: f64,
    splines_out_col: String,

    sql_splines_degree: usize,
    sql_splines_n_points: usize,
    sql_splines_out_col: String,
    sq_splines_n_internal_knots: usize,

    lowess_time_col: Option<String>,
    lowess_fraction: f64,
    lowess_iterations: usize,
    lowess_delta: f64,
    lowess_parallel: bool,
    lowess_selected_cols: Vec<String>,
    lowess_out_cols: HashMap<String, String>,
}

impl Mathematics {
    /// Estimate the selected curve's time spacing from its currently displayed points.
    /// Plot downsampling can make this differ from the underlying measurement grid.
    fn time_grid_summary(model: &PlotModel) -> Option<(f64, &'static str, f64)> {
        let curve = model
            .get_selected_curve_index()
            .and_then(|index| model.plots.get(index))
            .or_else(|| {
                if model.is_only_one_shown() {
                    model.plots.iter().find(|curve| curve.is_shown())
                } else {
                    None
                }
            })?;
        let experiment = model
            .series
            .get_experiment_by_id(&curve.experiment_id)
            .ok()?;
        let unit = experiment.dataset.schema.columns.get(&curve.x_name)?.unit;
        let unit_name = match unit {
            Unit::Second => "s",
            Unit::Minute => "min",
            Unit::Hour => "h",
            _ => return None,
        };

        let mut steps: Vec<f64> = curve
            .points
            .windows(2)
            .map(|pair| (pair[1][0] - pair[0][0]).abs())
            .filter(|step| step.is_finite() && *step > 0.0)
            .collect();
        if steps.is_empty() {
            return None;
        }
        steps.sort_by(f64::total_cmp);
        let median_step = steps[steps.len() / 2];
        let total_span = curve
            .points
            .iter()
            .map(|point| point[0])
            .filter(|x| x.is_finite())
            .fold((f64::INFINITY, f64::NEG_INFINITY), |(min, max), x| {
                (min.min(x), max.max(x))
            });
        let span = total_span.1 - total_span.0;
        (span.is_finite() && span > 0.0).then_some((median_step, unit_name, span))
    }

    fn show_window_duration(ui: &mut egui::Ui, model: &PlotModel, points: usize) {
        if let Some((median_step, unit, _)) = Self::time_grid_summary(model) {
            let covered = median_step * points.saturating_sub(1) as f64;
            ui.small(format!(
                "Plot-based estimate: {covered:.6} {unit} ({points} points; median plotted step {median_step:.6} {unit}). Downsampling may make this differ from the measurement grid."
            ));
        } else {
            ui.small("Time span unavailable: select a curve with a bound time unit and at least two distinct time points.");
        }
    }

    pub(crate) fn close_window_on_layer(&mut self, layer: Option<egui::LayerId>) -> bool {
        let Some(layer) = layer else { return false };
        let windows = [
            ("Hampel filter", &mut self.show_hampel_window),
            ("Savitzky Golay", &mut self.show_sg_window),
            ("rolling mean", &mut self.show_rolling_window),
            ("Splines", &mut self.show_splines_window),
            ("LSQ Splines", &mut self.show_sql_splines_window),
            ("LOWESS", &mut self.show_lowess_window),
        ];
        for (title, is_open) in windows {
            if layer.id == egui::Id::unique(Some(title)) && *is_open {
                *is_open = false;
                return true;
            }
        }
        false
    }

    pub fn new() -> Self {
        Self {
            show_hampel_window: false,
            show_sg_window: false,
            show_rolling_window: false,
            show_splines_window: false,
            show_sql_splines_window: false,
            show_lowess_window: false,

            hampel_col: None,
            hampel_window: 11,
            hampel_sigma: 3.0,
            hampel_strategy: HampelStrategy::ReplaceWithMedian,
            hampel_out_col: String::new(),

            sg_col: None,
            sg_window: 11,
            sg_poly_order: 3,
            sg_deriv: 0,
            sg_delta: 1.0,
            sg_out_col: String::new(),

            rolling_col: None,
            rolling_window: 5,
            rolling_out_col: String::new(),

            splines_new_time_col: "t_spline".to_string(),
            splines_n_points: 1000,
            splines_kind: SplineKindChoice::Linear,
            splines_bezier_tension: 0.5,
            splines_out_col: String::new(),

            sql_splines_degree: 3,
            sql_splines_n_points: 1000,
            sql_splines_out_col: String::new(),
            sq_splines_n_internal_knots: 0,

            lowess_time_col: None,
            lowess_fraction: 0.5,
            lowess_iterations: 3,
            lowess_delta: 0.01,
            lowess_parallel: true,
            lowess_selected_cols: Vec::new(),
            lowess_out_cols: HashMap::new(),
        }
    }

    pub fn show_menu(&mut self, ui: &mut egui::Ui) {
        if ui.button("Hampel filter").on_hover_text("Robust local-median outlier filter; opens settings for the selected column and output behavior.").clicked() {
            self.show_hampel_window = true;
        }

        if ui.button(" Savitzky Golay").on_hover_text("Local polynomial smoother or derivative estimator; opens column, window, order and output settings.").clicked() {
            self.show_sg_window = true;
        }

        if ui
            .button("rolling mean")
            .on_hover_text(
                "Sliding-window mean smoother; opens source, window and output settings.",
            )
            .clicked()
        {
            self.show_rolling_window = true;
        }

        if ui
            .button("Splines")
            .on_hover_text(
                "Interpolate the selected curve onto a new time grid; creates resampled columns.",
            )
            .clicked()
        {
            self.show_splines_window = true;
        }

        if ui
            .button("LSQ Splines")
            .on_hover_text(
                "Fit a least-squares B-spline and resample the selected curve into new columns.",
            )
            .clicked()
        {
            self.show_sql_splines_window = true;
        }

        if ui.button("LOWESS").on_hover_text("Locally weighted regression smoother for selected columns versus a chosen time column.").clicked() {
            self.show_lowess_window = true;
        }
    }

    pub fn show_windows(&mut self, ctx: &egui::Context, model: &mut PlotModel) {
        self.show_hampel_window(ctx, model);
        self.show_sg_window(ctx, model);
        self.show_rolling_window(ctx, model);
        self.show_splines_window(ctx, model);
        self.show_lsq_splines_window(ctx, model);
        self.show_lowess_window(ctx, model);
    }

    fn columns_for_selected(model: &PlotModel) -> Result<Vec<String>, TGAGUIError> {
        model.list_of_columns_for_selected()
    }

    fn selected_column_ui(
        ui: &mut egui::Ui,
        label_text: &str,
        combo_id: &str,
        selected_col: &mut Option<String>,
        columns: &[String],
        default_col: Option<&str>,
    ) {
        if columns.is_empty() {
            ui.label("No columns available for selected curve.");
            return;
        }

        if selected_col
            .as_ref()
            .map(|col| !columns.contains(col))
            .unwrap_or(true)
        {
            *selected_col = default_col
                .filter(|default| columns.iter().any(|column| column == default))
                .map(str::to_owned)
                .or_else(|| columns.first().cloned());
        }

        let selected_text = selected_col
            .clone()
            .unwrap_or_else(|| "Select column".to_owned());
        ui.horizontal(|ui| {
            ui.label(label_text);
            let combo = egui::ComboBox::from_id_salt(combo_id)
                .selected_text(selected_text)
                .show_ui(ui, |ui| {
                    for col in columns {
                        ui.selectable_value(selected_col, Some(col.clone()), col);
                    }
                });
            let _ = combo.response.accessible_name(label_text);
        });
    }

    fn set_action_error(model: &mut PlotModel, err: TGAGUIError) {
        let text = Self::error_text(err);
        let _ = model.push_message(&text);
    }

    fn error_text(err: TGAGUIError) -> String {
        err.to_string()
    }

    fn finalize_success(model: &mut PlotModel, message: &str) {
        if let Err(err) = model.create_points_for_selected_curve() {
            if let Ok(id) = model.get_experiment_by_selected_curve() {
                model.plots.retain(|curve| curve.experiment_id != id);
                model.reset_view();
            }
            Self::set_action_error(
                model,
                TGAGUIError::post_commit_warning(format!(
                    "the filter was applied, but its graph could not be refreshed; stale curves were removed ({err})"
                )),
            );
            return;
        }
        model.reset_view();
        let _ = model.push_message(message);
    }

    fn show_hampel_window(&mut self, ctx: &egui::Context, model: &mut PlotModel) {
        if !self.show_hampel_window {
            return;
        }

        let mut open = self.show_hampel_window;
        egui::Window::new("Hampel filter")
            .open(&mut open)
            .resizable(false)
            .show(ctx, |ui| {
                ui.label("Robust outlier filter: compares each point with a local median and MAD-based threshold. Apply reads the selected column and settings; a different output name creates a new column, while the source name replaces it.");
                ui.add_space(4.0);
                match Self::columns_for_selected(model) {
                    Ok(columns) => {
                        let default_col = model.this_is_y().ok();
                        Self::selected_column_ui(
                            ui,
                            "Column",
                            "hampel_filter_column_combo",
                            &mut self.hampel_col,
                            &columns,
                            default_col.as_deref(),
                        );
                    }
                    Err(err) => {
                        ui.colored_label(egui::Color32::RED, Self::error_text(err));
                    }
                }

                ui.horizontal(|ui| {
                    let label = ui.label("Window size");
                    let value = ui.add(egui::DragValue::new(&mut self.hampel_window).speed(1.0));
                    value.labelled_by(label.id);
                });
                Self::show_window_duration(ui, model, self.hampel_window);
                ui.horizontal(|ui| {
                    let label = ui.label("Sigma threshold");
                    let value = ui.add(egui::DragValue::new(&mut self.hampel_sigma).speed(0.1));
                    value.labelled_by(label.id);
                });

                let strategy_text = match self.hampel_strategy {
                    HampelStrategy::ReplaceWithMedian => "ReplaceWithMedian",
                    HampelStrategy::ReplaceWithNaN => "ReplaceWithNaN",
                    HampelStrategy::Drop => "Drop",
                };
                egui::ComboBox::from_id_salt("hampel_strategy_combo")
                    .selected_text(strategy_text)
                    .show_ui(ui, |ui| {
                        if ui
                            .selectable_label(
                                matches!(self.hampel_strategy, HampelStrategy::ReplaceWithMedian),
                                "ReplaceWithMedian",
                            )
                            .clicked()
                        {
                            self.hampel_strategy = HampelStrategy::ReplaceWithMedian;
                        }
                        if ui
                            .selectable_label(
                                matches!(self.hampel_strategy, HampelStrategy::ReplaceWithNaN),
                                "ReplaceWithNaN",
                            )
                            .clicked()
                        {
                            self.hampel_strategy = HampelStrategy::ReplaceWithNaN;
                        }
                        if ui
                            .selectable_label(
                                matches!(self.hampel_strategy, HampelStrategy::Drop),
                                "Drop",
                            )
                            .clicked()
                        {
                            self.hampel_strategy = HampelStrategy::Drop;
                        }
                    })
                    .response
                    .accessible_name("Strategy");

                ui.horizontal(|ui| {
                    let label = ui.label("Output column (optional)");
                    let input = ui.text_edit_singleline(&mut self.hampel_out_col);
                    input.labelled_by(label.id);
                });

                ui.horizontal(|ui| {
                    if ui.button("Apply").on_hover_text("Requires a source column, window size and sigma threshold. A different output name creates a column; using the source name edits it in place.").clicked() {
                        let Some(col) = self.hampel_col.clone().or_else(|| model.this_is_y().ok())
                        else {
                            let _ = model.push_message("Select a column first.");
                            return;
                        };
                        let out_col = if self.hampel_out_col.trim().is_empty() {
                            None
                        } else {
                            Some(self.hampel_out_col.trim())
                        };

                        match model.hampel_filter_for_selected_as(
                            &col,
                            self.hampel_window.max(1),
                            self.hampel_sigma,
                            self.hampel_strategy.clone(),
                            out_col,
                        ) {
                            Ok(()) => Self::finalize_success(model, "Hampel filter applied"),
                            Err(err) => Self::set_action_error(model, err),
                        }
                    }
                    if ui.button("Close").clicked() {
                        self.show_hampel_window = false;
                    }
                });
            });
        self.show_hampel_window = open;
    }

    fn show_sg_window(&mut self, ctx: &egui::Context, model: &mut PlotModel) {
        if !self.show_sg_window {
            return;
        }

        let mut open = self.show_sg_window;
        egui::Window::new("Savitzky Golay")
            .open(&mut open)
            .resizable(false)
            .show(ctx, |ui| {
                ui.label("Local polynomial smoothing (or derivative estimation) with the Savitzky-Golay method. Apply reads the selected column and numeric settings; a different output name creates a new column, while the source name replaces it.");
                ui.add_space(4.0);
                match Self::columns_for_selected(model) {
                    Ok(columns) => {
                        let default_col = model.this_is_y().ok();
                        Self::selected_column_ui(
                            ui,
                            "Column",
                            "sg_filter_column_combo",
                            &mut self.sg_col,
                            &columns,
                            default_col.as_deref(),
                        );
                    }
                    Err(err) => {
                        ui.colored_label(egui::Color32::RED, Self::error_text(err));
                    }
                }

                ui.horizontal(|ui| {
                    let label = ui.label("Window size");
                    let value = ui.add(egui::DragValue::new(&mut self.sg_window).speed(1.0));
                    value.labelled_by(label.id);
                });
                Self::show_window_duration(ui, model, self.sg_window);
                ui.horizontal(|ui| {
                    let label = ui.label("Polynomial order");
                    let value = ui.add(egui::DragValue::new(&mut self.sg_poly_order).speed(1.0));
                    value.labelled_by(label.id);
                });
                ui.horizontal(|ui| {
                    let label = ui.label("Derivative order");
                    let value = ui.add(egui::DragValue::new(&mut self.sg_deriv).speed(1.0));
                    value.labelled_by(label.id);
                });
                ui.horizontal(|ui| {
                    let label = ui.label("Time step");
                    let value = ui.add(egui::DragValue::new(&mut self.sg_delta).speed(0.1));
                    value.labelled_by(label.id);
                });

                ui.horizontal(|ui| {
                    let label = ui.label("Output column (optional)");
                    let input = ui.text_edit_singleline(&mut self.sg_out_col);
                    input.labelled_by(label.id);
                });

                ui.horizontal(|ui| {
                    if ui.button("Apply").on_hover_text("Requires a source column and window, polynomial order, derivative order and time step. A different output name creates a column; using the source name edits it in place.").clicked() {
                        let Some(col) = self.sg_col.clone().or_else(|| model.this_is_y().ok())
                        else {
                            let _ = model.push_message("Select a column first.");
                            return;
                        };
                        let out_col = if self.sg_out_col.trim().is_empty() {
                            None
                        } else {
                            Some(self.sg_out_col.trim())
                        };

                        match model.sg_filter_column_for_selected_as(
                            &col,
                            self.sg_window.max(1),
                            self.sg_poly_order,
                            self.sg_deriv,
                            self.sg_delta,
                            out_col,
                        ) {
                            Ok(()) => Self::finalize_success(model, "Savitzky Golay applied"),
                            Err(err) => Self::set_action_error(model, err),
                        }
                    }
                    if ui.button("Close").clicked() {
                        self.show_sg_window = false;
                    }
                });
            });
        self.show_sg_window = open;
    }

    fn show_rolling_window(&mut self, ctx: &egui::Context, model: &mut PlotModel) {
        if !self.show_rolling_window {
            return;
        }

        let mut open = self.show_rolling_window;
        egui::Window::new("rolling mean")
            .open(&mut open)
            .resizable(false)
            .show(ctx, |ui| {
                ui.label("Sliding-window arithmetic mean for smoothing short-scale fluctuations. Apply reads the selected column and window size; a different output name creates a new column, while the source name replaces it.");
                ui.add_space(4.0);
                match Self::columns_for_selected(model) {
                    Ok(columns) => {
                        let default_col = model.this_is_y().ok();
                        Self::selected_column_ui(
                            ui,
                            "Column",
                            "rolling_mean_column_combo",
                            &mut self.rolling_col,
                            &columns,
                            default_col.as_deref(),
                        );
                    }
                    Err(err) => {
                        ui.colored_label(egui::Color32::RED, Self::error_text(err));
                    }
                }

                ui.horizontal(|ui| {
                    let label = ui.label("Window size");
                    let value = ui.add(egui::DragValue::new(&mut self.rolling_window).speed(1.0));
                    value.labelled_by(label.id);
                });
                Self::show_window_duration(ui, model, self.rolling_window);

                ui.horizontal(|ui| {
                    let label = ui.label("Output column (optional)");
                    let input = ui.text_edit_singleline(&mut self.rolling_out_col);
                    input.labelled_by(label.id);
                });

                ui.horizontal(|ui| {
                    if ui.button("Apply").on_hover_text("Requires a source column and window size. A different output name creates a column; using the source name edits it in place.").clicked() {
                        let Some(col) = self.rolling_col.clone().or_else(|| model.this_is_y().ok())
                        else {
                            let _ = model.push_message("Select a column first.");
                            return;
                        };
                        let out_col = if self.rolling_out_col.trim().is_empty() {
                            None
                        } else {
                            Some(self.rolling_out_col.trim())
                        };

                        match model.rolling_mean_for_selected_as(
                            &col,
                            self.rolling_window.max(1),
                            out_col,
                        ) {
                            Ok(()) => Self::finalize_success(model, "Rolling mean applied"),
                            Err(err) => Self::set_action_error(model, err),
                        }
                    }
                    if ui.button("Close").clicked() {
                        self.show_rolling_window = false;
                    }
                });
            });
        self.show_rolling_window = open;
    }

    fn show_lowess_window(&mut self, ctx: &egui::Context, model: &mut PlotModel) {
        if !self.show_lowess_window {
            return;
        }

        let mut open = self.show_lowess_window;
        egui::Window::new("LOWESS")
            .open(&mut open)
            .resizable(false)
            .show(ctx, |ui| {
                ui.label("LOWESS fits local weighted regressions against the chosen time column. Select one or more source columns and numeric settings; outputs are separate columns by default, and an output name matching its source replaces that source.");
                ui.add_space(4.0);
                match Self::columns_for_selected(model) {
                    Ok(columns) => {
                        if columns.is_empty() {
                            ui.label("No columns available for selected curve.");
                            return;
                        }

                        let default_col = model.this_is_x().ok();
                        Self::selected_column_ui(
                            ui,
                            "Time column",
                            "lowess_time_col_combo",
                            &mut self.lowess_time_col,
                            &columns,
                            default_col.as_deref(),
                        );

                        self.lowess_selected_cols
                            .retain(|c| columns.iter().any(|x| x == c));
                        self.lowess_out_cols
                            .retain(|k, _| columns.iter().any(|x| x == k));

                        ui.separator();
                        ui.label("Columns to smooth");

                        let selected_time = self.lowess_time_col.clone();
                        for col in columns.iter() {
                            if Some(col.as_str()) == selected_time.as_deref() {
                                continue;
                            }

                            let mut checked = self.lowess_selected_cols.iter().any(|c| c == col);
                            if ui.checkbox(&mut checked, col).changed() {
                                if checked {
                                    if !self.lowess_selected_cols.iter().any(|c| c == col) {
                                        self.lowess_selected_cols.push(col.clone());
                                    }
                                } else {
                                    self.lowess_selected_cols.retain(|c| c != col);
                                }
                            }

                            if checked {
                                let out_name = self
                                    .lowess_out_cols
                                    .entry(col.clone())
                                    .or_insert_with(String::new);
                                ui.horizontal(|ui| {
                                    let label = ui.label(format!("Output for {} (optional)", col));
                                    let input = ui.text_edit_singleline(out_name);
                                    let _ = input.labelled_by(label.id);
                                });
                            }
                        }
                    }
                    Err(err) => {
                        ui.colored_label(egui::Color32::RED, Self::error_text(err));
                        return;
                    }
                }

                ui.separator();
                ui.horizontal(|ui| {
                    let label = ui.label("Fraction");
                    let input = ui.add(egui::DragValue::new(&mut self.lowess_fraction).speed(0.01));
                    let _ = input.labelled_by(label.id);
                });
                ui.horizontal(|ui| {
                    let label = ui.label("Iterations");
                    let input =
                        ui.add(egui::DragValue::new(&mut self.lowess_iterations).speed(1.0));
                    let _ = input.labelled_by(label.id);
                });
                ui.horizontal(|ui| {
                    let label = ui.label("Delta");
                    let input = ui.add(egui::DragValue::new(&mut self.lowess_delta).speed(0.001));
                    let _ = input.labelled_by(label.id);
                });
                ui.checkbox(&mut self.lowess_parallel, "parallel");

                ui.horizontal(|ui| {
                    if ui.button("Apply").on_hover_text("Requires a time column, at least one selected column, and fraction/iteration/delta settings. Creates output columns by default; matching an output name to its source edits it in place.").clicked() {
                        let Some(time_col) = self.lowess_time_col.clone() else {
                            let _ = model.push_message("Select time_col first.");
                            return;
                        };
                        if self.lowess_selected_cols.is_empty() {
                            let _ = model.push_message("Select at least one column to smooth.");
                            return;
                        }

                        let mut source_cols = self.lowess_selected_cols.clone();
                        source_cols.retain(|c| c != &time_col);
                        if source_cols.is_empty() {
                            let _ = model
                                .push_message("Selected columns cannot contain only time_col.");
                            return;
                        }

                        let out_names_owned: Vec<Option<String>> = source_cols
                            .iter()
                            .map(|col| {
                                self.lowess_out_cols.get(col).and_then(|s| {
                                    let t = s.trim();
                                    if t.is_empty() {
                                        None
                                    } else {
                                        Some(t.to_string())
                                    }
                                })
                            })
                            .collect();

                        let col_refs: Vec<&str> = source_cols.iter().map(String::as_str).collect();
                        let out_refs: Vec<Option<&str>> =
                            out_names_owned.iter().map(|o| o.as_deref()).collect();

                        let cfg = LowessConfig {
                            fraction: self.lowess_fraction,
                            iterations: self.lowess_iterations,
                            delta: self.lowess_delta,
                            parallel: self.lowess_parallel,
                            ..LowessConfig::default()
                        };

                        match model.lowess_smooth_columns_for_selected_as(
                            &time_col, &col_refs, &out_refs, cfg,
                        ) {
                            Ok(()) => Self::finalize_success(model, "LOWESS smoothing applied"),
                            Err(err) => Self::set_action_error(model, err),
                        }
                    }
                    if ui.button("Close").clicked() {
                        self.show_lowess_window = false;
                    }
                });
            });
        self.show_lowess_window = open;
    }

    fn show_splines_window(&mut self, ctx: &egui::Context, model: &mut PlotModel) {
        if !self.show_splines_window {
            return;
        }

        let mut open = self.show_splines_window;
        egui::Window::new("Splines")
            .open(&mut open)
            .resizable(false)
            .show(ctx, |ui| {
                ui.label("Resamples the selected curve onto a new time grid using the chosen linear, cosine or Bezier interpolation. It creates a new time column and a resampled Y column; the original columns are retained.");
                ui.add_space(4.0);
                ui.horizontal(|ui| {
                    let label = ui.label("New time column");
                    let input = ui.text_edit_singleline(&mut self.splines_new_time_col);
                    let _ = input.labelled_by(label.id);
                });
                ui.horizontal(|ui| {
                    let label = ui.label("Number of points");
                    let input =
                        ui.add(egui::DragValue::new(&mut self.splines_n_points).speed(10.0));
                    let _ = input.labelled_by(label.id);
                });

                let kind_text = match self.splines_kind {
                    SplineKindChoice::Linear => "Linear",
                    SplineKindChoice::Cosine => "Cosine",
                    SplineKindChoice::Bezier => "Bezier",
                };

                egui::ComboBox::from_id_salt("splines_kind_combo")
                    .selected_text(kind_text)
                    .show_ui(ui, |ui| {
                        ui.selectable_value(
                            &mut self.splines_kind,
                            SplineKindChoice::Linear,
                            "Linear",
                        );
                        ui.selectable_value(
                            &mut self.splines_kind,
                            SplineKindChoice::Cosine,
                            "Cosine",
                        );
                        ui.selectable_value(
                            &mut self.splines_kind,
                            SplineKindChoice::Bezier,
                            "Bezier",
                        );
                    })
                    .response
                    .accessible_name("Spline type");

                if self.splines_kind == SplineKindChoice::Bezier {
                    ui.horizontal(|ui| {
                        let label = ui.label("Bezier tension");
                        let input = ui.add(
                            egui::DragValue::new(&mut self.splines_bezier_tension).speed(0.05),
                        );
                        let _ = input.labelled_by(label.id);
                    });
                }

                ui.horizontal(|ui| {
                    let label = ui.label("Output column (optional)");
                    let input = ui.text_edit_singleline(&mut self.splines_out_col);
                    let _ = input.labelled_by(label.id);
                });

                ui.horizontal(|ui| {
                    if ui.button("Apply").on_hover_text("Requires a new time-column name and point count; Bezier also uses tension. Creates resampled time and Y columns; original columns are retained.").clicked() {
                        let new_time_col = self.splines_new_time_col.trim();
                        if new_time_col.is_empty() {
                            let _ = model.push_message("new_time_col cannot be empty.");
                            return;
                        }
                        let out_col = if self.splines_out_col.trim().is_empty() {
                            None
                        } else {
                            Some(self.splines_out_col.trim())
                        };

                        let kind = match self.splines_kind {
                            SplineKindChoice::Linear => SplineKind::Linear,
                            SplineKindChoice::Cosine => SplineKind::Cosine,
                            SplineKindChoice::Bezier => {
                                SplineKind::Bezier(self.splines_bezier_tension)
                            }
                        };

                        match model.splines_for_selected_as(
                            new_time_col,
                            self.splines_n_points.max(2),
                            kind,
                            out_col,
                        ) {
                            Ok(()) => Self::finalize_success(model, "Spline resampling applied"),
                            Err(err) => Self::set_action_error(model, err),
                        }
                    }
                    if ui.button("Close").clicked() {
                        self.show_splines_window = false;
                    }
                });
            });
        self.show_splines_window = open;
    }

    fn show_lsq_splines_window(&mut self, ctx: &egui::Context, model: &mut PlotModel) {
        if !self.show_sql_splines_window {
            return;
        }

        let mut open = self.show_sql_splines_window;
        egui::Window::new("LSQ Splines")
            .open(&mut open)
            .resizable(false)
            .show(ctx, |ui| {
                ui.label("Fits a least-squares B-spline to the selected curve and resamples it. Creates a new time column and fitted Y column; original columns are retained.");
                ui.add_space(4.0);
                ui.horizontal(|ui| {
                    let label = ui.label("Degree");
                    let input =
                        ui.add(egui::DragValue::new(&mut self.sql_splines_degree).speed(1.0));
                    let _ = input.labelled_by(label.id);
                });
                ui.horizontal(|ui| {
                    let label = ui.label("Number of points");
                    let input =
                        ui.add(egui::DragValue::new(&mut self.sql_splines_n_points).speed(10.0));
                    let _ = input.labelled_by(label.id);
                });
                ui.horizontal(|ui| {
                    let label = ui.label("Internal knots");
                    let input = ui.add(
                        egui::DragValue::new(&mut self.sq_splines_n_internal_knots).speed(10.0),
                    );
                    let _ = input.labelled_by(label.id);
                });
                if let Some((_, unit, span)) = Self::time_grid_summary(model) {
                    if self.sq_splines_n_internal_knots == 0 {
                        ui.small("No internal knots; knot spacing is not applicable.");
                    } else {
                        let spacing = span / (self.sq_splines_n_internal_knots + 1) as f64;
                        ui.small(format!(
                            "Plot-based estimate of internal-knot spacing: {spacing:.6} {unit}. Downsampling may affect the displayed time range."
                        ));
                    }
                } else {
                    ui.small("Knot spacing unavailable: select a curve with a bound time unit and at least two distinct time points.");
                }

                ui.horizontal(|ui| {
                    let label = ui.label("Output column (optional)");
                    let input = ui.text_edit_singleline(&mut self.sql_splines_out_col);
                    let _ = input.labelled_by(label.id);
                });

                ui.horizontal(|ui| {
                    if ui.button("Apply").on_hover_text("Requires degree, point count and internal-knot count. Uses the selected curve's X/Y columns; creates new time and fitted Y columns, retaining the originals.").clicked() {
                        let time_col = match model.this_is_x() {
                            Ok(col) => col,
                            Err(err) => {
                                Self::set_action_error(model, err);
                                return;
                            }
                        };
                        let y_col = match model.this_is_y() {
                            Ok(col) => col,
                            Err(err) => {
                                Self::set_action_error(model, err);
                                return;
                            }
                        };

                        let out_col = if self.sql_splines_out_col.trim().is_empty() {
                            None
                        } else {
                            Some(self.sql_splines_out_col.trim())
                        };

                        match model.lsq_spline_resample_columns_as(
                            &time_col,
                            "t_lsq_spline",
                            &[&y_col],
                            &[out_col],
                            self.sql_splines_n_points.max(2),
                            self.sql_splines_degree,
                            self.sq_splines_n_internal_knots,
                        ) {
                            Ok(()) => {
                                Self::finalize_success(model, "LSQ Spline resampling applied")
                            }
                            Err(err) => Self::set_action_error(model, err),
                        }
                    }
                    if ui.button("Close").clicked() {
                        self.show_sql_splines_window = false;
                    }
                });
            });
        self.show_sql_splines_window = open;
    }
}
