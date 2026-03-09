use crate::Kinetics::experimental_kinetics::exp_kinetics_smooth_filter::HampelStrategy;
use crate::Kinetics::experimental_kinetics::lowess_wrapper::LowessConfig;
use crate::Kinetics::experimental_kinetics::splines::SplineKind;
use crate::gui::experimental_kinetics_gui::model::{PlotModel, TGAGUIError};
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
            sq_splines_n_internal_knots: 10000,

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
        if ui.button("Hampel filter").clicked() {
            self.show_hampel_window = true;
        }

        if ui.button(" Savitzky Golay").clicked() {
            self.show_sg_window = true;
        }

        if ui.button("rolling mean").clicked() {
            self.show_rolling_window = true;
        }

        if ui.button("Splines").clicked() {
            self.show_splines_window = true;
        }

        if ui.button("LSQ Splines").clicked() {
            self.show_sql_splines_window = true;
        }
        if ui.button("Differentiate").clicked() {
            println!("Stub: Differentiate clicked");
        }

        if ui.button("LOWESS").clicked() {
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
        combo_id: &str,
        selected_col: &mut Option<String>,
        columns: &[String],
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
            *selected_col = columns.first().cloned();
        }

        let selected_text = selected_col.as_deref().unwrap_or("Select column");
        egui::ComboBox::from_id_salt(combo_id)
            .selected_text(selected_text)
            .show_ui(ui, |ui| {
                for col in columns {
                    ui.selectable_value(selected_col, Some(col.clone()), col);
                }
            });
    }

    fn set_action_error(model: &mut PlotModel, err: TGAGUIError) {
        let text = Self::error_text(err);
        let _ = model.push_message(&text);
    }

    fn error_text(err: TGAGUIError) -> String {
        match err {
            TGAGUIError::TGADomainError(domain) => format!("{:?}", domain),
            TGAGUIError::SettingsErrors(msg) => msg,
            TGAGUIError::BindingError(msg) => msg,
        }
    }

    fn finalize_success(model: &mut PlotModel, message: &str) {
        if let Err(err) = model.create_points_for_selected_curve() {
            Self::set_action_error(model, err);
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
                match Self::columns_for_selected(model) {
                    Ok(columns) => {
                        ui.label("col_name");
                        Self::selected_column_ui(
                            ui,
                            "hampel_filter_column_combo",
                            &mut self.hampel_col,
                            &columns,
                        );
                    }
                    Err(err) => {
                        ui.colored_label(egui::Color32::RED, Self::error_text(err));
                    }
                }

                ui.horizontal(|ui| {
                    ui.label("window");
                    ui.add(egui::DragValue::new(&mut self.hampel_window).speed(1.0));
                });
                ui.horizontal(|ui| {
                    ui.label("sigma");
                    ui.add(egui::DragValue::new(&mut self.hampel_sigma).speed(0.1));
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
                    });

                ui.horizontal(|ui| {
                    ui.label("out_col (optional)");
                    ui.text_edit_singleline(&mut self.hampel_out_col);
                });

                ui.horizontal(|ui| {
                    if ui.button("Apply").clicked() {
                        let Some(col) = self.hampel_col.clone() else {
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
                match Self::columns_for_selected(model) {
                    Ok(columns) => {
                        ui.label("col");
                        Self::selected_column_ui(
                            ui,
                            "sg_filter_column_combo",
                            &mut self.sg_col,
                            &columns,
                        );
                    }
                    Err(err) => {
                        ui.colored_label(egui::Color32::RED, Self::error_text(err));
                    }
                }

                ui.horizontal(|ui| {
                    ui.label("window");
                    ui.add(egui::DragValue::new(&mut self.sg_window).speed(1.0));
                });
                ui.horizontal(|ui| {
                    ui.label("poly_order");
                    ui.add(egui::DragValue::new(&mut self.sg_poly_order).speed(1.0));
                });
                ui.horizontal(|ui| {
                    ui.label("deriv");
                    ui.add(egui::DragValue::new(&mut self.sg_deriv).speed(1.0));
                });
                ui.horizontal(|ui| {
                    ui.label("delta");
                    ui.add(egui::DragValue::new(&mut self.sg_delta).speed(0.1));
                });

                ui.horizontal(|ui| {
                    ui.label("out_col (optional)");
                    ui.text_edit_singleline(&mut self.sg_out_col);
                });

                ui.horizontal(|ui| {
                    if ui.button("Apply").clicked() {
                        let Some(col) = self.sg_col.clone() else {
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
                match Self::columns_for_selected(model) {
                    Ok(columns) => {
                        ui.label("col_name");
                        Self::selected_column_ui(
                            ui,
                            "rolling_mean_column_combo",
                            &mut self.rolling_col,
                            &columns,
                        );
                    }
                    Err(err) => {
                        ui.colored_label(egui::Color32::RED, Self::error_text(err));
                    }
                }

                ui.horizontal(|ui| {
                    ui.label("window");
                    ui.add(egui::DragValue::new(&mut self.rolling_window).speed(1.0));
                });

                ui.horizontal(|ui| {
                    ui.label("out_col (optional)");
                    ui.text_edit_singleline(&mut self.rolling_out_col);
                });

                ui.horizontal(|ui| {
                    if ui.button("Apply").clicked() {
                        let Some(col) = self.rolling_col.clone() else {
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
                match Self::columns_for_selected(model) {
                    Ok(columns) => {
                        if columns.is_empty() {
                            ui.label("No columns available for selected curve.");
                            return;
                        }

                        ui.label("time_col");
                        Self::selected_column_ui(
                            ui,
                            "lowess_time_col_combo",
                            &mut self.lowess_time_col,
                            &columns,
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
                                    ui.label(format!("out for {} (optional)", col));
                                    ui.text_edit_singleline(out_name);
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
                    ui.label("fraction");
                    ui.add(egui::DragValue::new(&mut self.lowess_fraction).speed(0.01));
                });
                ui.horizontal(|ui| {
                    ui.label("iterations");
                    ui.add(egui::DragValue::new(&mut self.lowess_iterations).speed(1.0));
                });
                ui.horizontal(|ui| {
                    ui.label("delta");
                    ui.add(egui::DragValue::new(&mut self.lowess_delta).speed(0.001));
                });
                ui.checkbox(&mut self.lowess_parallel, "parallel");

                ui.horizontal(|ui| {
                    if ui.button("Apply").clicked() {
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
                ui.horizontal(|ui| {
                    ui.label("new_time_col");
                    ui.text_edit_singleline(&mut self.splines_new_time_col);
                });
                ui.horizontal(|ui| {
                    ui.label("n_points");
                    ui.add(egui::DragValue::new(&mut self.splines_n_points).speed(10.0));
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
                    });

                if self.splines_kind == SplineKindChoice::Bezier {
                    ui.horizontal(|ui| {
                        ui.label("bezier_tension");
                        ui.add(egui::DragValue::new(&mut self.splines_bezier_tension).speed(0.05));
                    });
                }

                ui.horizontal(|ui| {
                    ui.label("out_col (optional)");
                    ui.text_edit_singleline(&mut self.splines_out_col);
                });

                ui.horizontal(|ui| {
                    if ui.button("Apply").clicked() {
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
                ui.horizontal(|ui| {
                    ui.label("degree");
                    ui.add(egui::DragValue::new(&mut self.sql_splines_degree).speed(1.0));
                });
                ui.horizontal(|ui| {
                    ui.label("n_points");
                    ui.add(egui::DragValue::new(&mut self.sql_splines_n_points).speed(10.0));
                });
                ui.horizontal(|ui| {
                    ui.label("n_internal_knots");
                    ui.add(egui::DragValue::new(&mut self.sq_splines_n_internal_knots).speed(10.0));
                });

                ui.horizontal(|ui| {
                    ui.label("out_col (optional)");
                    ui.text_edit_singleline(&mut self.sql_splines_out_col);
                });

                ui.horizontal(|ui| {
                    if ui.button("Apply").clicked() {
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
