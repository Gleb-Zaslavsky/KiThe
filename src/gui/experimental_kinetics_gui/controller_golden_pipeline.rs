//! Golden Pipeline dialog for configuring and applying the golden pipeline transformation.

use crate::Kinetics::experimental_kinetics::column_provenance::{
    ColumnTransformKind, QualityStatus,
};
use crate::Kinetics::experimental_kinetics::exp_engine_api::GoldenPipelineConfig;
use crate::Kinetics::experimental_kinetics::math::ndarray_statistics::FilterQualityReport;
use crate::gui::experimental_kinetics_gui::model::PlotModel;
use eframe::egui;
use log::info;

/// State for the golden pipeline dialog.
#[derive(Clone, Debug)]
pub struct GoldenPipelineDialogState {
    pub open: bool,
    pub config: GoldenPipelineConfig,
    pub last_error: Option<String>,
    pub quality_report: Option<FilterQualityReport>,
    pub quality_note: Option<String>,
}

impl Default for GoldenPipelineDialogState {
    fn default() -> Self {
        Self {
            open: false,
            config: GoldenPipelineConfig::default(),
            last_error: None,
            quality_report: None,
            quality_note: None,
        }
    }
}

impl GoldenPipelineDialogState {
    /// Open the dialog with default config.
    pub fn open(&mut self) {
        self.open = true;
        self.config = GoldenPipelineConfig::default();
        self.last_error = None;
        self.quality_report = None;
        self.quality_note = None;
    }

    fn find_savitzky_golay_report(
        model: &PlotModel,
        experiment_id: &str,
    ) -> Option<FilterQualityReport> {
        let experiment = model.series.get_experiment_by_id(experiment_id).ok()?;
        for column in experiment.dataset.list_of_columns() {
            let Some(provenance) = experiment.dataset.column_provenance(&column) else {
                continue;
            };
            for step in provenance.steps.iter().rev() {
                if matches!(&step.kind, ColumnTransformKind::SavitzkyGolay { .. }) {
                    if let Some(quality) = &step.quality {
                        if matches!(&quality.status, QualityStatus::Computed) {
                            if let Some(report) = &quality.report {
                                return Some(report.clone());
                            }
                        }
                    }
                }
            }
        }
        None
    }

    /// Show the dialog window.
    pub fn show(&mut self, ctx: &egui::Context, model: &mut PlotModel) {
        if !self.open {
            return;
        }

        let mut dialog_open = self.open;
        let mut should_apply = false;
        let mut should_close = false;

        egui::Window::new("Golden Pipeline")
            .open(&mut dialog_open)
            .resizable(true)
            .default_size([500.0, 600.0])
            .show(ctx, |ui| {
                ui.heading("Golden Pipeline Configuration");
                ui.separator();

                egui::ScrollArea::vertical().show(ui, |ui| {
                    // Time cuts
                    ui.collapsing("Time cuts", |ui| {
                        ui.horizontal(|ui| {
                            ui.label("Cut before time (optional):");
                            if let Some(ref mut val) = self.config.time_cut_before {
                                ui.add(egui::DragValue::new(val).speed(0.1));
                            } else {
                                if ui.button("Set").clicked() {
                                    self.config.time_cut_before = Some(0.0);
                                }
                            }
                            if self.config.time_cut_before.is_some() && ui.button("Clear").clicked()
                            {
                                self.config.time_cut_before = None;
                            }
                        });
                        ui.horizontal(|ui| {
                            ui.label("Cut after time (optional):");
                            if let Some(ref mut val) = self.config.time_cut_after {
                                ui.add(egui::DragValue::new(val).speed(0.1));
                            } else {
                                if ui.button("Set").clicked() {
                                    self.config.time_cut_after = Some(0.0);
                                }
                            }
                            if self.config.time_cut_after.is_some() && ui.button("Clear").clicked()
                            {
                                self.config.time_cut_after = None;
                            }
                        });
                    });

                    // Savitzky–Golay smoothing
                    ui.collapsing("Savitzky–Golay Smoothing", |ui| {
                        ui.horizontal(|ui| {
                            ui.label("Window size:");
                            ui.add(
                                egui::DragValue::new(&mut self.config.sav_gol_config.window_size)
                                    .speed(1)
                                    .range(3..=101),
                            );
                        });
                        ui.horizontal(|ui| {
                            ui.label("Polynomial degree:");
                            ui.add(
                                egui::DragValue::new(&mut self.config.sav_gol_config.poly_degree)
                                    .speed(1)
                                    .range(1..=10),
                            );
                        });
                        ui.horizontal(|ui| {
                            ui.label("Derivative order:");
                            ui.add(
                                egui::DragValue::new(&mut self.config.sav_gol_config.deriv)
                                    .speed(1)
                                    .range(0..=5),
                            );
                        });
                        ui.horizontal(|ui| {
                            ui.label("Delta (step):");
                            ui.add(
                                egui::DragValue::new(&mut self.config.sav_gol_config.delta)
                                    .speed(0.01),
                            );
                        });
                    });

                    // Averaging time for conversion
                    ui.collapsing("Conversion", |ui| {
                        ui.horizontal(|ui| {
                            ui.label("Averaging time (for conversion) - IN SECONDS:");
                            ui.add(
                                egui::DragValue::new(&mut self.config.averaging_time).speed(0.1),
                            );
                        });
                    });

                    // Spline resampling
                    ui.collapsing("Spline Resampling", |ui| {
                        ui.horizontal(|ui| {
                            ui.label("Number of points:");
                            ui.add(
                                egui::DragValue::new(&mut self.config.spline_config.n_points)
                                    .speed(1)
                                    .range(10..=10000),
                            );
                        });
                        ui.horizontal(|ui| {
                            ui.label("Number of internal knots:");
                            ui.add(
                                egui::DragValue::new(
                                    &mut self.config.spline_config.n_internal_points,
                                )
                                .speed(1)
                                .range(10..=100000),
                            );
                        });
                        ui.horizontal(|ui| {
                            ui.label("Spline degree:");
                            ui.add(
                                egui::DragValue::new(&mut self.config.spline_config.degree)
                                    .speed(1)
                                    .range(1..=5),
                            );
                        });
                    });

                    // Output options
                    ui.collapsing("Output Options", |ui| {
                        ui.checkbox(
                            &mut self.config.save_to_new_experiment,
                            "Save to new experiment",
                        );
                        ui.checkbox(&mut self.config.del_old_experiment, "Delete old experiment");
                    });
                });

                ui.separator();

                if let Some(err) = &self.last_error {
                    ui.colored_label(egui::Color32::RED, err);
                }

                if let Some(report) = &self.quality_report {
                    ui.separator();
                    ui.strong("Savitzky-Golay mass smoothing quality");
                    ui.label(format!(
                        "Compared: {} -> {}",
                        report.column_raw, report.column_filtered
                    ));
                    egui::Grid::new("golden_pipeline_quality_report")
                        .num_columns(2)
                        .striped(true)
                        .show(ui, |ui| {
                            ui.label("First-difference roughness ratio");
                            ui.label(format!("{:.6}", report.roughness_ratio_1st));
                            ui.end_row();
                            ui.label("Second-difference roughness ratio");
                            ui.label(format!("{:.6}", report.roughness_ratio_2nd));
                            ui.end_row();
                            ui.label("RMSE vs input (mass units)");
                            ui.label(format!("{:.6}", report.rmse));
                            ui.end_row();
                            ui.label("Normalized RMSE");
                            ui.label(format!("{:.6}", report.normalized_rmse));
                            ui.end_row();
                            ui.label("Correlation");
                            ui.label(format!("{:.6}", report.correlation));
                            ui.end_row();
                        });
                    ui.label(
                        "Roughness ratios below 1 mean smoother differences; these metrics compare with the input, not with an unknown true signal.",
                    );
                }
                if let Some(note) = &self.quality_note {
                    ui.label(note);
                }

                ui.horizontal(|ui| {
                    if ui.button("Close").clicked() {
                        should_close = true;
                    }
                    if ui.button("Calculate").clicked() {
                        should_apply = true;
                    }
                });
            });

        if should_close {
            dialog_open = false;
        }

        self.open = dialog_open;

        if should_apply {
            info!("Applying golden pipeline with config: {:?}", self.config);
            self.last_error = None;
            self.quality_report = None;
            self.quality_note = None;
            // Get selected experiment ID
            match model.get_experiment_by_selected_curve() {
                Ok(id) => {
                    let output_id = if self.config.save_to_new_experiment {
                        format!("proceeded_{id}")
                    } else {
                        id.clone()
                    };
                    match model.apply_golden_pipeline(&id, self.config.clone()) {
                        Ok(()) => {
                            self.quality_report =
                                Self::find_savitzky_golay_report(model, &output_id);
                            if self.quality_report.is_none() {
                                self.quality_note = Some(
                                "Pipeline completed, but no comparable SG quality report was recorded."
                                    .to_string(),
                            );
                            }
                            ctx.request_repaint();
                            let _ = model.push_message(&format!(
                                "Golden Pipeline completed for '{output_id}'"
                            ));
                            info!("Golden pipeline applied successfully.");
                        }
                        Err(e) => {
                            let message = format!("Golden Pipeline failed: {e:?}");
                            self.last_error = Some(message.clone());
                            let _ = model.push_message(&message);
                        }
                    }
                }
                Err(e) => {
                    let message = format!("Golden Pipeline could not start: {e:?}");
                    self.last_error = Some(message.clone());
                    let _ = model.push_message(&message);
                }
            }
        }
    }
}
