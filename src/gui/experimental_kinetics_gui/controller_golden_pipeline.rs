//! Golden Pipeline dialog for configuring and applying the golden pipeline transformation.

use crate::Kinetics::experimental_kinetics::exp_engine_api::GoldenPipelineConfig;
use crate::gui::experimental_kinetics_gui::model::PlotModel;
use eframe::egui;
use log::info;

/// State for the golden pipeline dialog.
#[derive(Clone, Debug)]
pub struct GoldenPipelineDialogState {
    pub open: bool,
    pub config: GoldenPipelineConfig,
    pub last_error: Option<String>,
}

impl Default for GoldenPipelineDialogState {
    fn default() -> Self {
        Self {
            open: false,
            config: GoldenPipelineConfig::default(),
            last_error: None,
        }
    }
}

impl GoldenPipelineDialogState {
    /// Open the dialog with default config.
    pub fn open(&mut self) {
        self.open = true;
        self.config = GoldenPipelineConfig::default();
        self.last_error = None;
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

                ui.horizontal(|ui| {
                    if ui.button("Cancel").clicked() {
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
            // Get selected experiment ID
            match model.get_experiment_by_selected_curve() {
                Ok(id) => match model.apply_golden_pipeline(&id, self.config.clone()) {
                    Ok(()) => {
                        self.last_error = None;
                        self.open = false;
                        info!("Golden pipeline applied successfully.");
                    }
                    Err(e) => {
                        self.last_error = Some(format!("Failed to apply golden pipeline: {:?}", e));
                    }
                },
                Err(e) => {
                    self.last_error = Some(format!("No selected experiment: {:?}", e));
                }
            }
        }
    }
}
