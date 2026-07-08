use crate::gui::gui_plot::PlotWindow;
use crate::gui::ivp_common::IvpTaskState;

pub struct SolidIVPApp {
    pub task: IvpTaskState,
    plot_window: Option<PlotWindow>,
    last_error: Option<String>,
}

impl Default for SolidIVPApp {
    fn default() -> Self {
        Self::new()
    }
}

impl SolidIVPApp {
    pub fn new() -> Self {
        Self {
            task: IvpTaskState::new(),
            plot_window: None,
            last_error: None,
        }
    }

    pub fn show(&mut self, ctx: &egui::Context, open: &mut bool) {
        ctx.global_style_mut(|style| {
            style.text_styles.insert(
                egui::TextStyle::Body,
                egui::FontId::new(16.0, egui::FontFamily::Proportional),
            );
            style.text_styles.insert(
                egui::TextStyle::Button,
                egui::FontId::new(16.0, egui::FontFamily::Proportional),
            );
            style.text_styles.insert(
                egui::TextStyle::Heading,
                egui::FontId::new(24.0, egui::FontFamily::Proportional),
            );
        });

        egui::Window::new("Solid state kinetic models")
            .open(open)
            .default_size([1200.0, 800.0])
            .show(ctx, |ui| {
                self.task.show_model_table(ui);
                ui.separator();
                self.task.show_solver_combo(ui);
                ui.separator();
                self.task.show_parameter_inputs(ui);
                ui.separator();
                self.task.show_problem_inputs(ui);
                ui.separator();
                self.task.show_summary(ui);
                ui.add_space(10.0);

                if ui.button("Run Model").clicked() {
                    match self.task.run_kinetic_model() {
                        Ok(result) => {
                            let plot = PlotWindow::new(
                                "time".to_string(),
                                vec![result.primary_series_name().to_string()],
                                result.time.clone(),
                                result.values.clone(),
                            );
                            self.plot_window = Some(plot);
                            self.last_error = None;
                        }
                        Err(err) => {
                            self.last_error = Some(err);
                        }
                    }
                }

                if let Some(err) = &self.last_error {
                    ui.separator();
                    ui.label(format!("Error: {}", err));
                }
            });

        if let Some(plot) = &mut self.plot_window {
            plot.show(ctx);
        }
    }
}
