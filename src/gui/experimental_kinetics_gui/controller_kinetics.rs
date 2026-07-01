use crate::gui::experimental_kinetics_gui::model::{Colours, PlotModel};
use crate::gui::ivp_common::{IvpTaskState, ivp_result_to_experiment};

#[derive(Debug, Clone)]
pub struct DirectProblemDialogState {
    pub open: bool,
    pub task: IvpTaskState,
    pub output_experiment_id: String,
    pub last_error: Option<String>,
    pub last_success: Option<String>,
}

impl Default for DirectProblemDialogState {
    fn default() -> Self {
        Self::new()
    }
}

impl DirectProblemDialogState {
    pub fn new() -> Self {
        Self {
            open: false,
            task: IvpTaskState::new(),
            output_experiment_id: String::from("ivp_result"),
            last_error: None,
            last_success: None,
        }
    }

    pub fn open(&mut self) {
        self.open = true;
        self.last_error = None;
        self.last_success = None;
    }

    fn resolved_experiment_id(&self, model: &PlotModel) -> String {
        let base = self.output_experiment_id.trim();
        let candidate = if base.is_empty() {
            "ivp_result".to_string()
        } else {
            base.to_string()
        };
        if !model.series.ids().iter().any(|id| id == &candidate) {
            return candidate;
        }

        let mut suffix = 1usize;
        loop {
            let next = format!("{}_{}", candidate, suffix);
            if !model.series.ids().iter().any(|id| id == &next) {
                return next;
            }
            suffix += 1;
        }
    }

    pub fn show(&mut self, ctx: &egui::Context, model: &mut PlotModel) {
        if !self.open {
            return;
        }

        let mut open = self.open;
        egui::Window::new("Direct Problem")
            .open(&mut open)
            .default_size([1100.0, 760.0])
            .show(ctx, |ui| {
                self.task.show_model_table(ui);
                ui.separator();
                self.task.show_solver_combo(ui);
                ui.separator();
                self.task.show_parameter_inputs(ui);
                ui.separator();
                self.task.show_problem_inputs(ui);
                ui.separator();

                ui.horizontal(|ui| {
                    ui.label("output experiment:");
                    ui.text_edit_singleline(&mut self.output_experiment_id);
                });

                self.task.show_summary(ui);
                ui.add_space(10.0);

                if ui.button("Solve IVP").clicked() {
                    match self.task.run_kinetic_model() {
                        Ok(result) => {
                            let experiment_id = self.resolved_experiment_id(model);
                            match ivp_result_to_experiment(&result, experiment_id.clone()) {
                                Ok(experiment) => match model
                                    .push_experiment_and_create_plot(experiment, Colours::Red)
                                {
                                    Ok(()) => {
                                        self.last_error = None;
                                        self.last_success =
                                            Some(format!("Created experiment '{}'", experiment_id));
                                    }
                                    Err(err) => {
                                        self.last_error = Some(format!("{:?}", err));
                                    }
                                },
                                Err(err) => {
                                    self.last_error = Some(format!("{:?}", err));
                                }
                            }
                        }
                        Err(err) => {
                            self.last_error = Some(err);
                        }
                    }
                }

                if let Some(msg) = &self.last_success {
                    ui.separator();
                    ui.label(msg);
                }
                if let Some(err) = &self.last_error {
                    ui.separator();
                    ui.label(format!("Error: {}", err));
                }
            });
        self.open = open;
    }
}
