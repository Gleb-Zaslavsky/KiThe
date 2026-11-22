use crate::gui::experimental_kinetics_gui::model::{Colours, PlotModel};
use crate::gui::ivp_common::{IvpTaskState, ivp_result_to_experiment};
use crate::gui::network_ivp_gui::NetworkIvpApp;

pub struct DirectProblemDialogState {
    pub open: bool,
    pub task: IvpTaskState,
    pub output_experiment_id: String,
    pub last_error: Option<String>,
    pub last_success: Option<String>,
    pub network_open: bool,
    pub network_app: NetworkIvpApp,
}

impl Default for DirectProblemDialogState {
    fn default() -> Self {
        Self::new()
    }
}

impl DirectProblemDialogState {
    pub(crate) fn close_window_on_layer(&mut self, layer: Option<egui::LayerId>) -> bool {
        if self.open
            && layer.is_some_and(|layer| layer.id == egui::Id::unique(Some("Direct Problem")))
        {
            self.open = false;
            return true;
        }
        false
    }

    pub fn new() -> Self {
        Self {
            open: false,
            task: IvpTaskState::new(),
            output_experiment_id: String::from("ivp_result"),
            last_error: None,
            last_success: None,
            network_open: false,
            network_app: NetworkIvpApp::new_for_experimental_kinetics(),
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
        if self.open {
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
                    let label = ui.label("Output experiment:");
                    let input = ui.text_edit_singleline(&mut self.output_experiment_id);
                    let _ = input.labelled_by(label.id);
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
                                        let message = format!("Created experiment '{experiment_id}'");
                                        self.last_success = Some(message.clone());
                                        let _ = model.push_message(&message);
                                    }
                                    Err(crate::gui::experimental_kinetics_gui::model::TGAGUIError::PostCommitWarning(warning)) => {
                                        self.last_error = None;
                                        let message = format!("Warning: {warning}");
                                        self.last_success = Some(message.clone());
                                        let _ = model.push_message(&message);
                                    }
                                    Err(err) => {
                                        let message = format!("Error: {err}");
                                        self.last_error = Some(message.clone());
                                        let _ = model.push_message(&message);
                                    }
                                },
                                Err(err) => {
                                    let message = format!("Error: {err}");
                                    self.last_error = Some(message.clone());
                                    let _ = model.push_message(&message);
                                }
                            }
                        }
                        Err(err) => {
                            let message = format!("Error: {err}");
                            self.last_error = Some(message.clone());
                            let _ = model.push_message(&message);
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
        self.network_app
            .set_existing_experiment_ids(model.list_of_experiments());
        self.network_app.show(ctx, &mut self.network_open);
        if let Some((prefix, network_description, solutions)) =
            self.network_app.take_publish_request()
        {
            let mut runs = Vec::with_capacity(solutions.len());
            let mut seen = std::collections::HashSet::new();
            let existing: std::collections::HashSet<_> =
                model.list_of_experiments().into_iter().collect();
            let mut id_error = None;
            for solution in solutions {
                let id = solution.run_metadata.temperature_binding.map_or_else(
                    || {
                        format!(
                            "{}_symbolic",
                            if prefix.is_empty() {
                                "network"
                            } else {
                                &prefix
                            }
                        )
                    },
                    |binding| {
                        crate::gui::network_ivp_gui::network_condition_experiment_id(
                            &prefix, binding,
                        )
                    },
                );
                if existing.contains(&id) || !seen.insert(id.clone()) {
                    id_error = Some(format!(
                        "Experiment ID '{id}' already exists or repeats in this batch"
                    ));
                    break;
                }
                let y_columns = crate::gui::ivp_common::network_solution_species_columns(&solution);
                match crate::gui::ivp_common::network_solution_to_experiment(
                    &solution,
                    id.clone(),
                    &network_description,
                ) {
                    Ok(experiment) => runs.push((experiment, y_columns)),
                    Err(error) => {
                        id_error = Some(error.to_string());
                        break;
                    }
                }
            }
            if let Some(error) = id_error {
                let _ =
                    model.push_message(&format!("Reaction-network series was not added: {error}"));
            } else {
                match model.push_network_series_and_create_plots(runs) {
                    Ok(()) => {
                    let _ = model.push_message(
                        "Reaction-network condition series added to the experiment list",
                    );
                    }
                    Err(crate::gui::experimental_kinetics_gui::model::TGAGUIError::PostCommitWarning(warning)) => {
                        let _ = model.push_message(&format!("Reaction-network series added; plot warning: {warning}"));
                    }
                    Err(error) => {
                        let _ = model.push_message(&format!("Reaction-network series was not added: {error}"));
                    }
                }
            }
        }
    }
}
