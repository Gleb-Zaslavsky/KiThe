//! Column Manager Module
//!
//! This module provides a table-based UI for managing experiment columns.
//!
//! ## Architecture
//! - `ColumnManagerState`: Maintains persistent cache of column samples and UI state
//! - Checkboxes track which columns are selected for operations
//! - Operations: drop (multiple) and rename (single only)
//!
//! ## Data Flow
//! 1. `update_cache()` fetches new samples and merges into existing cache
//! 2. Old column data persists until explicitly updated
//! 3. User selects columns via checkboxes
//! 4. Buttons trigger operations on selected columns

use crate::gui::experimental_kinetics_gui::model::PlotModel;
use egui_extras::{Column, TableBuilder};
use log::info;
use std::collections::HashMap;
/// State for Column Manager window
/// Stores cached column samples and tracks user selections
#[derive(Default)]
pub struct ColumnManagerState {
    /// Persistent cache: experiment_id -> column_name -> sample_data
    /// Updated incrementally; old columns remain until replaced
    cached_samples: HashMap<String, HashMap<String, Vec<f64>>>,
    /// Tracks checked columns: (experiment_id, column_name) -> checked
    column_checkboxes: HashMap<(String, String), bool>,
    /// Input field for new column name (used in rename operation)
    new_column_name: String,
    /// System message label for errors/status
    system_message: String,
    /// Metadata editor window state
    metadata_window_open: bool,
    selected_experiment: Option<String>,
    comment_input: String,
    heating_rate_input: String,
    temperature_input: String,
    /// Confirmation dialog state for experiment dropping
    confirm_drop_experiment: Option<String>,
}

impl ColumnManagerState {
    pub fn new() -> Self {
        Self::default()
    }

    /// Updates cache with new samples from model
    /// Merges new data into existing cache, preserving old columns
    fn update_cache(&mut self, model: &mut PlotModel) -> Result<(), String> {
        let new_samples = model
            .column_samples_for_all_experiments(50)
            .map_err(|e| format!("{:?}", e))?;
        info!("new samples for experiment {:?}", &new_samples.keys());
        for (exp_id, cols) in &new_samples {
            info!("for experiment {} created columns", exp_id);
            for col in cols.keys() {
                info!("{:?}", col);
            }
        }
        for (exp_id, cols) in new_samples {
            let exp_cache = self
                .cached_samples
                .entry(exp_id)
                .or_insert_with(HashMap::new);
            for (col_name, sampled_col) in cols {
                exp_cache.insert(col_name, sampled_col.col);
            }
        }
        Ok(())
    }

    /// Renders control buttons and input field
    fn render_controls(&mut self, ui: &mut egui::Ui, model: &mut PlotModel) {
        ui.horizontal(|ui| {
            ui.label("New name:");
            ui.text_edit_singleline(&mut self.new_column_name);
        });

        ui.horizontal(|ui| {
            if ui.button("Drop Column").clicked() {
                self.handle_drop_columns(model);
            }
            if ui.button("Rename Column").clicked() {
                self.handle_rename_column(model);
            }
            if ui.button("Edit Metadata").clicked() {
                self.metadata_window_open = true;
            }
            if ui.button("Drop experiment").clicked() {
                self.handle_drop_experiment(model);
            }
        });

        if !self.system_message.is_empty() {
            ui.colored_label(egui::Color32::RED, &self.system_message);
        }
        ui.separator();
    }

    /// Handles drop operation: allows multiple selected columns
    fn handle_drop_columns(&mut self, model: &mut PlotModel) {
        let selected: Vec<_> = self
            .column_checkboxes
            .iter()
            .filter(|&(_, &checked)| checked)
            .map(|((exp_id, col_name), _)| (exp_id.clone(), col_name.clone()))
            .collect();

        if selected.is_empty() {
            self.system_message = "No columns selected for dropping".to_string();
            return;
        }

        for (exp_id, col_name) in selected {
            if let Err(e) = model.drop_column_for_experiment(&exp_id, &col_name) {
                self.system_message = format!("Error dropping {}: {:?}", col_name, e);
                return;
            }
            self.cached_samples
                .get_mut(&exp_id)
                .and_then(|cols| cols.remove(&col_name));
            self.column_checkboxes.remove(&(exp_id, col_name));
        }
        self.system_message = "Columns dropped successfully".to_string();
    }

    /// Handles rename operation: requires exactly one selected column
    fn handle_rename_column(&mut self, model: &mut PlotModel) {
        let selected: Vec<_> = self
            .column_checkboxes
            .iter()
            .filter(|&(_, &checked)| checked)
            .map(|((exp_id, col_name), _)| (exp_id.clone(), col_name.clone()))
            .collect();

        if selected.len() != 1 {
            self.system_message = "Select exactly one column for renaming".to_string();
            return;
        }

        if self.new_column_name.trim().is_empty() {
            self.system_message = "Enter new column name".to_string();
            return;
        }

        let (exp_id, old_name) = &selected[0];
        if let Err(e) =
            model.rename_column_for_experiment(exp_id, old_name, self.new_column_name.trim())
        {
            self.system_message = format!("Error renaming: {:?}", e);
            return;
        }

        // Update cache with new name
        if let Some(exp_cols) = self.cached_samples.get_mut(exp_id) {
            if let Some(data) = exp_cols.remove(old_name) {
                exp_cols.insert(self.new_column_name.trim().to_string(), data);
            }
        }
        self.column_checkboxes
            .remove(&(exp_id.clone(), old_name.clone()));
        self.system_message = "Column renamed successfully".to_string();
        self.new_column_name.clear();
    }

    /// Handles drop experiment operation with confirmation dialog
    fn handle_drop_experiment(&mut self, model: &mut PlotModel) {
        // First, find which experiments have selected columns
        let selected_experiments: Vec<_> = self
            .column_checkboxes
            .iter()
            .filter(|&(_, &checked)| checked)
            .map(|((exp_id, _), _)| exp_id.clone())
            .collect();

        if selected_experiments.is_empty() {
            self.system_message = "No columns selected for experiment dropping".to_string();
            return;
        }

        // For simplicity, we'll use the first selected experiment
        // In a more complex UI, we might show a dialog to select which experiment to drop
        let exp_id = &selected_experiments[0];

        // Show confirmation dialog
        self.confirm_drop_experiment = Some(exp_id.clone());
    }

    /// Renders the data table with checkboxes in header
    /// Each column represents one experiment-column pair
    fn render_table(&mut self, ui: &mut egui::Ui, model: &PlotModel) {
        let experiments = model.list_of_experiments();
        if experiments.is_empty() {
            ui.label("No experiments loaded");
            return;
        }

        // Build list of all (exp_id, col_name) pairs for table columns
        let mut all_columns = Vec::new();
        for exp_id in &experiments {
            if let Ok(cols) = model.list_of_columns(exp_id) {
                for col_name in cols {
                    all_columns.push((exp_id.clone(), col_name));
                }
            }
        }

        egui::ScrollArea::both().show(ui, |ui| {
            TableBuilder::new(ui)
                .striped(true)
                .column(Column::auto().resizable(true))
                .columns(Column::auto().resizable(true), all_columns.len())
                .header(40.0, |mut header| {
                    header.col(|ui| {
                        ui.heading("Row");
                    });
                    // Each column header: checkbox + "exp_id col_name"
                    for (exp_id, col_name) in &all_columns {
                        header.col(|ui| {
                            ui.vertical(|ui| {
                                let key = (exp_id.clone(), col_name.clone());
                                let checked = self.column_checkboxes.entry(key).or_insert(false);
                                ui.checkbox(checked, "");
                                ui.label(format!("{} {}", exp_id, col_name));
                            });
                        });
                    }
                })
                .body(|mut body| {
                    let max_rows = self
                        .cached_samples
                        .values()
                        .flat_map(|cols| cols.values().map(|v| v.len()))
                        .max()
                        .unwrap_or(0);

                    for row_idx in 0..max_rows {
                        body.row(18.0, |mut row| {
                            row.col(|ui| {
                                ui.label(row_idx.to_string());
                            });
                            for (exp_id, col_name) in &all_columns {
                                row.col(|ui| {
                                    if let Some(exp_cols) = self.cached_samples.get(exp_id) {
                                        if let Some(col_data) = exp_cols.get(col_name) {
                                            if let Some(val) = col_data.get(row_idx) {
                                                ui.label(format!("{:.4}", val));
                                            } else {
                                                ui.label("-");
                                            }
                                        } else {
                                            ui.label("-");
                                        }
                                    } else {
                                        ui.label("-");
                                    }
                                });
                            }
                        });
                    }
                });
        });
    }

    /// Renders metadata editor window
    fn show_metadata_window(&mut self, ctx: &egui::Context, model: &mut PlotModel) {
        if !self.metadata_window_open {
            return;
        }

        egui::Window::new("Edit Metadata")
            .open(&mut self.metadata_window_open.clone())
            .resizable(true)
            .default_size([400.0, 300.0])
            .show(ctx, |ui| {
                let experiments = model.list_of_experiments();

                ui.label("Select Experiment:");
                egui::ComboBox::from_id_salt("exp_selector")
                    .selected_text(self.selected_experiment.as_deref().unwrap_or("None"))
                    .show_ui(ui, |ui| {
                        for exp_id in &experiments {
                            ui.selectable_value(
                                &mut self.selected_experiment,
                                Some(exp_id.clone()),
                                exp_id,
                            );
                        }
                    });

                ui.separator();

                ui.label("Comment:");
                ui.text_edit_multiline(&mut self.comment_input);

                ui.label("Heating Rate (K/min):");
                ui.text_edit_singleline(&mut self.heating_rate_input);

                ui.label("Temperature (K):");
                ui.text_edit_singleline(&mut self.temperature_input);

                ui.separator();

                if ui.button("Apply").clicked() {
                    if let Some(exp_id) = &self.selected_experiment {
                        let mut success = true;

                        if !self.comment_input.trim().is_empty() {
                            if let Err(e) = model.set_comment(exp_id, self.comment_input.trim()) {
                                self.system_message = format!("Error setting comment: {:?}", e);
                                success = false;
                            }
                        }

                        if !self.heating_rate_input.trim().is_empty() {
                            if let Ok(rate) = self.heating_rate_input.trim().parse::<f64>() {
                                if let Err(e) = model.set_heating_rate(exp_id, rate) {
                                    self.system_message =
                                        format!("Error setting heating rate: {:?}", e);
                                    success = false;
                                }
                            } else {
                                self.system_message = "Heating rate must be a number".to_string();
                                success = false;
                            }
                        }

                        if !self.temperature_input.trim().is_empty() {
                            if let Ok(temp) = self.temperature_input.trim().parse::<f64>() {
                                if let Err(e) = model.set_experiment_temperature(exp_id, temp) {
                                    self.system_message =
                                        format!("Error setting temperature: {:?}", e);
                                    success = false;
                                }
                            } else {
                                self.system_message = "Temperature must be a number".to_string();
                                success = false;
                            }
                        }

                        if success {
                            self.system_message = "Metadata updated successfully".to_string();
                            self.metadata_window_open = false;
                        }
                    } else {
                        self.system_message = "Select an experiment first".to_string();
                    }
                }
            });
    }
}

/// Main window function for Column Manager
pub fn show_column_manager_window(
    ctx: &egui::Context,
    model: &mut PlotModel,
    state: &mut ColumnManagerState,
    open: &mut bool,
) {
    state.show_metadata_window(ctx, model);
    state.system_message.clear();
    egui::Window::new("Column Manager")
        .open(open)
        .resizable(true)
        .default_size([800.0, 600.0])
        .show(ctx, |ui| {
            if let Err(e) = state.update_cache(model) {
                ui.colored_label(egui::Color32::RED, format!("Error: {}", e));
                return;
            }

            // Handle experiment drop confirmation dialog
            if let Some(ref exp_id) = state.confirm_drop_experiment {
                let mut confirmed = false;
                let mut cancelled = false;

                egui::Window::new("Confirm Drop Experiment")
                    .collapsible(false)
                    .resizable(false)
                    .show(ctx, |ui| {
                        ui.label(format!(
                            "Are you sure you want to drop experiment {}?",
                            exp_id
                        ));
                        ui.horizontal(|ui| {
                            if ui.button("Yes").clicked() {
                                confirmed = true;
                            }
                            if ui.button("No").clicked() {
                                cancelled = true;
                            }
                        });
                    });

                if confirmed {
                    if let Err(e) = model.drop_experiment_by_id(exp_id) {
                        state.system_message = format!("Error dropping experiment: {:?}", e);
                    } else {
                        state.system_message =
                            format!("Experiment {} dropped successfully", exp_id);
                        // Clear cached data for this experiment
                        state.cached_samples.remove(exp_id);
                        // Remove all checkboxes for this experiment
                        state.column_checkboxes.retain(|(eid, _), _| eid != exp_id);
                    }
                    state.confirm_drop_experiment = None;
                }

                if cancelled {
                    state.confirm_drop_experiment = None;
                }
            }

            state.render_controls(ui, model);
            state.render_table(ui, model);
        });
}
