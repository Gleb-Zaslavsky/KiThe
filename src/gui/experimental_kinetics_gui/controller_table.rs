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

use crate::Kinetics::experimental_kinetics::experiment_series_main::{
    ExperimentThermalMode, TGASeries,
};
use crate::Kinetics::experimental_kinetics::one_experiment_dataset::PreviewRows;
use crate::gui::experimental_kinetics_gui::controller_methods::KineticMethodsWindowState;
use crate::gui::experimental_kinetics_gui::model::PlotModel;
use egui_extras::{Column, TableBuilder};
use log::info;
use std::collections::{HashMap, HashSet};

fn format_sample_value(value: f64) -> String {
    if value.is_finite() && value != 0.0 && (value.abs() < 0.001 || value.abs() >= 1_000_000.0) {
        format!("{value:.4e}")
    } else {
        format!("{value:.4}")
    }
}

fn display_experiment_id(id: &str) -> &str {
    if id.is_empty() { "(unnamed)" } else { id }
}

fn data_column_width(available_width: f32, column_count: usize) -> f32 {
    if column_count == 0 {
        return 120.0;
    }
    (available_width / column_count as f32).clamp(88.0, 220.0)
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum MessageKind {
    Info,
    Success,
    Warning,
    Error,
}

#[derive(Clone, Debug)]
struct SystemMessage {
    kind: MessageKind,
    text: String,
}

impl SystemMessage {
    fn new(kind: MessageKind, text: impl Into<String>) -> Self {
        Self {
            kind,
            text: text.into(),
        }
    }

    fn accessible_text(&self) -> String {
        let level = match self.kind {
            MessageKind::Info => "Info",
            MessageKind::Success => "Success",
            MessageKind::Warning => "Warning",
            MessageKind::Error => "Error",
        };
        format!("{level}: {}", self.text)
    }
}

#[cfg(test)]
mod sample_format_tests {
    use super::{MessageKind, SystemMessage, data_column_width, format_sample_value};

    #[test]
    fn table_column_width_adapts_to_number_of_columns() {
        let few_columns = data_column_width(800.0, 2);
        let many_columns = data_column_width(800.0, 12);

        assert!(few_columns > many_columns);
        assert_eq!(many_columns, 88.0);
    }

    #[test]
    fn small_nonzero_values_remain_visible() {
        assert_eq!(format_sample_value(0.000_012_3), "1.2300e-5");
        assert_eq!(format_sample_value(-0.000_012_3), "-1.2300e-5");
    }

    #[test]
    fn ordinary_values_keep_compact_decimal_format() {
        assert_eq!(format_sample_value(12.345_67), "12.3457");
        assert_eq!(format_sample_value(0.0), "0.0000");
    }

    #[test]
    fn system_messages_have_distinct_severity() {
        assert_eq!(
            SystemMessage::new(MessageKind::Success, "Column renamed").accessible_text(),
            "Success: Column renamed"
        );
        assert_eq!(
            SystemMessage::new(MessageKind::Warning, "Time already uses hours").kind,
            MessageKind::Warning
        );
        assert_eq!(
            SystemMessage::new(MessageKind::Error, "Could not rename").kind,
            MessageKind::Error
        );
        assert_eq!(
            SystemMessage::new(MessageKind::Info, "Ready").accessible_text(),
            "Info: Ready"
        );
    }

    #[test]
    fn column_manager_status_is_published_to_the_main_message_label() {
        let mut state = super::ColumnManagerState::default();
        let mut model = crate::gui::experimental_kinetics_gui::model::PlotModel::new();
        state.set_message(MessageKind::Error, "Could not drop column: busy");

        state.publish_message_to_main_status(&mut model);

        assert_eq!(model.message, "Error: Could not drop column: busy");
        state.set_message(MessageKind::Success, "Column dropped");
        state.publish_message_to_main_status(&mut model);
        assert_eq!(model.message, "Success: Column dropped");

        model.message = "Saved another file".to_string();
        state.set_message(MessageKind::Error, "Could not drop column: busy");
        state.publish_message_to_main_status(&mut model);
        assert_eq!(model.message, "Error: Could not drop column: busy");

        state.set_message_if_changed(MessageKind::Error, "Column sample refresh failed");
        let revision = state.message_revision;
        state.set_message_if_changed(MessageKind::Error, "Column sample refresh failed");
        assert_eq!(state.message_revision, revision);
    }

    #[test]
    fn replacing_document_clears_all_pending_delete_confirmations() {
        let mut state = super::ColumnManagerState::default();
        state.confirm_drop_experiment = Some("same_id".to_string());
        state.confirm_drop_columns = vec![("same_id".to_string(), "mass".to_string())];
        state.confirm_delete_all = true;

        state.clear_pending_confirmations();

        assert!(state.confirm_drop_experiment.is_none());
        assert!(state.confirm_drop_columns.is_empty());
        assert!(!state.confirm_delete_all);
    }
}
/// State for Column Manager window
/// Stores cached column samples and tracks user selections
#[derive(Default)]
pub struct ColumnManagerState {
    /// Persistent cache: experiment_id -> column_name -> sample_data
    /// Updated incrementally; old columns remain until replaced
    cached_samples: HashMap<String, PreviewRows>,
    /// Tracks checked columns: (experiment_id, column_name) -> checked
    column_checkboxes: HashMap<(String, String), bool>,
    /// Input field for new column name (used in rename operation)
    new_column_name: String,
    /// System message label for errors/status
    system_message: Option<SystemMessage>,
    message_revision: u64,
    last_published_revision: u64,
    /// Metadata editor window state
    metadata_window_open: bool,
    selected_experiment: Option<String>,
    table_experiment: Option<String>,
    compare_mode: bool,
    compare_experiments: HashSet<String>,
    selected_column: Option<(String, String)>,
    comment_input: String,
    heating_rate_input: String,
    temperature_input: String,
    /// Confirmation dialog state for experiment dropping
    confirm_drop_experiment: Option<String>,
    confirm_drop_columns: Vec<(String, String)>,
    confirm_delete_all: bool,
    details_tab: DetailsTab,
}

#[derive(Default, Clone, Copy, PartialEq, Eq)]
enum DetailsTab {
    #[default]
    Overview,
    Processing,
    Quality,
    History,
}

impl ColumnManagerState {
    pub(crate) fn clear_pending_confirmations(&mut self) {
        self.confirm_drop_experiment = None;
        self.confirm_drop_columns.clear();
        self.confirm_delete_all = false;
    }

    fn set_message(&mut self, kind: MessageKind, text: impl Into<String>) {
        self.system_message = Some(SystemMessage::new(kind, text));
        self.message_revision = self.message_revision.wrapping_add(1);
    }

    fn set_message_if_changed(&mut self, kind: MessageKind, text: impl Into<String>) {
        let text = text.into();
        if self
            .system_message
            .as_ref()
            .is_some_and(|current| current.kind == kind && current.text == text)
        {
            return;
        }
        self.set_message(kind, text);
    }

    fn publish_message_to_main_status(&mut self, model: &mut PlotModel) {
        if self.message_revision == self.last_published_revision {
            return;
        }
        let Some(message) = &self.system_message else {
            return;
        };
        let text = message.accessible_text();
        model.message = text;
        self.last_published_revision = self.message_revision;
    }

    pub(crate) fn close_window_on_layer(&mut self, layer: Option<egui::LayerId>) -> bool {
        let Some(layer) = layer else { return false };
        if layer.id == egui::Id::unique(Some("Confirm Drop Experiment"))
            && self.confirm_drop_experiment.is_some()
        {
            self.confirm_drop_experiment = None;
            return true;
        }
        if layer.id == egui::Id::unique(Some("Confirm Drop Columns"))
            && !self.confirm_drop_columns.is_empty()
        {
            self.confirm_drop_columns.clear();
            return true;
        }
        if layer.id == egui::Id::unique(Some("Confirm Delete All Workspace"))
            && self.confirm_delete_all
        {
            self.confirm_delete_all = false;
            return true;
        }
        if layer.id == egui::Id::unique(Some("Edit Metadata")) && self.metadata_window_open {
            self.metadata_window_open = false;
            return true;
        }
        false
    }

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
        for (exp_id, preview) in &new_samples {
            info!("for experiment {} created columns", exp_id);
            for col in preview.columns.keys() {
                info!("{:?}", col);
            }
        }
        for (exp_id, preview) in new_samples {
            self.cached_samples.insert(exp_id, preview);
        }

        let experiments = model.list_of_experiments();
        self.cached_samples
            .retain(|exp_id, _| experiments.contains(exp_id));

        for exp_id in experiments {
            let columns = model
                .list_of_columns(&exp_id)
                .map_err(|e| format!("{:?}", e))?;
            if let Some(preview) = self.cached_samples.get_mut(&exp_id) {
                preview.columns.retain(|name, _| columns.contains(name));
            }
        }
        Ok(())
    }

    /// Renders control buttons and input field
    fn render_controls(&mut self, ui: &mut egui::Ui, model: &mut PlotModel) {
        let experiments = model.list_of_experiments();
        if self
            .table_experiment
            .as_ref()
            .is_none_or(|id| !experiments.contains(id))
        {
            self.table_experiment = experiments.first().cloned();
        }
        self.compare_experiments
            .retain(|id| experiments.contains(id));
        ui.horizontal(|ui| {
            ui.checkbox(&mut self.compare_mode, "Compare experiments")
                .on_hover_text("Read-only view of equally sized samples from selected experiments. Source grids are not interpolated.");
            if !self.compare_mode {
                egui::ComboBox::from_label("Experiment")
                    .selected_text(
                        self.table_experiment
                            .as_deref()
                            .unwrap_or("Select experiment"),
                    )
                    .show_ui(ui, |ui| {
                        for id in &experiments {
                            ui.selectable_value(&mut self.table_experiment, Some(id.clone()), id);
                        }
                    });
            }
        });
        if self.compare_mode {
            egui::CollapsingHeader::new("Experiments to compare")
                .default_open(true)
                .show(ui, |ui| {
                    ui.horizontal_wrapped(|ui| {
                        for id in &experiments {
                            let mut selected = self.compare_experiments.contains(id);
                            let label = if id.is_empty() {
                                "(unnamed experiment)"
                            } else {
                                id
                            };
                            if ui.checkbox(&mut selected, label).changed() {
                                if selected {
                                    self.compare_experiments.insert(id.clone());
                                } else {
                                    self.compare_experiments.remove(id);
                                }
                            }
                        }
                    });
                });
        }
        let mutation_enabled = !self.compare_mode && self.table_experiment.is_some();
        ui.horizontal(|ui| {
            let selected_count = self.selected_column_checkboxes().len();
            if ui
                .add_enabled(
                    mutation_enabled && selected_count > 0,
                    egui::Button::new(format!("Delete columns ({selected_count})")),
                )
                .on_hover_text(
                    "Deletes the checked columns from the current experiment after confirmation.",
                )
                .clicked()
            {
                self.handle_drop_columns();
            }
            if let Some(id) = self.table_experiment.clone() {
                if ui
                    .add_enabled(
                        mutation_enabled,
                        egui::Button::new(format!(
                            "Delete experiment: {}",
                            display_experiment_id(&id)
                        )),
                    )
                    .on_hover_text(
                        "Deletes the current experiment and its plots after confirmation.",
                    )
                    .clicked()
                {
                    self.confirm_drop_experiment = Some(id);
                }
            }
            ui.menu_button("More", |ui| {
                if ui
                    .add_enabled(mutation_enabled, egui::Button::new("Edit Metadata"))
                    .on_hover_text("Edit the current experiment's comment and thermal metadata. This does not change columns.")
                    .clicked()
                {
                    self.selected_experiment = self.table_experiment.clone();
                    self.metadata_window_open = true;
                    ui.close();
                }
                ui.separator();
                if ui
                    .add_enabled(
                        mutation_enabled,
                        egui::Button::new("Check time monotonicity"),
                    )
                    .on_hover_text("Checks the current experiment's bound time column. Enter no values; data is not modified.")
                    .clicked()
                {
                    self.handle_check_time_monotonity(model);
                    ui.close();
                }
                if ui
                    .add_enabled(mutation_enabled, egui::Button::new("Sort by time"))
                    .on_hover_text("Sorts every row in the current experiment by its bound time column. No new column is created.")
                    .clicked()
                {
                    self.handle_sort_by_time(model);
                    ui.close();
                }
            });
            ui.with_layout(
                egui::Layout::right_to_left(egui::Align::Center),
                |ui| {
                    if ui
                        .button("Delete All")
                        .on_hover_text("Clear every experiment, plot, and selection from the workspace after an explicit confirmation.")
                        .clicked()
                    {
                        self.confirm_delete_all = true;
                    }
                },
            );
        });

        if let Some(message) = &self.system_message {
            let color = match message.kind {
                MessageKind::Info => ui.visuals().text_color(),
                MessageKind::Success => egui::Color32::LIGHT_GREEN,
                MessageKind::Warning => egui::Color32::YELLOW,
                MessageKind::Error => egui::Color32::LIGHT_RED,
            };
            ui.colored_label(color, message.accessible_text());
        }
        ui.separator();
    }

    /// Handles drop operation: allows multiple selected columns
    fn selected_column_checkboxes(&self) -> Vec<(String, String)> {
        let Some(active_id) = self.table_experiment.as_deref() else {
            return Vec::new();
        };
        let mut selected = self
            .column_checkboxes
            .iter()
            .filter(|((exp_id, _), checked)| exp_id == active_id && **checked)
            .map(|((exp_id, name), _)| (exp_id.clone(), name.clone()))
            .collect::<Vec<_>>();
        selected.sort();
        selected
    }

    fn handle_drop_columns(&mut self) {
        if self.compare_mode {
            self.set_message(
                MessageKind::Error,
                "Column deletion is unavailable in comparison mode",
            );
            return;
        }
        let selected = self.selected_column_checkboxes();
        if selected.is_empty() {
            self.set_message(MessageKind::Error, "No columns selected for dropping");
            return;
        }
        self.confirm_drop_columns = selected;
    }

    fn commit_drop_columns(&mut self, model: &mut PlotModel) {
        let selected = std::mem::take(&mut self.confirm_drop_columns);
        let Some((exp_id, _)) = selected.first() else {
            self.set_message(MessageKind::Error, "No columns selected for dropping");
            return;
        };
        if selected
            .iter()
            .any(|(selected_exp, _)| selected_exp != exp_id)
        {
            self.set_message(
                MessageKind::Error,
                "A column-delete batch must target one experiment",
            );
            return;
        }
        let col_names = selected
            .iter()
            .map(|(_, name)| name.clone())
            .collect::<Vec<_>>();
        if let Err(error) = model.drop_columns_for_experiment(exp_id, &col_names) {
            self.set_message(
                MessageKind::Error,
                format!("Could not drop columns: {error:?}"),
            );
            return;
        }
        if let Some(preview) = self.cached_samples.get_mut(exp_id) {
            for name in &col_names {
                preview.columns.remove(name);
            }
        }
        for name in col_names {
            self.column_checkboxes.remove(&(exp_id.clone(), name));
        }
        self.set_message(MessageKind::Success, "Columns dropped successfully");
    }

    /// Handles rename operation: requires exactly one selected column
    fn handle_rename_column(&mut self, model: &mut PlotModel, exp_id: &str, old_name: &str) {
        if self.compare_mode {
            self.set_message(
                MessageKind::Error,
                "Column renaming is unavailable in comparison mode",
            );
            return;
        }

        if self.new_column_name.trim().is_empty() {
            self.set_message(MessageKind::Error, "Enter new column name");
            return;
        }

        if let Err(e) =
            model.rename_column_for_experiment(exp_id, old_name, self.new_column_name.trim())
        {
            self.set_message(
                MessageKind::Error,
                format!("Could not rename column: {e:?}"),
            );
            return;
        }

        if let Some(preview) = self.cached_samples.get_mut(exp_id) {
            preview.columns.remove(old_name);
        }
        self.column_checkboxes
            .remove(&(exp_id.to_string(), old_name.to_string()));
        if self
            .selected_column
            .as_ref()
            .map(|(selected_exp, selected_col)| selected_exp == exp_id && selected_col == old_name)
            .unwrap_or(false)
        {
            self.selected_column =
                Some((exp_id.to_string(), self.new_column_name.trim().to_string()));
        }
        self.set_message(MessageKind::Success, "Column renamed successfully");
        self.new_column_name.clear();
    }

    /// Handles drop experiment operation with confirmation dialog
    fn handle_check_time_monotonity(&mut self, model: &mut PlotModel) {
        let selected_experiments = self.table_experiment.iter().cloned().collect::<Vec<_>>();

        if selected_experiments.is_empty() {
            self.set_message(MessageKind::Error, "Select an experiment first");
            info!("Select an experiment first");
            return;
        }

        let mut per_experiment_messages = Vec::new();

        for exp_id in selected_experiments {
            match model.time_monotonicity_report_for_experiment(&exp_id) {
                Ok(report) if report.is_strictly_increasing => {
                    let msg = format!(
                        "Experiment '{}': time monotonicity OK for '{}'",
                        exp_id, report.time_column
                    );
                    info!("{}", msg);
                    per_experiment_messages.push(msg);
                }
                Ok(report) => {
                    let mut details = Vec::new();
                    if !report.null_indices.is_empty() {
                        details.push(format!("nulls at {:?}", report.null_indices));
                    }
                    if !report.nan_indices.is_empty() {
                        details.push(format!("NaNs at {:?}", report.nan_indices));
                    }
                    if !report.infinite_indices.is_empty() {
                        details.push(format!("infinities at {:?}", report.infinite_indices));
                    }
                    if !report.duplicate_pairs.is_empty() {
                        details.push(format!("duplicates at {:?}", report.duplicate_pairs));
                    }
                    if !report.decreasing_pairs.is_empty() {
                        details.push(format!("decreasing at {:?}", report.decreasing_pairs));
                    }
                    let msg = format!(
                        "Experiment '{}': time monotonicity FAILED for '{}' ({})",
                        exp_id,
                        report.time_column,
                        details.join("; ")
                    );
                    info!("{}", msg);
                    per_experiment_messages.push(msg);
                }
                Err(e) => {
                    let msg = format!(
                        "Experiment '{}': time monotonicity check error: {}",
                        exp_id, e
                    );
                    info!("{}", msg);
                    per_experiment_messages.push(msg);
                }
            }
        }

        let kind = if per_experiment_messages
            .iter()
            .any(|message| message.contains("FAILED") || message.contains("error:"))
        {
            MessageKind::Warning
        } else {
            MessageKind::Success
        };
        self.set_message(kind, per_experiment_messages.join(" | "));
    }

    fn handle_sort_by_time(&mut self, model: &mut PlotModel) {
        let selected_experiments = self.table_experiment.iter().cloned().collect::<Vec<_>>();

        if selected_experiments.is_empty() {
            self.set_message(MessageKind::Error, "Select an experiment first");
            info!("Select an experiment first");
            return;
        }

        let mut per_experiment_messages = Vec::new();

        for exp_id in selected_experiments {
            match model.time_monotonicity_report_for_experiment(&exp_id) {
                Ok(report) if !report.sortable_without_data_loss => {
                    let msg = format!(
                        "Experiment '{}': cannot sort time column '{}' because it contains null/NaN/inf values",
                        exp_id, report.time_column
                    );
                    info!("{}", msg);
                    per_experiment_messages.push(msg);
                }
                Ok(report) if report.is_strictly_increasing => {
                    let msg = format!("Experiment '{}': time already monotonic", exp_id);
                    info!("{}", msg);
                    per_experiment_messages.push(msg);
                }
                Ok(report) => match model.sort_experiment_by_time(&exp_id) {
                    Ok(()) => {
                        let msg = if !report.duplicate_pairs.is_empty() {
                            format!(
                                "Experiment '{}': sorted by time; duplicate time values remain in '{}'",
                                exp_id, report.time_column
                            )
                        } else {
                            format!(
                                "Experiment '{}': sorted by time column '{}'",
                                exp_id, report.time_column
                            )
                        };
                        info!("{}", msg);
                        per_experiment_messages.push(msg);
                    }
                    Err(e) => {
                        let msg = format!("Experiment '{}': sort by time failed: {}", exp_id, e);
                        info!("{}", msg);
                        per_experiment_messages.push(msg);
                    }
                },
                Err(e) => {
                    let msg = format!(
                        "Experiment '{}': time monotonicity report error: {}",
                        exp_id, e
                    );
                    info!("{}", msg);
                    per_experiment_messages.push(msg);
                }
            }
        }

        let kind = if per_experiment_messages
            .iter()
            .any(|message| message.contains("failed") || message.contains("cannot sort"))
        {
            MessageKind::Error
        } else if per_experiment_messages
            .iter()
            .any(|message| message.contains("duplicates remain"))
        {
            MessageKind::Warning
        } else if per_experiment_messages
            .iter()
            .any(|message| message.contains("sorted by time"))
        {
            MessageKind::Success
        } else {
            MessageKind::Info
        };
        self.set_message(kind, per_experiment_messages.join(" | "));
    }

    /// Renders the data table with checkboxes in header
    /// Each column represents one experiment-column pair
    fn render_table(&mut self, ui: &mut egui::Ui, model: &PlotModel) {
        let experiments = if self.compare_mode {
            let mut selected = self.compare_experiments.iter().cloned().collect::<Vec<_>>();
            selected.sort();
            selected
        } else {
            self.table_experiment.iter().cloned().collect()
        };
        if experiments.is_empty() {
            ui.label(if self.compare_mode {
                "Select experiments to compare"
            } else {
                "No experiments loaded"
            });
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

        let show_row_index = !self.compare_mode;
        let column_width = data_column_width(ui.available_width(), all_columns.len());
        ui.push_id("column_manager_table", |ui| {
            egui::ScrollArea::both()
                .id_salt("column_manager_table_scroll")
                .auto_shrink([false, false])
                .max_height(ui.available_height().max(1.0))
                .show(ui, |ui| {
                    let mut table = TableBuilder::new(ui).striped(true);
                    if show_row_index {
                        table = table.column(Column::initial(58.0).resizable(true));
                    }
                    table
                        .columns(
                            Column::initial(column_width).resizable(false),
                            all_columns.len(),
                        )
                        .header(40.0, |mut header| {
                            if show_row_index {
                                header.col(|ui| {
                                    ui.heading("Row");
                                });
                            }
                            // Each column header: checkbox + clickable title.
                            for (exp_id, col_name) in &all_columns {
                                header.col(|ui| {
                                    ui.vertical(|ui| {
                                        let key = (exp_id.clone(), col_name.clone());
                                        let checked =
                                            self.column_checkboxes.entry(key).or_insert(false);
                                        let is_selected = self
                                            .selected_column
                                            .as_ref()
                                            .is_some_and(|(selected_exp, selected_col)| {
                                                selected_exp == exp_id && selected_col == col_name
                                            });
                                        let title = ui.selectable_label(
                                            is_selected,
                                            format!("{} {}", exp_id, col_name),
                                        );
                                        let checkbox = ui.add_enabled(
                                            !self.compare_mode,
                                            egui::Checkbox::new(checked, ""),
                                        ).on_hover_text(if self.compare_mode {
                                            "Column selection and data-changing actions are disabled in comparison mode."
                                        } else {
                                            "Select this column for the confirmed bulk-delete action. Clicking the column name only opens its details."
                                        });
                                        let _ = checkbox.labelled_by(title.id);
                                        if title.clicked() {
                                            self.selected_column =
                                                Some((exp_id.clone(), col_name.clone()));
                                        }
                                        if let Ok(Some(provenance)) =
                                            model.column_provenance_for_experiment(exp_id, col_name)
                                        {
                                            title.on_hover_text(provenance);
                                        }
                                    });
                                });
                            }
                        })
                        .body(|mut body| {
                            let preview_lengths = experiments
                                .iter()
                                .filter_map(|id| {
                                    self.cached_samples
                                        .get(id)
                                        .map(|preview| preview.row_indices.len())
                                })
                                .collect::<Vec<_>>();
                            let row_count = if self.compare_mode {
                                preview_lengths.into_iter().min().unwrap_or(0)
                            } else {
                                preview_lengths.into_iter().max().unwrap_or(0)
                            };

                            for row_idx in 0..row_count {
                                body.row(18.0, |mut row| {
                                    if show_row_index {
                                        row.col(|ui| {
                                            let indices = experiments
                                                .iter()
                                                .filter_map(|id| {
                                                    let preview = self.cached_samples.get(id)?;
                                                    preview
                                                        .row_indices
                                                        .get(row_idx)
                                                        .map(|source_row| source_row.to_string())
                                                })
                                                .collect::<Vec<_>>()
                                                .join(" | ");
                                            ui.label(indices);
                                        });
                                    }
                                    for (exp_id, col_name) in &all_columns {
                                        row.col(|ui| {
                                            if let Some(preview) = self.cached_samples.get(exp_id) {
                                                let sample_idx =
                                                    if self.compare_mode && row_count > 1 {
                                                        row_idx
                                                            * preview
                                                                .row_indices
                                                                .len()
                                                                .saturating_sub(1)
                                                            / (row_count - 1)
                                                    } else {
                                                        row_idx
                                                    };
                                                if let Some(values) = preview.columns.get(col_name)
                                                {
                                                    match values.get(sample_idx).copied().flatten()
                                                    {
                                                        Some(val) => {
                                                            ui.label(format_sample_value(val))
                                                        }
                                                        None => ui.label("null"),
                                                    };
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
        });
    }

    fn render_selected_column_details(&mut self, ui: &mut egui::Ui, model: &mut PlotModel) {
        ui.push_id("column_manager_details", |ui| {
            ui.separator();
            ui.horizontal(|ui| {
                ui.heading("Column details");
                if ui.small_button("Clear").on_hover_text("Close the column inspector; this does not change data.").clicked() {
                    self.selected_column = None;
                }
            });

            let Some((exp_id, col_name)) = self.selected_column.clone() else {
                ui.label("Click a column title to inspect its provenance and operation history.");
                return;
            };

            ui.horizontal(|ui| {
                for (tab, label) in [(DetailsTab::Overview, "Overview"), (DetailsTab::Processing, "Processing"), (DetailsTab::Quality, "Quality"), (DetailsTab::History, "History")] {
                    ui.selectable_value(&mut self.details_tab, tab, label);
                }
            });

            let details_height = ui.available_height().max(140.0);
            egui::ScrollArea::vertical()
                .id_salt(format!(
                    "column_manager_details_scroll:{}:{}",
                    exp_id, col_name
                ))
                .max_height(details_height)
                .show(ui, |ui| {
                    if self.details_tab == DetailsTab::Overview {
                        ui.strong(format!("{} / {}", display_experiment_id(&exp_id), col_name));
                        if let Ok(experiment) = model.series.get_experiment_by_id(&exp_id) {
                            if let Some(meta) = experiment.dataset.schema.columns.get(&col_name) {
                                ui.label(format!("Unit: {:?} | Nature: {:?}", meta.unit, meta.nature));
                                let bindings = [("time", experiment.dataset.schema.time.as_deref()), ("temperature", experiment.dataset.schema.temperature.as_deref()), ("mass", experiment.dataset.schema.mass.as_deref()), ("conversion", experiment.dataset.schema.alpha.as_deref()), ("dimensionless mass", experiment.dataset.schema.eta.as_deref())];
                                let bound = bindings.iter().filter_map(|(role, name)| (*name == Some(col_name.as_str())).then_some(*role)).collect::<Vec<_>>();
                                ui.label(format!("Schema binding: {}", if bound.is_empty() { "none".to_string() } else { bound.join(", ") }));
                            }
                        }
                        ui.label(format!("On plot: {}", model.plots.iter().any(|curve| curve.experiment_id == exp_id && (curve.x_name == col_name || curve.y_name == col_name))));
                        let mass_active = model.active_mass_source_for_experiment(&exp_id).is_ok_and(|source| source == col_name);
                        let temp_active = model.active_temperature_source_for_experiment(&exp_id).is_ok_and(|source| source == col_name);
                        ui.label(format!("Calculation source: {}", match (mass_active, temp_active) { (true, true) => "mass and temperature", (true, false) => "mass", (false, true) => "temperature", _ => "none" }));
                        if self.compare_mode { ui.label("Comparison mode is read-only."); }
                        else {
                            ui.horizontal(|ui| {
                                let label = ui.label("Rename to:");
                                let input = ui.text_edit_singleline(&mut self.new_column_name);
                                let _ = input.labelled_by(label.id);
                                if ui.add_enabled(!self.new_column_name.trim().is_empty(), egui::Button::new("Rename"))
                                    .on_hover_text("Renames this column; no new column is created.")
                                    .clicked() { self.handle_rename_column(model, &exp_id, &col_name); }
                            });
                        }
                    }

                    let provenance = model
                .series
                .get_experiment_by_id(&exp_id)
                .ok()
                .and_then(|experiment| experiment.dataset.column_provenance(&col_name))
                .cloned();
                    let history = model
                        .column_history_feed_for_experiment(&exp_id, &col_name)
                        .unwrap_or_else(|e| format!("Failed to read history: {:?}", e));

                    if self.details_tab == DetailsTab::Processing { ui.label("Provenance"); }
            if let Some(provenance) = provenance {
                if self.details_tab == DetailsTab::Processing { ui.label(format!("Source columns: {}", provenance.root_columns.join(", "))); }
                for (index, step) in provenance.steps.iter().enumerate() {
                    if self.details_tab == DetailsTab::Processing {
                    ui.separator();
                    ui.label(format!(
                        "Step {}: {} -> {}",
                        index + 1,
                        step.input_columns.join(", "),
                        step.output_column
                    ));
                    ui.label(step.kind.to_string());
                    }

                    if self.details_tab == DetailsTab::Quality { if let Some(quality) = &step.quality {
                        ui.label(format!("Quality: {}", step.kind));
                        if let Some(report) = &quality.report {
                        ui.label(format!("Compared {} vs {} ({}, {})", report.column_raw, report.column_filtered, quality.compared_to, quality.status));
                        egui::Grid::new(("column_quality", &exp_id, &col_name, index))
                            .striped(true)
                            .show(ui, |ui| {
                                ui.strong("Metric");
                                ui.strong("Value");
                                ui.strong("Interpretation");
                                ui.end_row();

                                for (name, value, meaning) in [
                                    (
                                        "RMSE",
                                        report.rmse,
                                        "Difference from the comparison column; lower is closer.",
                                    ),
                                    (
                                        "Normalized RMSE",
                                        report.normalized_rmse,
                                        "Scale-normalized difference; lower is closer.",
                                    ),
                                    (
                                        "Correlation",
                                        report.correlation,
                                        "Shape similarity; values closer to 1 indicate stronger agreement.",
                                    ),
                                    (
                                        "First-difference roughness ratio",
                                        report.roughness_ratio_1st,
                                        "Output/input first-difference roughness; below 1 is smoother.",
                                    ),
                                    (
                                        "Second-difference roughness ratio",
                                        report.roughness_ratio_2nd,
                                        "Output/input second-difference roughness; below 1 is smoother.",
                                    ),
                                ] {
                                    ui.label(name);
                                    ui.label(format!("{value:.6}"));
                                    ui.label(meaning);
                                    ui.end_row();
                                }
                            });
                        ui.label(
                            "These metrics describe change relative to the comparison column; they do not establish accuracy against a true signal.",
                        );
                        } else {
                            ui.label(format!("Quality metrics unavailable: {}", quality.status));
                        }
                    }}
                }
            } else {
                if matches!(self.details_tab, DetailsTab::Processing | DetailsTab::Quality) { ui.label("No provenance available"); }
            }

                    if self.details_tab == DetailsTab::History {
                        ui.label("Operation history");
                        ui.monospace(history);
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

                egui::ComboBox::from_label("Select Experiment")
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

                let comment_label = ui.label("Comment:").on_hover_text(
                    "Optional experiment note. Leave blank to keep the current comment.",
                );
                ui.text_edit_multiline(&mut self.comment_input)
                    .labelled_by(comment_label.id);

                let rate_label = ui.label("Heating Rate (K/min):").on_hover_text(
                    "Optional numeric metadata in K/min; leave blank to keep the current value.",
                );
                ui.text_edit_singleline(&mut self.heating_rate_input)
                    .labelled_by(rate_label.id);

                let temperature_label = ui.label("Temperature (K):").on_hover_text(
                    "Optional numeric metadata in K; leave blank to keep the current value.",
                );
                ui.text_edit_singleline(&mut self.temperature_input)
                    .labelled_by(temperature_label.id);

                ui.separator();

                if ui
                    .button("Apply")
                    .on_hover_text(
                        "Saves experiment metadata only; no data columns are created or changed.",
                    )
                    .clicked()
                {
                    if let Some(exp_id) = self.selected_experiment.clone() {
                        match apply_metadata_draft(
                            &mut model.series,
                            &exp_id,
                            &self.comment_input,
                            &self.heating_rate_input,
                            &self.temperature_input,
                        ) {
                            Ok(()) => {
                                self.set_message(
                                    MessageKind::Success,
                                    "Metadata updated successfully",
                                );
                                self.metadata_window_open = false;
                            }
                            Err(error) => self.set_message(
                                MessageKind::Error,
                                format!("Metadata update failed: {error}"),
                            ),
                        }
                    } else {
                        self.set_message(MessageKind::Error, "Select an experiment first");
                    }
                }
            });
    }
}

fn parse_optional_finite_positive(input: &str, label: &str) -> Result<Option<f64>, String> {
    if input.trim().is_empty() {
        return Ok(None);
    }
    let value = input
        .trim()
        .parse::<f64>()
        .map_err(|_| format!("{label} must be a number"))?;
    if !value.is_finite() || value <= 0.0 {
        return Err(format!("{label} must be finite and greater than zero"));
    }
    Ok(Some(value))
}

fn apply_metadata_draft(
    series: &mut TGASeries,
    experiment_id: &str,
    comment: &str,
    heating_rate_input: &str,
    temperature_input: &str,
) -> Result<(), String> {
    // Parse the complete draft before touching the experiment, then commit once.
    let heating_rate = parse_optional_finite_positive(heating_rate_input, "Heating rate")?;
    let temperature = parse_optional_finite_positive(temperature_input, "Temperature")?;
    if heating_rate.is_some() && temperature.is_some() {
        return Err("Set either heating rate or isothermal temperature, not both".to_string());
    }

    series
        .try_mutate_by_id(experiment_id, |experiment| {
            if !comment.trim().is_empty() {
                experiment.meta.comment = Some(comment.trim().to_string());
            }
            if let Some(rate) = heating_rate {
                experiment.meta.heating_rate = Some(rate);
                experiment.meta.thermal_mode = ExperimentThermalMode::NonIsothermal;
                experiment.meta.isothermal_temperature = None;
            }
            if let Some(temperature) = temperature {
                experiment.meta.isothermal_temperature = Some(temperature);
                experiment.meta.thermal_mode = ExperimentThermalMode::Isothermal;
                experiment.meta.heating_rate = None;
            }
            Ok(())
        })
        .map_err(|error| error.to_string())
}

#[cfg(test)]
mod metadata_input_tests {
    use super::{apply_metadata_draft, parse_optional_finite_positive};
    use crate::Kinetics::experimental_kinetics::experiment_series_main::{
        ExperimentMeta, ExperimentThermalMode, TGAExperiment, TGASeries,
    };
    use crate::Kinetics::experimental_kinetics::testing_mod::VirtualTGA;
    use crate::gui::experimental_kinetics_gui::model::PlotModel;
    use egui_kittest::Harness;
    use egui_kittest::kittest::Queryable;

    struct MetadataUiState {
        model: PlotModel,
        manager: super::ColumnManagerState,
    }

    fn series_with_metadata() -> TGASeries {
        let mut series = TGASeries::new();
        let data = VirtualTGA {
            time: vec![0.0, 1.0, 2.0],
            temperature: vec![500.0; 3],
            mass: vec![10.0, 9.0, 8.0],
        };
        let dataset = crate::Kinetics::experimental_kinetics::one_experiment_dataset::TGADataset::create_from_synthetic_data(&data).unwrap();
        series
            .push(TGAExperiment::new(dataset).with_id("meta"))
            .unwrap();
        series
            .try_mutate_by_id("meta", |experiment| {
                experiment.meta = ExperimentMeta {
                    id: "meta".into(),
                    heating_rate: Some(4.0),
                    isothermal_temperature: None,
                    thermal_mode: ExperimentThermalMode::NonIsothermal,
                    comment: Some("keep me".into()),
                };
                Ok(())
            })
            .unwrap();
        series
    }

    #[test]
    fn metadata_numeric_inputs_reject_invalid_values_before_commit() {
        assert_eq!(
            parse_optional_finite_positive("", "Temperature").unwrap(),
            None
        );
        assert_eq!(
            parse_optional_finite_positive("523.15", "Temperature").unwrap(),
            Some(523.15)
        );
        for input in ["abc", "NaN", "inf", "0", "-2"] {
            assert!(
                parse_optional_finite_positive(input, "Temperature").is_err(),
                "accepted {input}"
            );
        }
    }

    #[test]
    fn invalid_metadata_draft_preserves_all_existing_fields() {
        let mut series = series_with_metadata();
        let before = series.get_experiment_by_id("meta").unwrap().meta.clone();

        let error = apply_metadata_draft(&mut series, "meta", "changed", "5", "NaN").unwrap_err();

        assert!(error.contains("finite and greater than zero"));
        let after = &series.get_experiment_by_id("meta").unwrap().meta;
        assert_eq!(after.id, before.id);
        assert_eq!(after.heating_rate, before.heating_rate);
        assert_eq!(after.isothermal_temperature, before.isothermal_temperature);
        assert_eq!(after.thermal_mode, before.thermal_mode);
        assert_eq!(after.comment, before.comment);
    }

    #[test]
    fn valid_metadata_draft_commits_mode_and_comment_together() {
        let mut series = series_with_metadata();
        apply_metadata_draft(&mut series, "meta", "  updated  ", "", "523.15").unwrap();

        let meta = &series.get_experiment_by_id("meta").unwrap().meta;
        assert_eq!(meta.comment.as_deref(), Some("updated"));
        assert_eq!(meta.thermal_mode, ExperimentThermalMode::Isothermal);
        assert_eq!(meta.isothermal_temperature, Some(523.15));
        assert!(meta.heating_rate.is_none());
    }

    #[test]
    fn egui_apply_with_invalid_metadata_keeps_the_original_record() {
        let mut model = PlotModel::default();
        model.series = series_with_metadata();
        let mut manager = super::ColumnManagerState::new();
        manager.metadata_window_open = true;
        manager.selected_experiment = Some("meta".into());
        manager.comment_input = "must not partially apply".into();
        manager.heating_rate_input = "5".into();
        manager.temperature_input = "NaN".into();
        let before = model
            .series
            .get_experiment_by_id("meta")
            .unwrap()
            .meta
            .clone();
        let mut harness = Harness::new_ui_state(
            |ui, state: &mut MetadataUiState| {
                state
                    .manager
                    .show_metadata_window(ui.ctx(), &mut state.model);
            },
            MetadataUiState { model, manager },
        );

        harness.get_by_label("Apply").click();
        harness.run();

        let state = harness.state();
        let after = &state
            .model
            .series
            .get_experiment_by_id("meta")
            .unwrap()
            .meta;
        assert_eq!(after.comment, before.comment);
        assert_eq!(after.heating_rate, before.heating_rate);
        assert_eq!(after.isothermal_temperature, before.isothermal_temperature);
        assert_eq!(after.thermal_mode, before.thermal_mode);
        assert!(
            state
                .manager
                .system_message
                .as_ref()
                .unwrap()
                .text
                .contains("Metadata update failed")
        );
    }

    #[test]
    fn egui_apply_commits_valid_metadata_and_comment_together() {
        let mut model = PlotModel::default();
        model.series = series_with_metadata();
        let mut manager = super::ColumnManagerState::new();
        manager.metadata_window_open = true;
        manager.selected_experiment = Some("meta".into());
        manager.comment_input = "updated from GUI".into();
        manager.temperature_input = "523.15".into();
        let mut harness = Harness::new_ui_state(
            |ui, state: &mut MetadataUiState| {
                state
                    .manager
                    .show_metadata_window(ui.ctx(), &mut state.model);
            },
            MetadataUiState { model, manager },
        );

        harness.get_by_label("Apply").click();
        harness.run();

        let state = harness.state();
        let meta = &state
            .model
            .series
            .get_experiment_by_id("meta")
            .unwrap()
            .meta;
        assert_eq!(meta.comment.as_deref(), Some("updated from GUI"));
        assert_eq!(meta.thermal_mode, ExperimentThermalMode::Isothermal);
        assert_eq!(meta.isothermal_temperature, Some(523.15));
        assert!(meta.heating_rate.is_none());
        assert!(
            state
                .manager
                .system_message
                .as_ref()
                .unwrap()
                .text
                .contains("Metadata updated successfully")
        );
    }
}

/// Main window function for Column Manager
pub fn show_column_manager_window(
    ctx: &egui::Context,
    model: &mut PlotModel,
    state: &mut ColumnManagerState,
    methods_state: &mut KineticMethodsWindowState,
    open: &mut bool,
) -> bool {
    state.show_metadata_window(ctx, model);
    let mut clear_workspace = false;
    // state.system_message.clear();
    egui::Window::new("Column Manager")
        .open(open)
        .resizable(true)
        .default_size([800.0, 600.0])
        .show(ctx, |ui| {
            let mut commit_reconstruction = false;
            let mut discard_reconstruction = false;
            let mut confirm_columns = false;
            let mut cancel_columns = false;
            if let Err(e) = state.update_cache(model) {
                let text = format!("Could not refresh column samples: {e}");
                state.set_message_if_changed(MessageKind::Error, text.clone());
                ui.colored_label(egui::Color32::LIGHT_RED, text);
                return;
            }

            if let Some(batch) = methods_state.fit_pending_batch.as_ref() {
                ui.separator();
                ui.heading("Fit reconstruction preview");
                ui.label(format!(
                    "Source: {} | X: {} | {} rows | {} fitted column(s) | {:?}",
                    batch.source_experiment_id,
                    batch.x_col,
                    batch.x_values.len(),
                    batch.columns.len(),
                    batch.mode,
                ));
                for column in &batch.columns {
                    let coefficients = column.fit.coefficients.iter()
                        .map(|(name, value)| format!("{name}={value:.5e}"))
                        .collect::<Vec<_>>().join(", ");
                    ui.label(format!(
                        "{}: R2={:.6}, model={}, Y={}, parameters: {}",
                        column.name,
                        column.fit.r2,
                        column.fit.kinetic_model.display_name(),
                        column.fit.y_col,
                        coefficients,
                    ));
                }
                ui.horizontal(|ui| {
                    let label = ui.label("New experiment ID:");
                    let input = ui.text_edit_singleline(&mut methods_state.fit_reconstructed_experiment_id);
                    let _ = input.labelled_by(label.id);
                    if ui.button("Add reconstructed experiment")
                        .on_hover_text("Creates a separate complete experiment and plots its fitted columns. The source experiment is not modified.")
                        .clicked()
                    {
                        commit_reconstruction = true;
                    }
                    if ui.button("Discard preview")
                        .on_hover_text("Discard this evaluation and its temporary fit sessions; source experiment data remains unchanged.")
                        .clicked() {
                        discard_reconstruction = true;
                    }
                });
                egui::CollapsingHeader::new("Sample values")
                    .default_open(true)
                    .show(ui, |ui| {
                        let (start, end) = match batch.mode {
                            crate::Kinetics::experimental_kinetics::math::fitting::FitEvaluationMode::Extrapolate { start, end, .. }
                            | crate::Kinetics::experimental_kinetics::math::fitting::FitEvaluationMode::FillGap { start, end, .. } => (start, end),
                        };
                        let preview_indices = batch.x_values.iter().enumerate()
                            .filter(|(_, x)| **x >= start && **x <= end)
                            .map(|(index, _)| index)
                            .take(8)
                            .collect::<Vec<_>>();
                        egui::ScrollArea::both().max_height(140.0).show(ui, |ui| {
                            egui::Grid::new("fit_reconstruction_preview_values")
                                .striped(true)
                                .show(ui, |ui| {
                                    ui.strong(&batch.x_col);
                                    for column in &batch.columns { ui.strong(&column.name); }
                                    ui.end_row();
                                    for index in &preview_indices {
                                        ui.label(format!("{:.6}", batch.x_values[*index]));
                                        for column in &batch.columns {
                                            ui.label(format!("{:.6}", column.values[*index]));
                                        }
                                        ui.end_row();
                                    }
                                });
                        });
                    });
            }
            if commit_reconstruction {
                if let Some(batch) = methods_state.fit_pending_batch.as_ref() {
                    let id = methods_state.fit_reconstructed_experiment_id.trim().to_string();
                    match model.push_fitted_evaluation_batch_and_plot(batch, id) {
                        Ok(()) => {
                            methods_state.fit_pending_batch = None;
                            methods_state.fit_sessions.clear();
                            state.set_message(MessageKind::Success, "Reconstructed experiment added; source data unchanged");
                        }
                        Err(crate::gui::experimental_kinetics_gui::model::TGAGUIError::PostCommitWarning(warning)) => {
                            methods_state.fit_pending_batch = None;
                            methods_state.fit_sessions.clear();
                            state.set_message(MessageKind::Warning, warning);
                        }
                        Err(error) => state.set_message(MessageKind::Error, format!("Reconstruction was not added: {error:?}")),
                    }
                }
            }
            if discard_reconstruction {
                methods_state.fit_pending_batch = None;
                methods_state.fit_sessions.clear();
                methods_state.fit_model_result = None;
            }

            if !state.confirm_drop_columns.is_empty() {
                egui::Window::new("Confirm Drop Columns")
                    .collapsible(false)
                    .resizable(false)
                    .show(ctx, |ui| {
                        ui.label(format!("Delete {} column(s) from {}?", state.confirm_drop_columns.len(), display_experiment_id(&state.confirm_drop_columns[0].0)));
                        for (_, name) in &state.confirm_drop_columns {
                            ui.label(name);
                        }
                        ui.horizontal(|ui| {
                            if ui.button("Delete columns")
                                .on_hover_text("Permanently removes the listed columns from this experiment; undo may be available from History.")
                                .clicked() { confirm_columns = true; }
                            if ui.button("Cancel").on_hover_text("Close this confirmation without changing data.").clicked() { cancel_columns = true; }
                        });
                    });
            }

            if state.confirm_delete_all {
                egui::Window::new("Confirm Delete All Workspace")
                    .collapsible(false)
                    .resizable(false)
                    .show(ctx, |ui| {
                        ui.heading(
                            egui::RichText::new("DELETE ALL WORKSPACE DATA?")
                                .size(24.0)
                                .strong()
                                .color(egui::Color32::LIGHT_RED),
                        );
                        ui.label(format!(
                            "This will permanently remove {} experiment(s) and {} plotted curve(s).",
                            model.list_of_experiments().len(),
                            model.plots.len(),
                        ));
                        ui.label("This action cannot be undone.");
                        ui.horizontal(|ui| {
                            if ui
                                .add_sized(
                                    [230.0, 48.0],
                                    egui::Button::new(
                                        egui::RichText::new("DELETE EVERYTHING")
                                            .strong()
                                            .color(egui::Color32::WHITE),
                                    )
                                    .fill(egui::Color32::DARK_RED),
                                )
                                .clicked()
                            {
                                clear_workspace = true;
                            }
                            if ui.button("Cancel").clicked() {
                                state.confirm_delete_all = false;
                            }
                        });
                    });
            }
            if confirm_columns {
                state.commit_drop_columns(model);
            }
            if cancel_columns {
                state.confirm_drop_columns.clear();
            }

            // Handle experiment drop confirmation dialog
            if let Some(exp_id) = state.confirm_drop_experiment.clone() {
                let mut confirmed = false;
                let mut cancelled = false;

                egui::Window::new("Confirm Drop Experiment")
                    .collapsible(false)
                    .resizable(false)
                    .show(ctx, |ui| {
                        ui.label(format!("Delete experiment {} and its plots?", display_experiment_id(&exp_id)));
                        ui.horizontal(|ui| {
                            if ui.button("Yes").on_hover_text("Deletes this experiment and all its plots.").clicked() {
                                confirmed = true;
                            }
                            if ui.button("No").on_hover_text("Keep the experiment unchanged.").clicked() {
                                cancelled = true;
                            }
                        });
                    });

                if confirmed {
                    if let Err(e) = model.drop_experiment_by_id(&exp_id) {
                        state.set_message(MessageKind::Error, format!("Could not drop experiment: {e:?}"));
                    } else {
                        state.set_message(MessageKind::Success, format!("Experiment {} dropped successfully", exp_id));
                        // Clear cached data for this experiment
                        state.cached_samples.remove(&exp_id);
                        // Remove all checkboxes for this experiment
                        state.column_checkboxes.retain(|(eid, _), _| eid != &exp_id);
                        if state
                            .selected_column
                            .as_ref()
                            .map(|(selected_exp, _)| selected_exp == &exp_id)
                            .unwrap_or(false)
                        {
                            state.selected_column = None;
                        }
                    }
                    state.confirm_drop_experiment = None;
                }

                if cancelled {
                    state.confirm_drop_experiment = None;
                }
            }

            state.render_controls(ui, model);
            let max_samples_height = (ui.available_height() - 120.0).clamp(180.0, 760.0);
            egui::Resize::default()
                .id_salt("column_manager_samples_resize")
                .default_height((max_samples_height * 0.72).max(180.0))
                .min_size(egui::vec2(300.0, 110.0))
                .max_height(max_samples_height)
                .resizable(egui::Vec2b::new(false, true))
                .show(ui, |ui| state.render_table(ui, model));
            state.render_selected_column_details(ui, model);
        });
    if clear_workspace {
        model.clear_workspace();
        model.message = "Workspace cleared. All experiments and plots were removed.".to_string();
        *state = ColumnManagerState::default();
        *methods_state = KineticMethodsWindowState::default();
        *open = false;
        return true;
    }
    state.publish_message_to_main_status(model);
    false
}
