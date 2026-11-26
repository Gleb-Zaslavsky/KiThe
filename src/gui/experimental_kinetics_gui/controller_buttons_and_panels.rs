//! Top menu bar: File Manager, Math, Kinetic Methods,
//! Direct Problem and Test Options with stub handlers.
#![allow(non_camel_case_types)]
#![allow(non_snake_case)]
use crate::Kinetics::experimental_kinetics::exp_engine_api::XY;
use crate::Kinetics::experimental_kinetics::experiment_series_main::{
    ExperimentMeta, ExperimentThermalMode, TGAExperiment, TGASeries,
};
use crate::Kinetics::experimental_kinetics::one_experiment_dataset::{
    ColumnNature, TGADataset, TGADomainError, Unit,
};
use crate::gui::experimental_kinetics_gui::model::PlotModel;
use crate::gui::experimental_kinetics_gui::model::TGAGUIError;

use crate::gui::experimental_kinetics_gui::controller_filters::Mathematics;
use crate::gui::experimental_kinetics_gui::controller_golden_pipeline::GoldenPipelineDialogState;
use crate::gui::experimental_kinetics_gui::controller_kinetics::DirectProblemDialogState;
use crate::gui::experimental_kinetics_gui::controller_methods::{
    KineticMethods, KineticMethodsWindowState,
};
use crate::gui::experimental_kinetics_gui::kitheplot_wrapper::KiThePlotWindowState;

use crate::gui::experimental_kinetics_gui::test_options::TestOptions;
use eframe::egui;
use log::{info, warn};
use polars::prelude::col;
use std::collections::{HashMap, HashSet};
use std::path::Path;
use std::path::PathBuf;
use std::{
    fs::File,
    io::{BufRead, BufReader},
};

fn write_atomically<T>(
    path: &Path,
    write: impl FnOnce(&Path) -> Result<T, String>,
) -> Result<T, String> {
    let parent = path
        .parent()
        .filter(|parent| !parent.as_os_str().is_empty())
        .unwrap_or_else(|| Path::new("."));
    let temporary_file = tempfile::NamedTempFile::new_in(parent)
        .map_err(|error| format!("Could not create temporary file: {error}"))?;
    let result = write(temporary_file.path())?;
    temporary_file
        .as_file()
        .sync_all()
        .map_err(|error| format!("Could not flush temporary file: {error}"))?;
    temporary_file
        .persist(path)
        .map_err(|error| format!("Could not replace '{}': {}", path.display(), error.error))?;
    Ok(result)
}

//===================================================================================
//  DROP DOWN MENUE AT THE TOP OF WINDOW
//====================================================================================
/// High-level entrypoint for top menus and panels.
pub struct TopDropDownMenues;

#[derive(Clone, Copy, PartialEq, Eq)]
enum TopMenuShortcut {
    File,
    Math,
    KineticMethods,
}

impl TopDropDownMenues {
    /// Render the top menus directly inside a parent `Ui`.
    pub fn top_menus(
        ui: &mut egui::Ui,
        model: &mut PlotModel,
        mathematics: &mut Mathematics,
        new_experiment_dialog: &mut NewExperimentDialogState,
        test_options: &mut TestOptions,
        kithe_plot_window: &mut KiThePlotWindowState,
        kinetic_methods: &mut KineticMethodsWindowState,
        direct_problem_state: &mut DirectProblemDialogState,
    ) {
        let ctx = ui.ctx().clone();
        let mut requested_menu = None;
        if !ctx.egui_wants_keyboard_input() && !egui::Popup::is_any_open(&ctx) {
            let (open_import, save, save_as, undo, open_column_manager, requested_menu_key) = ctx
                .input_mut(|input| {
                    let open_import = input.consume_key(egui::Modifiers::COMMAND, egui::Key::O);
                    let save_as = input.consume_key(
                        egui::Modifiers::COMMAND.plus(egui::Modifiers::SHIFT),
                        egui::Key::S,
                    );
                    let save =
                        !save_as && input.consume_key(egui::Modifiers::COMMAND, egui::Key::S);
                    let undo = input.consume_key(egui::Modifiers::COMMAND, egui::Key::Z);
                    let open_column_manager =
                        input.consume_key(egui::Modifiers::NONE, egui::Key::F1);
                    let file = input.consume_key(egui::Modifiers::NONE, egui::Key::F2);
                    let math = input.consume_key(egui::Modifiers::NONE, egui::Key::F3);
                    let kinetics = input.consume_key(egui::Modifiers::NONE, egui::Key::F4);
                    let requested = if file {
                        Some(TopMenuShortcut::File)
                    } else if math {
                        Some(TopMenuShortcut::Math)
                    } else if kinetics {
                        Some(TopMenuShortcut::KineticMethods)
                    } else {
                        None
                    };
                    (
                        open_import,
                        save,
                        save_as,
                        undo,
                        open_column_manager,
                        requested,
                    )
                });
            requested_menu = requested_menu_key;

            if open_import {
                new_experiment_dialog.open_document_dialog(model);
            }
            if save || save_as {
                new_experiment_dialog.save_series_shortcut(model, save_as);
            }
            if undo {
                Self::undo_selected(model);
            }
            if open_column_manager {
                new_experiment_dialog.column_manager_open = true;
            }
        }

        ui.horizontal_wrapped(|ui| {
            Self::menu_button_with_shortcut(
                ui,
                "File Manager",
                "F2",
                requested_menu == Some(TopMenuShortcut::File),
                |ui| {
                    FileManager::show(ui, model, new_experiment_dialog, kithe_plot_window);
                },
            );

            Self::menu_button_with_shortcut(
                ui,
                "Math",
                "F3",
                requested_menu == Some(TopMenuShortcut::Math),
                |ui| {
                    mathematics.show_menu(ui);
                },
            );

            Self::menu_button_with_shortcut(
                ui,
                "Kinetic Methods",
                "F4",
                requested_menu == Some(TopMenuShortcut::KineticMethods),
                |ui| {
                    KineticMethods::show_menu(ui, kinetic_methods);
                },
            );

            ui.menu_button("Direct Problem", |ui| {
                if ui.button("Simple One-Step Kinetics").clicked() {
                    direct_problem_state.open();
                    ui.close();
                }
                if ui.button("Complex Direct Problem").on_hover_text("Solve a reaction network for a series of temperatures or heating rates and add the runs as experiments.").clicked() {
                    direct_problem_state.network_open = true;
                    ui.close();
                }
            });

            ui.menu_button("Test Options", |ui| {
                TestOptions::show(ui, model, test_options);
            });
        });
        ui.separator();
        mathematics.show_windows(ui.ctx(), model);
        // Show any open Kinetic Methods windows
        let _ = KineticMethods::show_windows(ui, model, kinetic_methods);
    }

    fn menu_button_with_shortcut(
        ui: &mut egui::Ui,
        label: &str,
        shortcut: &str,
        requested_open: bool,
        add_contents: impl FnOnce(&mut egui::Ui),
    ) {
        let response = ui
            .button(label)
            .on_hover_text(format!("Open {label} menu ({shortcut})"));
        let open_command = if requested_open {
            Some(egui::SetOpenCommand::Bool(true))
        } else if response.clicked() {
            Some(egui::SetOpenCommand::Toggle)
        } else {
            None
        };
        egui::Popup::menu(&response)
            .open_memory(open_command)
            .show(add_contents);
    }

    fn undo_selected(model: &mut PlotModel) {
        match model.reverse_selected() {
            Ok(()) => {
                let _ = model.push_message("Last operation undone");
            }
            Err(err) => {
                warn!("Undo failed: {err}");
                let _ = model.push_message(&format!("Undo failed: {err}"));
            }
        }
    }
}

// Group: File Manager
struct FileManager;

impl FileManager {
    fn show(
        ui: &mut egui::Ui,
        model: &mut PlotModel,
        new_experiment_dialog: &mut NewExperimentDialogState,
        kithe_plot_window: &mut KiThePlotWindowState,
    ) {
        if ui
            .button("Open Document")
            .on_hover_text("Choose a typed KiThe series CSV to replace the current document after the unsaved-changes check. Raw files and typed datasets belong in Import Data. Shortcut: Ctrl/Cmd+O.")
            .clicked()
        {
            new_experiment_dialog.open_document_dialog(model);
        }
        if ui
            .button("Import Data")
            .on_hover_text(
                "Import raw TXT/CSV with manual bindings or add experiments from a typed KiThe series.",
            )
            .clicked()
        {
            new_experiment_dialog.choose_import_file();
        }
        if ui.button("Append to experiment")
            .on_hover_text("Import a continuation, preview its time alignment, then append it or replace the selected experiment.")
            .clicked()
        {
            new_experiment_dialog.choose_append_file(model);
        }

        if ui
            .button("Manage Plot")
            .on_hover_text("Show or hide plot curves and manage their selection.")
            .clicked()
        {
            new_experiment_dialog.manage_plot_open = true;
            ui.ctx().request_repaint();
        }
        ui.separator();
        let save_label = if new_experiment_dialog.is_document_dirty(model) {
            "Save *"
        } else {
            "Save"
        };
        if ui
            .button(save_label)
            .on_hover_text(
                "Save the complete current document, including every experiment and column. Shortcut: Ctrl/Cmd+S.",
            )
            .clicked()
        {
            new_experiment_dialog.save_series_shortcut(model, false);
        }
        if ui
            .button("Save As")
            .on_hover_text("Choose a destination for the complete typed KiThe document.")
            .clicked()
        {
            new_experiment_dialog.open_save_document_dialog(model);
            ui.ctx().request_repaint();
        }
        if ui
            .button("Export CSV")
            .on_hover_text("Export the full series or one experiment as a typed KiThe CSV.")
            .clicked()
        {
            new_experiment_dialog.open_export_dialog(SaveSeriesFormat::Csv, model);
            ui.ctx().request_repaint();
        }
        if ui
            .button("Export TXT")
            .on_hover_text("Export selected columns of one experiment as plain tab-delimited text.")
            .clicked()
        {
            new_experiment_dialog.open_export_dialog(SaveSeriesFormat::Txt, model);
            ui.ctx().request_repaint();
        }
        if ui
            .button("Export Excel")
            .on_hover_text("Export selected experiments and columns to an Excel workbook (.xlsx), one worksheet per experiment.")
            .clicked()
        {
            new_experiment_dialog.open_export_dialog(SaveSeriesFormat::Xlsx, model);
            ui.ctx().request_repaint();
        }
        if ui
            .button("Plot Editor")
            .on_hover_text(
                "Open the print-ready plot editor for journal/report figures. This does not change the live interactive plots or experiment data.",
            )
            .clicked()
        {
            Self::save_as_image(ui.ctx(), model, kithe_plot_window);
        }
    }

    fn save_as_image(
        ctx: &egui::Context,
        model: &mut PlotModel,
        kithe_plot_window: &mut KiThePlotWindowState,
    ) {
        match kithe_plot_window.open_from_model(ctx, model) {
            Ok(()) => {
                let _ = model.push_message("KiThe Plot Redactor window opened");
            }
            Err(err) => {
                let message = err.to_string();
                let _ = model.push_message(&format!("Opening Plot Editor failed: {message}"));
            }
        }
    }
}

//=================================================================================================
//  QUICK ACTION PANEL
//=================================================================================================
#[derive(Default)]
pub struct QuickActionPanelState {
    pub input_value: f64,
    pub apply_only_to_selected_chart: bool,
    pub input_string: Option<String>,
    pub x_or_y: Option<XY>, //XY
    pub plot_recreation_required: bool,
    pub golden_pipeline_dialog: GoldenPipelineDialogState,
}

pub struct WrightPanelControllers;

impl WrightPanelControllers {
    fn quick_action_tooltip(action: &str) -> &'static str {
        match action {
            "Manage Plots" => "Open plot visibility and selection controls; no data changes.",
            "Column Manager" => {
                "Inspect columns, samples, provenance and filter-quality reports; no data changes. Shortcut: F1."
            }
            "golden pipeline" => {
                "Open the preprocessing pipeline; applying it may add or modify columns as configured."
            }
            "Clear Selected" => "Clear the current plot selection; no data changes.",
            "Reset View" => "Reset plot axes; no data changes.",
            "Refresh" => "Redraw plots; no data changes.",
            "cut before time" => {
                "Uses Input value as the time threshold. Removes earlier rows from the selected experiment in place; no new column."
            }
            "cut selected" => {
                "Choose X or Y, then select a plot range. Removes rows inside that range (including its boundaries) from the selected experiment in place; no new column."
            }
            "Zoom To Selection" => "Fit the plot view to the selected region; no data changes.",
            "move time to zero" => {
                "Subtracts the first time value from the bound time column in place; no new column."
            }
            "Delete Plot" => {
                "Remove the selected curve from the plot; experiment data is unchanged."
            }
            "Undo" => {
                "Undo the latest reversible operation for the selected experiment. Shortcut: Ctrl/Cmd+Z."
            }
            "calc relative mass" => {
                "Uses Input value as the positive baseline end time and the active mass source. Adds a relative-mass column; optional output name is in the text field."
            }
            "calc conversion" => {
                "Uses Input value as the positive baseline end time and the active mass source. Adds a conversion column; optional output name is in the text field."
            }
            "from mV to mg" => {
                "Applies configured voltage-to-mass calibration. By default replaces the mass column; enter another output name to create a separate calibrated column."
            }
            "from s to h" => {
                "Converts the bound time column from seconds to hours in place; no new column."
            }
            "from C° to K" => {
                "Converts the bound temperature column from Celsius to Kelvin in place; no new column."
            }
            "average on column" => {
                "Reports the mean of the selected curve's Y column over the selected range (or full column); no data changes."
            }
            "conversion rate" => {
                "Adds d(conversion)/dt using the active mass source; optional output name is in the text field. Source is preserved."
            }
            "dim-less mass rate" => {
                "Adds d(relative mass)/dt using the active mass source; optional output name is in the text field. Source is preserved."
            }
            "mass rate" => {
                "Adds d(mass)/dt using the active mass source; optional output name is in the text field. Source is preserved."
            }
            "T rate" => {
                "Adds d(temperature)/dt using the active temperature source; optional output name is in the text field. Source is preserved."
            }
            "Add" => {
                "Requires X or Y selection and Input value. Adds to that existing column in place; no new column."
            }
            "Sub" => {
                "Requires X or Y selection and Input value. Subtracts from that existing column in place; no new column."
            }
            "Mul" => {
                "Requires X or Y selection and Input value. Multiplies that existing column in place; no new column."
            }
            "Divide" => {
                "Requires X or Y selection and a non-zero Input value. Divides that existing column in place; no new column."
            }
            "ln" => {
                "Requires X or Y selection. Replaces that existing column with its natural logarithm; no new column."
            }
            "exp" => {
                "Requires X or Y selection. Replaces that existing column with its exponential; no new column."
            }
            _ => "Run the selected action.",
        }
    }

    pub fn quick_action_panel(
        ui: &mut egui::Ui,
        model: &mut PlotModel,
        state: &mut QuickActionPanelState,
        manage_plot_dialog: &mut NewExperimentDialogState,
    ) -> Result<(), TGAGUIError> {
        state.golden_pipeline_dialog.show(ui.ctx(), model);
        let labels = [
            "Manage Plots",
            "Column Manager",
            "Undo",
            "Clear Selected",
            "Reset View",
            "Refresh",
            "cut before time",
            "cut selected",
            "Zoom To Selection",
            "move time to zero",
            "Delete Plot",
            "golden pipeline",
            // moved arithmetic buttons to be near input field below
            "calc relative mass",
            "calc conversion",
            "from mV to mg",
            "from s to h",
            "from C° to K",
            "average on column",
            "conversion rate",
            "dim-less mass rate",
            "mass rate",
            "T rate",
        ];

        let button_size = egui::vec2(108.0, 22.0);
        egui::Grid::new("exp_kinetics_quick_action_grid")
            .num_columns(4)
            .spacing([4.0, 4.0])
            .show(ui, |ui| {
                for (index, label) in labels.iter().enumerate() {
                    let button = ui.add_sized(button_size, egui::Button::new(*label));
                    if button
                        .on_hover_text(Self::quick_action_tooltip(label))
                        .clicked()
                    {
                        if let Err(err) = Self::handle_quick_action(
                            label,
                            model,
                            state,
                            manage_plot_dialog,
                            ui.ctx(),
                        ) {
                            Self::set_action_error(model, err);
                        }
                    }
                    if (index + 1) % 4 == 0 {
                        ui.end_row();
                    }
                }
            });

        ui.separator();

        // Controls row: Input value + X/Y selector + new column entry
        ui.horizontal_wrapped(|ui| {
            // Input value
            let input_value_label = ui.label("Input value:");
            let input_value = ui
                .add(egui::DragValue::new(&mut state.input_value).speed(0.1))
                .on_hover_text("Numeric parameter used by time cuts, baseline normalization and arithmetic operations.");
            let _ = input_value.labelled_by(input_value_label.id);

            ui.separator();

            // X/Y selector: sets state.x_or_y
            if ui.button("X").on_hover_text("Select the X column for cut/arithmetic operations.").clicked() {
                state.x_or_y = Some(XY::X);
            }
            if ui.button("Y").on_hover_text("Select the Y column for cut/arithmetic operations.").clicked() {
                state.x_or_y = Some(XY::Y);
            }

            ui.separator();

            // Optional output column name for transforms that materialize a derived column.
            let mut tmp = state.input_string.clone().unwrap_or_default();
            let output_column_label = ui.label("output column:");
            let output_column = ui
                .text_edit_singleline(&mut tmp)
                .on_hover_text("Optional output name. Leave empty for the default; using the source name replaces it, another name creates a new column.");
            let _ = output_column.clone().labelled_by(output_column_label.id);
            if output_column.changed() {
                if tmp.trim().is_empty() {
                    state.input_string = None;
                } else {
                    state.input_string = Some(tmp);
                }
            }
        });

        ui.add_space(4.0);
        if let Ok(id) = model.get_experiment_by_selected_curve() {
            if let Ok(columns) = model.list_of_columns(&id) {
                ui.horizontal_wrapped(|ui| {
                    let mass_source_label = ui.label("Mass source:");
                    let mut mass_source = model
                        .active_mass_source_for_experiment(&id)
                        .unwrap_or_else(|_| String::new());
                    let mass_source_combo =
                        egui::ComboBox::from_id_salt("exp_kinetics_mass_source")
                            .selected_text(if mass_source.is_empty() {
                                "choose".to_string()
                            } else {
                                mass_source.clone()
                            })
                            .show_ui(ui, |ui| {
                                for col in &columns {
                                    ui.selectable_value(&mut mass_source, col.clone(), col);
                                }
                            });
                    let _ = mass_source_combo.response.labelled_by(mass_source_label.id);
                    if !mass_source.is_empty()
                        && model
                            .active_mass_sources
                            .get(&id)
                            .map(|current| current != &mass_source)
                            .unwrap_or(true)
                    {
                        model.set_active_mass_source_for_experiment(&id, &mass_source);
                        ui.ctx().request_repaint();
                    }

                    if model.series.get_temperature_col(&id).is_ok() {
                        ui.separator();

                        let temperature_source_label = ui.label("Temperature source:");
                        let mut temperature_source = model
                            .active_temperature_source_for_experiment(&id)
                            .unwrap_or_else(|_| String::new());
                        let temperature_source_combo =
                            egui::ComboBox::from_id_salt("exp_kinetics_temperature_source")
                                .selected_text(if temperature_source.is_empty() {
                                    "choose".to_string()
                                } else {
                                    temperature_source.clone()
                                })
                                .show_ui(ui, |ui| {
                                    for col in &columns {
                                        ui.selectable_value(
                                            &mut temperature_source,
                                            col.clone(),
                                            col,
                                        );
                                    }
                                });
                        let _ = temperature_source_combo
                            .response
                            .labelled_by(temperature_source_label.id);
                        if !temperature_source.is_empty()
                            && model
                                .active_temperature_sources
                                .get(&id)
                                .map(|current| current != &temperature_source)
                                .unwrap_or(true)
                        {
                            model.set_active_temperature_source_for_experiment(
                                &id,
                                &temperature_source,
                            );
                            ui.ctx().request_repaint();
                        }
                    } else if let Ok(experiment) = model.series.get_experiment_by_id(&id) {
                        if let Some(temperature) = experiment.meta.isothermal_temperature {
                            ui.separator();
                            ui.label(format!("Isothermal setpoint: {temperature:.3} K"));
                        }
                    }
                });
            }
        }

        // Arithmetic operation buttons below the optional output column field.
        egui::Grid::new("exp_kinetics_arith_grid")
            .num_columns(3)
            .spacing([4.0, 4.0])
            .show(ui, |ui| {
                let button_size_small = egui::vec2(64.0, 22.0);
                for (idx, label) in ["Divide", "Sub", "Mul", "Add", "ln", "exp"]
                    .iter()
                    .enumerate()
                {
                    let button = ui.add_sized(button_size_small, egui::Button::new(*label));
                    if button
                        .on_hover_text(Self::quick_action_tooltip(label))
                        .clicked()
                    {
                        if let Err(err) = Self::handle_quick_action(
                            label,
                            model,
                            state,
                            manage_plot_dialog,
                            ui.ctx(),
                        ) {
                            Self::set_action_error(model, err);
                        }
                    }
                    if (idx + 1) % 3 == 0 {
                        ui.end_row();
                    }
                }
            });

        ui.checkbox(
            &mut state.apply_only_to_selected_chart,
            "apply only to the selected chart.",
        )
        .on_hover_text("When enabled, Add/Sub/Mul/Divide affect only the selected column inside the selected plot's X range. ln and exp always transform the entire selected column; without the toggle, arithmetic also affects the entire column.");
        ui.separator();
        Ok(())
    }

    fn handle_quick_action(
        action: &str,
        model: &mut PlotModel,
        state: &mut QuickActionPanelState,
        manage_plot_dialog: &mut NewExperimentDialogState,
        ctx: &egui::Context,
    ) -> Result<(), TGAGUIError> {
        match action {
            "Manage Plots" => {
                manage_plot_dialog.manage_plot_open = true;
                ctx.request_repaint();
            }
            "Column Manager" => {
                manage_plot_dialog.column_manager_open = true;
                ctx.request_repaint();
            }
            "golden pipeline" => {
                state.golden_pipeline_dialog.open();
                ctx.request_repaint();
            }
            "Clear Selected" => {
                model.clear_selection_rect();
                model.clear_selection();
            }
            "Reset View" => {
                model.reset_view();
                ctx.request_repaint();
            }
            "Refresh" => {
                ctx.request_repaint();
            }
            "Zoom To Selection" => {
                model.fit_to_selection();
                ctx.request_repaint();
            }
            "cut before time" => {
                model.cut_before_time_for_selected(state.input_value)?;
                ctx.request_repaint();
            }

            "cut selected" => {
                let xy = state.x_or_y.ok_or_else(|| {
                    TGAGUIError::BindingError("Select X or Y before applying operation".to_string())
                })?;
                model.cut_range_x_or_y_for_selected(xy)?;
                ctx.request_repaint();
            }

            "from C° to K" => {
                model.from_C_to_K_of_selected()?;
                ctx.request_repaint();
            }
            "from s to h" => {
                model.from_s_to_h_of_selected()?;
                ctx.request_repaint();
            }
            "Delete Plot" => {
                model.delete_selected_curve()?;
                ctx.request_repaint();
            }
            "Undo" => {
                TopDropDownMenues::undo_selected(model);
                ctx.request_repaint();
            }
            "calc relative mass" => {
                model
                    .relative_mass_for_selected(state.input_value, state.input_string.as_deref())?;
                ctx.request_repaint();
            }
            "calc conversion" => {
                model.conversion_for_selected(state.input_value, state.input_string.as_deref())?;
                ctx.request_repaint();
            }
            "from mV to mg" => {
                model.calibrate_mass_from_voltage_with_new_optional_column_for_selected(
                    state.input_string.as_deref(),
                )?;
                ctx.request_repaint();
            }
            "move time to zero" => {
                model.move_time_to_zero_for_selected()?;
                ctx.request_repaint();
            }

            "Add" => {
                let column = Self::resolve_selected_column(model, state)?;
                info!("adding to column {}", column);
                if state.apply_only_to_selected_chart {
                    model.add_column_in_its_range_for_selected(&column, state.input_value)?;
                } else {
                    model.add_column_for_selected(&column, state.input_value)?;
                }
                ctx.request_repaint();
            }
            "Sub" => {
                let column = Self::resolve_selected_column(model, state)?;
                if state.apply_only_to_selected_chart {
                    model.sub_column_in_its_range_for_selected(&column, state.input_value)?;
                } else {
                    model.sub_column_for_selected(&column, state.input_value)?;
                }
                ctx.request_repaint();
            }
            "Mul" => {
                let column = Self::resolve_selected_column(model, state)?;
                if state.apply_only_to_selected_chart {
                    model.mul_column_in_its_range_for_selected(&column, state.input_value)?;
                } else {
                    model.mul_column_of_selected(&column, state.input_value)?;
                }
                ctx.request_repaint();
            }
            "Divide" => {
                let column = Self::resolve_selected_column(model, state)?;
                if state.apply_only_to_selected_chart {
                    model.div_column_in_its_range_for_selected(&column, state.input_value)?;
                } else {
                    model.div_column_of_selected(&column, state.input_value)?;
                }
                ctx.request_repaint();
            }
            "ln" => {
                let column = Self::resolve_selected_column(model, state)?;
                model.ln_column_for_selected(&column)?;
                ctx.request_repaint();
            }
            "exp" => {
                let column = Self::resolve_selected_column(model, state)?;
                model.exp_column_for_selected(&column)?;
                ctx.request_repaint();
            }

            "average on column" => {
                let column = Self::resolve_selected_column(model, state)?;
                let r = if state.apply_only_to_selected_chart {
                    let id = model.get_experiment_by_selected_curve()?;
                    let range = model.interaction.selection_rect.ok_or_else(|| {
                        TGAGUIError::BindingError(
                            "Select a plot region before averaging a range.".into(),
                        )
                    })?;
                    let bounds = range.bounds();
                    let x_col = model.selected_curve_axis_name(XY::X)?;
                    model.mean_on_interval(&id, &column, &x_col, bounds.0, bounds.1)?
                } else {
                    let id = model.get_experiment_by_selected_curve()?;
                    model.mean_on_column(&id, &column)?
                };
                let message = format!("Mean of '{column}' = {r:.8}");
                info!("{message}");
                let _ = model.push_message(&message);
            }

            "conversion rate" => {
                model.derive_deta_dt_for_selected(state.input_string.as_deref())?;
                ctx.request_repaint();
            }
            "dim-less mass rate" => {
                model.derive_dalpha_dt_for_selected(state.input_string.as_deref())?;
                ctx.request_repaint();
            }
            "mass rate" => {
                model.derive_mass_rate_for_selected(state.input_string.as_deref())?;
                ctx.request_repaint();
            }
            "T rate" => {
                model.derive_temperature_rate_for_selected(state.input_string.as_deref())?;
                ctx.request_repaint();
            }
            _ => {
                return Err(TGAGUIError::BindingError(format!(
                    "Unsupported quick action: {action}"
                )));
            }
        }
        Ok(())
    }

    fn set_action_error(model: &mut PlotModel, err: TGAGUIError) {
        let _ = model.push_message(&err.to_string());
    }
    /// The functions for processing button presses for div, mul, add, sub, ln, exp.
    /// Pressing the "X" or "Y" button calls the this_is_x_or_y function, which specifies the name of the column
    /// to convert. Then, for a specific operation, like add_column_for_selected and sub_column_for_selected
    ///  functions are called, and the value for addition, multiplication, etc. is taken from the input value field.
    fn resolve_selected_column(
        model: &PlotModel,
        state: &QuickActionPanelState,
    ) -> Result<String, TGAGUIError> {
        let xy = state.x_or_y.ok_or_else(|| {
            TGAGUIError::BindingError("Select X or Y before applying operation".to_string())
        })?;
        model.this_is_x_or_y(xy)
    }

    //==================================================================================
    //     COLOUR AND SHOW PANEL
    //==================================================================================
    /// Обрабатывает взаимодействия с пользовательским интерфейсом (кнопки, слайдеры)
    ///
    /// Эта функция отвечает за:
    /// - Обработку нажатий кнопок (сброс вида, очистка выделения)
    /// - Обработку элементов управления параметрами графиков (амплитуда, частота, фаза, цвет)
    ///
    /// # Параметры
    /// * `ui` - изменяемая ссылка на пользовательский интерфейс egui
    /// * `model` - изменяемая ссылка на модель данных графика
    pub fn handle_ui_interactions(
        ui: &mut egui::Ui,
        model: &mut PlotModel,
        _manage_plot_dialog: &mut NewExperimentDialogState,
    ) {
        model.prune_stale_plots();
        ui.separator();

        for curve in model.plots.iter_mut() {
            ui.collapsing(format!("Configure {}", curve.get_label()), |ui| {
                ui.checkbox(&mut curve.shown, "Show");
                ui.horizontal(|ui| {
                    ui.label("Color:");
                    let mut color_rgb = [
                        curve.color[0] as f32 / 255.0,
                        curve.color[1] as f32 / 255.0,
                        curve.color[2] as f32 / 255.0,
                    ];
                    if ui.color_edit_button_rgb(&mut color_rgb).changed() {
                        curve.color = [
                            (color_rgb[0] * 255.0) as u8,
                            (color_rgb[1] * 255.0) as u8,
                            (color_rgb[2] * 255.0) as u8,
                        ];
                    }
                });
            });
        }
    }
}

//========================================================================================
// NEW EXPERIMENT DIALOGUE WINNDOW
//============================================================================================

#[derive(Copy, Clone, Debug, Eq, PartialEq)]
enum SaveSeriesFormat {
    Csv,
    Txt,
    Xlsx,
}

impl SaveSeriesFormat {
    fn extension(self) -> &'static str {
        match self {
            Self::Csv => "csv",
            Self::Txt => "txt",
            Self::Xlsx => "xlsx",
        }
    }

    fn dialog_title(self) -> &'static str {
        match self {
            Self::Csv => "Export Typed CSV",
            Self::Txt => "Export Plain TXT",
            Self::Xlsx => "Export Excel Workbook",
        }
    }
}

fn excel_sheet_name(experiment_id: &str, used_names: &mut HashSet<String>) -> String {
    let sanitized = experiment_id
        .chars()
        .filter(|character| {
            !matches!(character, ':' | '\\' | '/' | '?' | '*' | '[' | ']')
                && !character.is_control()
        })
        .collect::<String>();
    let base = sanitized.trim_matches('\'').trim();
    let base = if base.is_empty() { "Experiment" } else { base };

    for suffix_index in 1.. {
        let suffix = if suffix_index == 1 {
            String::new()
        } else {
            format!(" ({suffix_index})")
        };
        let max_base_chars = 31usize.saturating_sub(suffix.chars().count());
        let prefix = base.chars().take(max_base_chars).collect::<String>();
        let candidate = format!("{prefix}{suffix}");
        if used_names.insert(candidate.to_lowercase()) {
            return candidate;
        }
    }
    unreachable!("the suffix index always produces a unique worksheet name")
}

fn write_xlsx_export(
    model: &PlotModel,
    path: &Path,
    experiment_ids: &[String],
    columns_by_experiment: &HashMap<String, Vec<String>>,
) -> Result<(), String> {
    let mut workbook = rust_xlsxwriter::Workbook::new();
    let mut used_sheet_names = HashSet::with_capacity(experiment_ids.len());

    for experiment_id in experiment_ids {
        let columns = columns_by_experiment
            .get(experiment_id)
            .ok_or_else(|| format!("No columns are selected for experiment '{experiment_id}'"))?;
        let experiment = model
            .series
            .get_experiment_by_id(experiment_id)
            .map_err(|error| error.to_string())?;
        let expressions = columns
            .iter()
            .map(|name| col(name.as_str()))
            .collect::<Vec<_>>();
        let frame = experiment
            .dataset
            .frame
            .clone()
            .select(expressions)
            .collect()
            .map_err(|error| format!("Could not read experiment '{experiment_id}': {error}"))?;
        let numeric_columns = columns
            .iter()
            .map(|name| {
                frame
                    .column(name)
                    .and_then(|column| column.f64())
                    .map(|values| (name.as_str(), values))
                    .map_err(|error| {
                        format!("Column '{name}' in experiment '{experiment_id}' is not numeric: {error}")
                    })
            })
            .collect::<Result<Vec<_>, _>>()?;
        if numeric_columns.len() > 16_384 {
            return Err(format!(
                "Experiment '{experiment_id}' has more than Excel's 16,384-column worksheet limit"
            ));
        }
        if frame.height() > 1_048_575 {
            return Err(format!(
                "Experiment '{experiment_id}' has more rows than Excel's worksheet limit allows with a header"
            ));
        }

        let worksheet = workbook.add_worksheet();
        let sheet_name = excel_sheet_name(experiment_id, &mut used_sheet_names);
        worksheet
            .set_name(&sheet_name)
            .map_err(|error| format!("Could not name worksheet '{sheet_name}': {error}"))?;
        for (column_index, (name, _)) in numeric_columns.iter().enumerate() {
            let column_index = u16::try_from(column_index).map_err(|_| {
                "Excel export supports at most 16,384 columns per experiment".to_string()
            })?;
            worksheet
                .write_string(0, column_index, *name)
                .map_err(|error| error.to_string())?;
        }

        for row_index in 0..frame.height() {
            let row = u32::try_from(row_index + 1)
                .map_err(|_| "Excel export exceeds the worksheet row limit".to_string())?;
            for (column_index, (_, values)) in numeric_columns.iter().enumerate() {
                let column_index = u16::try_from(column_index).map_err(|_| {
                    "Excel export supports at most 16,384 columns per experiment".to_string()
                })?;
                if let Some(value) = values.get(row_index) {
                    if value.is_finite() {
                        worksheet
                            .write_number(row, column_index, value)
                            .map_err(|error| error.to_string())?;
                    } else {
                        let label = if value.is_nan() {
                            "NaN".to_string()
                        } else if value.is_sign_positive() {
                            "Infinity".to_string()
                        } else {
                            "-Infinity".to_string()
                        };
                        worksheet
                            .write_string(row, column_index, label)
                            .map_err(|error| error.to_string())?;
                    }
                }
            }
        }
    }

    workbook
        .save(path)
        .map_err(|error| format!("Could not write Excel workbook: {error}"))
}

fn column_nature_label(nature: ColumnNature) -> &'static str {
    match nature {
        ColumnNature::Time => "Time",
        ColumnNature::Temperature => "Temperature",
        ColumnNature::Mass => "Mass",
        ColumnNature::Concentration => "Concentration",
        ColumnNature::Conversion => "Conversion",
        ColumnNature::DimensionlessMass => "Relative mass",
        ColumnNature::MassRate => "Mass rate",
        ColumnNature::ConversionRate => "Conversion rate",
        ColumnNature::TemperatureRate => "Temperature rate",
        ColumnNature::DimensionlessMassRate => "Relative-mass rate",
        ColumnNature::Unknown => "Unclassified",
        ColumnNature::ActivationEnergy => "Activation energy",
        ColumnNature::PredexFactor => "Pre-exponential factor",
        ColumnNature::R2 => "R2",
        ColumnNature::FittedCurve => "Fitted curve",
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
enum SaveScope {
    #[default]
    Series,
    SelectedExperiments,
    Experiment,
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
enum SaveDialogPurpose {
    Document,
    #[default]
    Export,
}

#[derive(Clone, Debug)]
struct PendingFileWrite {
    path: PathBuf,
    purpose: SaveDialogPurpose,
    export: Option<ExportRequest>,
}

#[derive(Clone, Debug)]
struct ExportRequest {
    format: SaveSeriesFormat,
    experiment_ids: Vec<String>,
    columns_by_experiment: HashMap<String, Vec<String>>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum ImportFileFormat {
    Raw,
    TypedDataset,
    TypedSeries,
}

fn detect_import_file_format(path: &Path) -> Result<ImportFileFormat, String> {
    let file = File::open(path).map_err(|error| format!("Could not read file: {error}"))?;
    let lines = BufReader::new(file)
        .lines()
        .take(16)
        .collect::<Result<Vec<_>, _>>()
        .map_err(|error| format!("Could not read file header: {error}"))?;
    let first = lines
        .iter()
        .find(|line| !line.trim().is_empty())
        .map(String::as_str)
        .unwrap_or_default();

    if first.starts_with("# format_version=")
        || lines.iter().any(|line| line.starts_with("#META\t"))
    {
        Ok(ImportFileFormat::TypedSeries)
    } else if first == "# KiThe TGA Dataset" {
        Ok(ImportFileFormat::TypedDataset)
    } else {
        Ok(ImportFileFormat::Raw)
    }
}

/// State for the new experiment file selection dialog
#[derive(Clone, Debug, Default)]
pub struct PlotConfig {
    pub experiment: Option<String>,
    pub x_column: Option<String>,
    pub y_column: Option<String>,
}

#[derive(Clone, Debug)]
pub struct NewExperimentDialogState {
    pub open: bool,
    pub manage_plot_open: bool,
    pub column_manager_open: bool,
    pub save_series_open: bool,
    save_dialog_purpose: SaveDialogPurpose,
    pub current_dir: PathBuf,
    pub selected_file: Option<PathBuf>,
    file_path_input: String,
    pub last_error: Option<String>,
    save_series_format: SaveSeriesFormat,
    save_series_filename_input: String,
    save_series_last_status: Option<String>,
    save_scope: SaveScope,
    save_experiment_id: Option<String>,
    save_selected_experiments: HashSet<String>,
    save_selected_columns: HashMap<String, HashSet<String>>,
    current_document_path: Option<PathBuf>,
    document_generation: u64,
    saved_data_revision: u64,
    pending_file_write: Option<PendingFileWrite>,
    pending_open_series: Option<PathBuf>,
    unsaved_changes_open: bool,
    heating_rate_input: String,
    isothermal_temp_input: String,
    thermal_mode: ExperimentThermalMode,
    use_measured_temperature_mean: bool,
    bind_mass_input: String,
    bind_mass_unit: Unit,
    bind_temperature_input: String,
    bind_temperature_unit: Unit,
    bind_time_input: String,
    bind_time_unit: Unit,
    meta: ExperimentMeta,
    plot_configs: Vec<PlotConfig>,
    append_mode: bool,
    append_target_id: Option<String>,
    append_replace_target: bool,
    append_result_id: String,
    append_preview: Option<TGADataset>,
    append_preview_rows: Option<usize>,
    append_preview_signature: Option<String>,
    import_file_format: Option<ImportFileFormat>,
    replace_series_on_import: bool,
    #[cfg(test)]
    pub(crate) test_next_open_file: Option<PathBuf>,
    #[cfg(test)]
    pub(crate) test_next_save_file: Option<PathBuf>,
}

impl NewExperimentDialogState {
    pub(crate) fn close_window_on_layer(&mut self, layer: Option<egui::LayerId>) -> bool {
        let Some(layer) = layer else { return false };
        if self.append_mode
            && layer.id == egui::Id::unique(Some("Append to Experiment"))
            && self.open
        {
            self.open = false;
            return true;
        }
        let windows = [
            ("Import Data", &mut self.open),
            ("Manage Plots", &mut self.manage_plot_open),
            ("Column Manager", &mut self.column_manager_open),
        ];
        for (title, is_open) in windows {
            if layer.id == egui::Id::unique(Some(title)) && *is_open {
                *is_open = false;
                return true;
            }
        }
        if self.save_series_open
            && [
                "Export Typed CSV",
                "Export Plain TXT",
                "Export Excel Workbook",
            ]
            .into_iter()
            .any(|title| layer.id == egui::Id::unique(Some(title)))
        {
            self.save_series_open = false;
            return true;
        }
        false
    }

    /// Create a new dialog state with the current directory
    pub fn new() -> Self {
        Self {
            open: false,
            manage_plot_open: false,
            column_manager_open: false,
            save_series_open: false,
            save_dialog_purpose: SaveDialogPurpose::Export,
            current_dir: std::env::current_dir().unwrap_or_else(|_| PathBuf::from(".")),
            selected_file: None,
            file_path_input: String::new(),
            last_error: None,
            save_series_format: SaveSeriesFormat::Csv,
            save_series_filename_input: "series_export.csv".to_string(),
            save_series_last_status: None,
            save_scope: SaveScope::Series,
            save_experiment_id: None,
            save_selected_experiments: HashSet::new(),
            save_selected_columns: HashMap::new(),
            current_document_path: None,
            document_generation: 0,
            saved_data_revision: 0,
            pending_file_write: None,
            pending_open_series: None,
            unsaved_changes_open: false,
            heating_rate_input: String::new(),
            isothermal_temp_input: String::new(),
            thermal_mode: ExperimentThermalMode::Unspecified,
            use_measured_temperature_mean: false,
            bind_mass_input: String::new(),
            bind_mass_unit: Unit::MilliVolt,
            bind_temperature_input: String::new(),
            bind_temperature_unit: Unit::Celsius,
            bind_time_input: String::new(),
            bind_time_unit: Unit::Second,
            meta: ExperimentMeta::new(),
            plot_configs: vec![PlotConfig::default()],
            append_mode: false,
            append_target_id: None,
            append_replace_target: false,
            append_result_id: String::new(),
            append_preview: None,
            append_preview_rows: None,
            append_preview_signature: None,
            import_file_format: None,
            replace_series_on_import: false,
            #[cfg(test)]
            test_next_open_file: None,
            #[cfg(test)]
            test_next_save_file: None,
        }
    }

    fn select_thermal_mode(&mut self, mode: ExperimentThermalMode) {
        self.thermal_mode = mode;
        if mode != ExperimentThermalMode::Isothermal {
            self.use_measured_temperature_mean = false;
            self.isothermal_temp_input.clear();
        }
        if mode != ExperimentThermalMode::NonIsothermal {
            self.heating_rate_input.clear();
        }
    }

    fn time_unit_options() -> &'static [Unit] {
        &[Unit::Second, Unit::Minute, Unit::Hour]
    }

    fn temperature_unit_options() -> &'static [Unit] {
        &[Unit::Kelvin, Unit::Celsius]
    }

    fn mass_unit_options() -> &'static [Unit] {
        &[Unit::MilliVolt, Unit::Milligram, Unit::Gram]
    }

    fn show_bind_input_block(
        ui: &mut egui::Ui,
        label: &str,
        input: &mut String,
        unit: &mut Unit,
        combo_id: &str,
        options: &[Unit],
    ) {
        ui.vertical(|ui| {
            let input_label = ui.label(label);
            let input = ui.text_edit_singleline(input);
            let _ = input.labelled_by(input_label.id);
            let unit_label = ui.label(format!("{label} unit"));
            let combo = egui::ComboBox::from_id_salt(combo_id)
                .selected_text(unit.to_string())
                .show_ui(ui, |ui| {
                    for option in options {
                        ui.selectable_value(unit, *option, option.to_string());
                    }
                });
            let _ = combo.response.labelled_by(unit_label.id);
        });
    }

    fn first_matching_column(columns: &[String], candidates: &[&str]) -> Option<String> {
        candidates.iter().find_map(|candidate| {
            columns
                .iter()
                .find(|col| col.as_str() == *candidate)
                .cloned()
        })
    }

    /// 1) label "bind m" - under it an entry field for text input and under it a drop down menu field with Unit enum one_experiment_dataset.rs
    ///  with Unit::mV as default 2) a label "bind T" - under it an entry field for text input and under it a drop down menu field with Unit
    /// enum with Unit::Celsius as default 3) a label "bind t" - under it an entry field for text input and under it a drop down menu field with
    /// Unit enum with Unit::Second as default.
    fn show_bindings_row(&mut self, ui: &mut egui::Ui) {
        ui.horizontal(|ui| {
            Self::show_bind_input_block(
                ui,
                "bind m",
                &mut self.bind_mass_input,
                &mut self.bind_mass_unit,
                "bind_m_unit",
                Self::mass_unit_options(),
            );

            Self::show_bind_input_block(
                ui,
                "bind T",
                &mut self.bind_temperature_input,
                &mut self.bind_temperature_unit,
                "bind_t_unit",
                Self::temperature_unit_options(),
            );
            Self::show_bind_input_block(
                ui,
                "bind t",
                &mut self.bind_time_input,
                &mut self.bind_time_unit,
                "bind_time_unit",
                Self::time_unit_options(),
            );
        });
    }

    fn autofill_bind_inputs_from_last_experiment(&mut self, model: &PlotModel) {
        let Some(last_exp) = model.series.experiments().last() else {
            return;
        };
        let Ok(columns) = model.list_of_columns(&last_exp.meta.id) else {
            return;
        };

        if self.bind_time_input.trim().is_empty() {
            if let Some(col) = Self::first_matching_column(&columns, &["t", "time"]) {
                self.bind_time_input = col;
            }
        }
        if self.bind_temperature_input.trim().is_empty() {
            if let Some(col) = Self::first_matching_column(&columns, &["T", "temperature"]) {
                self.bind_temperature_input = col;
            }
        }
        if self.bind_mass_input.trim().is_empty() {
            if let Some(col) = Self::first_matching_column(&columns, &["m", "mass"]) {
                self.bind_mass_input = col;
            }
        }
    }

    fn apply_bindings_to_last_experiment(
        &mut self,
        model: &mut PlotModel,
    ) -> Result<(), TGAGUIError> {
        let Some(last_exp) = model.series.experiments().last() else {
            return Ok(());
        };
        let exp_id = last_exp.meta.id.clone();

        let mass = self.bind_mass_input.trim().to_owned();
        let temperature = self.bind_temperature_input.trim().to_owned();
        let time = self.bind_time_input.trim().to_owned();
        let mass_unit = self.bind_mass_unit;
        let temperature_unit = self.bind_temperature_unit;
        let time_unit = self.bind_time_unit;

        model
            .series
            .try_mutate_by_id(&exp_id, |exp| {
                if !mass.is_empty() {
                    exp.dataset.bind_mass_inplace(&mass, mass_unit)?;
                }
                if !temperature.is_empty() {
                    exp.dataset
                        .bind_temperature_inplace(&temperature, temperature_unit)?;
                }
                if !time.is_empty() {
                    exp.dataset.bind_time_inplace(&time, time_unit)?;
                }
                Ok(())
            })
            .map_err(|e| TGAGUIError::BindingError(format!("binding failed: {:?}", e)))?;

        if !mass.is_empty() {
            info!(
                "bound mass column {}",
                model
                    .series
                    .get_mass_col(&exp_id)
                    .map_err(|e| TGAGUIError::BindingError(format!(
                        "failed to read bound mass: {:?}",
                        e
                    )))?
            );
        }
        if !temperature.is_empty() {
            info!(
                "bound temperature column {}",
                model.series.get_temperature_col(&exp_id).map_err(
                    |e| TGAGUIError::BindingError(format!(
                        "failed to read bound temperature: {:?}",
                        e
                    ))
                )?
            );
        }
        if !time.is_empty() {
            info!(
                "bound time column {}",
                model
                    .series
                    .get_time_col(&exp_id)
                    .map_err(|e| TGAGUIError::BindingError(format!(
                        "failed to read bound time: {:?}",
                        e
                    )))?
            );
        }

        Ok(())
    }

    fn finish_new_experiment_import(
        &mut self,
        model: &mut PlotModel,
        exp_id: &str,
    ) -> Result<(), TGAGUIError> {
        self.finish_new_experiment_import_with_bindings(model, exp_id, true)
    }

    fn finish_new_experiment_import_with_bindings(
        &mut self,
        model: &mut PlotModel,
        exp_id: &str,
        apply_user_bindings: bool,
    ) -> Result<(), TGAGUIError> {
        let mut adopted_sensor_temperature = None;
        let result = (|| {
            model.series.try_mutate_by_id(exp_id, |exp| {
                exp.check_nulls_for_operation_borrowed("parsing file")?;
                exp.sample_all_columns_even_layers_table(50, Vec::new())
            })?;
            if apply_user_bindings {
                self.autofill_bind_inputs_from_last_experiment(model);
                self.apply_bindings_to_last_experiment(model)?;
            }
            if self.thermal_mode == ExperimentThermalMode::Isothermal
                && self.use_measured_temperature_mean
            {
                let column = self.bind_temperature_input.trim();
                if column.is_empty() {
                    return Err(TGAGUIError::BindingError(
                        "Choose and bind a measured temperature column before using its mean"
                            .into(),
                    ));
                }
                let dataset = &model.series.get_experiment_by_id(exp_id)?.dataset;
                let preview = dataset
                    .preview_rows(&[column], usize::MAX)
                    .map_err(TGADomainError::from)?;
                let values = preview.columns.get(column).ok_or_else(|| {
                    TGAGUIError::BindingError(format!(
                        "Temperature column '{column}' is unavailable"
                    ))
                })?;
                let mut valid_count = 0usize;
                let mut measured_mean = 0.0_f64;
                for value in values
                    .iter()
                    .flatten()
                    .copied()
                    .filter(|value| value.is_finite())
                {
                    valid_count += 1;
                    measured_mean += (value - measured_mean) / valid_count as f64;
                }
                if valid_count == 0 {
                    return Err(TGAGUIError::BindingError(format!(
                        "Temperature column '{column}' has no finite measurements"
                    )));
                }
                let mean_kelvin = match self.bind_temperature_unit {
                    Unit::Kelvin => measured_mean,
                    Unit::Celsius => measured_mean + 273.15,
                    unit => {
                        return Err(TGAGUIError::BindingError(format!(
                            "Cannot use {:?} as a temperature unit",
                            unit
                        )));
                    }
                };
                if !mean_kelvin.is_finite() || mean_kelvin <= 0.0 {
                    return Err(TGAGUIError::BindingError(
                        "Mean measured temperature must be finite and above 0 K".into(),
                    ));
                }
                model.series.try_mutate_by_id(exp_id, |experiment| {
                    experiment.meta.isothermal_temperature = Some(mean_kelvin);
                    experiment.meta.thermal_mode = ExperimentThermalMode::Isothermal;
                    Ok(())
                })?;
                adopted_sensor_temperature = Some((valid_count, mean_kelvin));
            }
            Ok(())
        })();

        if let Err(operation_error) = result {
            return match model.drop_experiment_by_id(exp_id) {
                Ok(()) => Err(operation_error),
                Err(cleanup_error) => Err(TGAGUIError::BindingError(format!(
                    "Import failed: {operation_error:?}; rollback of experiment '{exp_id}' also failed: {cleanup_error:?}"
                ))),
            };
        }

        if let Err(error) = self.create_default_plots_for_last_experiment(model) {
            model.plots.retain(|plot| plot.experiment_id != exp_id);
            model.reset_view();
            return Err(TGAGUIError::post_commit_warning(format!(
                "experiment '{exp_id}' was imported and remains available, but its plots could not be created ({error})"
            )));
        }
        if let Some((valid_count, mean_kelvin)) = adopted_sensor_temperature {
            let _ = model.push_message(&format!(
                "Isothermal temperature set from {valid_count} valid measurements: {mean_kelvin:.3} K"
            ));
        }
        Ok(())
    }

    fn create_default_plots_for_last_experiment(
        &self,
        model: &mut PlotModel,
    ) -> Result<(), TGAGUIError> {
        let Some(last_exp) = model.series.experiments().last() else {
            return Ok(());
        };
        let exp_id = last_exp.meta.id.clone();
        self.create_default_plots_for_experiment(model, &exp_id)
    }

    fn create_default_plots_for_experiment(
        &self,
        model: &mut PlotModel,
        exp_id: &str,
    ) -> Result<(), TGAGUIError> {
        let time_col = model.series.get_time_col(&exp_id).or_else(|_| {
            model
                .list_of_columns(exp_id)
                .ok()
                .and_then(|columns| Self::first_matching_column(&columns, &["t", "time"]))
                .ok_or_else(|| {
                    TGAGUIError::BindingError(format!(
                        "no time column is bound or inferable for '{}'",
                        exp_id
                    ))
                })
        })?;

        let mut created_any = false;
        for y_col in [
            model.series.get_temperature_col(exp_id).ok(),
            model.series.get_mass_col(exp_id).ok(),
        ]
        .into_iter()
        .flatten()
        {
            if y_col == time_col {
                continue;
            }
            model.set_x(exp_id, &time_col)?;
            model.set_y(exp_id, &y_col)?;
            model.create_points_for_curve(exp_id)?;
            created_any = true;
        }

        if created_any {
            model.reset_view();
        }

        Ok(())
    }

    fn import_typed_series_file(
        &mut self,
        model: &mut PlotModel,
        path: &Path,
    ) -> Result<Vec<String>, TGAGUIError> {
        let imported = TGASeries::from_csv_series(path)?;
        let imported_ids = imported.ids();
        if imported_ids.is_empty() {
            return Err(TGAGUIError::BindingError(
                "Typed CSV contains no experiments".into(),
            ));
        }

        if self.replace_series_on_import {
            model.series = imported;
            self.document_generation = self.document_generation.wrapping_add(1);
            self.pending_file_write = None;
            self.append_preview = None;
            self.append_preview_rows = None;
            self.append_preview_signature = None;
            model.plots.clear();
            model.active_mass_sources.clear();
            model.active_temperature_sources.clear();
            self.current_document_path = Some(path.to_path_buf());
            self.saved_data_revision = model.series.data_revision();
        } else {
            let mut candidate = model.series.clone();
            for experiment in imported.experiments().iter().cloned() {
                candidate.push(experiment)?;
            }
            model.series = candidate;
        }

        let mut plot_warnings = Vec::new();
        for id in &imported_ids {
            if let Err(error) = self.create_default_plots_for_experiment(model, id) {
                model.plots.retain(|plot| plot.experiment_id != *id);
                plot_warnings.push(format!("{id}: {error}"));
            }
        }
        model.reset_view();
        if plot_warnings.is_empty() {
            let _ = model.push_message(&format!(
                "Imported {} experiment(s) from typed KiThe CSV",
                imported_ids.len()
            ));
        } else {
            let _ = model.push_message(&format!(
                "Imported {} experiment(s); some plots could not be created, but the data and metadata remain available: {}",
                imported_ids.len(),
                plot_warnings.join("; ")
            ));
        }
        Ok(imported_ids)
    }

    fn request_open_typed_series(&mut self, model: &mut PlotModel, path: &Path) {
        if self.is_document_dirty(model) {
            self.pending_open_series = Some(path.to_path_buf());
            self.unsaved_changes_open = true;
            self.open = false;
        } else {
            self.load_typed_series_document(model, path);
        }
    }

    fn load_typed_series_document(&mut self, model: &mut PlotModel, path: &Path) {
        self.replace_series_on_import = true;
        match self.import_typed_series_file(model, path) {
            Ok(ids) => {
                info!("Opened typed KiThe document with {} experiments", ids.len());
                self.open = false;
                self.last_error = None;
            }
            Err(error) => {
                self.last_error = Some(error.to_string());
                let _ = model.push_message(&format!("Open document failed: {error}"));
            }
        }
        self.replace_series_on_import = false;
    }

    pub fn show_unsaved_changes_dialog(&mut self, ctx: &egui::Context, model: &mut PlotModel) {
        if !self.unsaved_changes_open {
            return;
        }
        let mut save = false;
        let mut discard = false;
        let mut cancel = false;
        egui::Window::new("Unsaved document changes")
            .collapsible(false)
            .resizable(false)
            .show(ctx, |ui| {
                ui.label("The current document has unsaved changes. Save them before opening another document?");
                ui.horizontal(|ui| {
                    save = ui.button("Save").clicked();
                    discard = ui.button("Discard").clicked();
                    cancel = ui.button("Cancel").clicked();
                });
            });

        if cancel {
            self.pending_open_series = None;
            self.unsaved_changes_open = false;
        } else if discard {
            self.unsaved_changes_open = false;
            if let Some(path) = self.pending_open_series.take() {
                self.load_typed_series_document(model, &path);
            }
        } else if save {
            self.unsaved_changes_open = false;
            if let Some(path) = self.current_document_path.clone() {
                self.write_document(model, &path);
            } else {
                self.open_save_document_dialog(model);
            }
        }
    }

    /// Show the new experiment dialog
    pub fn show_new_experiment_dialogue(&mut self, ctx: &egui::Context, model: &mut PlotModel) {
        let mut should_load_experiment = false;
        let mut commit_append = false;
        let mut choose_another_file = false;

        if self.open {
            let mut dialog_open = self.open;
            egui::Window::new(if self.append_mode { "Append to Experiment" } else { "Import Data" })
                .open(&mut dialog_open)
                .resizable(true)
                .default_size(if self.append_mode { [760.0, 520.0] } else { [700.0, 480.0] })
                .show(ctx, |ui| {
                    ui.horizontal(|ui| {
                        ui.label("Selected file:");
                        ui.label(
                            self.selected_file
                                .as_deref()
                                .map(|path| path.display().to_string())
                                .unwrap_or_else(|| "No file selected".to_string()),
                        );
                        if ui.button("Choose another...").clicked() {
                            choose_another_file = true;
                        }
                    });

                    ui.separator();

                    if self.append_mode && self.append_preview.is_some() {
                        ui.horizontal(|ui| {
                            ui.strong(format!("Preview: {} combined rows", self.append_preview_rows.unwrap_or_default()));
                            let commit_label = if self.append_replace_target { "Replace source with combined experiment" } else { "Add combined experiment" };
                            if ui.button(commit_label)
                                .on_hover_text("Commit the validated preview. Replacement is explicit and happens only after validation.")
                                .clicked()
                            {
                                commit_append = true;
                            }
                        });
                    }

                    // Metadata panel (right side)
                    egui::Panel::right("metadata_panel")
                        .resizable(true)
                        .default_size(250.0)
                        .show(ui, |ui| {
                            ui.heading("Metadata");
                            ui.separator();
                            if self.append_mode {
                                let experiments = model.list_of_experiments();
                                if self.append_target_id.as_ref().is_none_or(|id| !experiments.contains(id)) {
                                    self.append_target_id = experiments.first().cloned();
                                }
                                let target_changed = egui::ComboBox::from_label("Append to")
                                    .selected_text(self.append_target_id.as_deref().unwrap_or("Select experiment"))
                                    .show_ui(ui, |ui| {
                                        for id in &experiments {
                                            ui.selectable_value(&mut self.append_target_id, Some(id.clone()), id);
                                        }
                                    }).response.changed();
                                if target_changed {
                                    self.append_preview = None;
                                    self.append_preview_rows = None;
                                    self.append_preview_signature = None;
                                }
                                ui.checkbox(&mut self.append_replace_target, "Replace the original experiment");
                                if !self.append_replace_target {
                                    let label = ui.label("New experiment ID:");
                                    let input = ui.text_edit_singleline(&mut self.append_result_id);
                                    let _ = input.labelled_by(label.id);
                                }
                                ui.label("Incoming columns must match the selected experiment's column names and order. Time is preserved if already continuous; otherwise the segment is shifted to the next typical time step.");
                                ui.separator();
                            }
                            let footer_height = 64.0 + if self.last_error.is_some() { 40.0 } else { 0.0 };
                            let content_height = (ui.available_height() - footer_height).max(100.0);
                            egui::ScrollArea::vertical()
                                .auto_shrink([false, false])
                                .max_height(content_height)
                                .show(ui, |ui| {
                            ui.heading("📝 Metadata");
                            ui.separator();

                            let label = ui.label("File path:");
                            let input = ui.text_edit_singleline(&mut self.file_path_input);
                            let path_changed = input.changed();
                            let _ = input.labelled_by(label.id);
                            if path_changed && !self.append_mode {
                                let path = PathBuf::from(self.file_path_input.trim());
                                self.update_import_file_format(&path);
                            }

                            if self.import_file_format == Some(ImportFileFormat::TypedSeries)
                                && !self.append_mode
                            {
                                ui.label("KiThe series CSV detected. Experiment metadata, bindings, units, provenance, operation history, and the archived calibration line are read from the file. Archived coefficients do not replace current application Settings.");
                                ui.horizontal(|ui| {
                                    ui.selectable_value(
                                        &mut self.replace_series_on_import,
                                        true,
                                        "Open as document",
                                    );
                                    ui.selectable_value(
                                        &mut self.replace_series_on_import,
                                        false,
                                        "Add to current series",
                                    );
                                });
                            } else if !self.append_mode {
                                if self.import_file_format == Some(ImportFileFormat::TypedDataset) {
                                    ui.label("KiThe dataset CSV detected. Column metadata, units, and bindings will be restored from the file.");
                                }
                                ui.add_space(5.0);
                                let label = ui.label("Experiment ID:");
                                let input = ui.text_edit_singleline(&mut self.meta.id);
                                let _ = input.labelled_by(label.id);

                                ui.add_space(5.0);
                                ui.horizontal(|ui| {
                                    ui.label("Thermal regime:");
                                    if ui.radio_value(&mut self.thermal_mode, ExperimentThermalMode::Unspecified, "Infer from metadata").changed() {
                                        self.select_thermal_mode(ExperimentThermalMode::Unspecified);
                                    }
                                    if ui.radio_value(&mut self.thermal_mode, ExperimentThermalMode::Isothermal, "Isothermal").changed() {
                                        self.select_thermal_mode(ExperimentThermalMode::Isothermal);
                                    }
                                    if ui.radio_value(&mut self.thermal_mode, ExperimentThermalMode::NonIsothermal, "Non-isothermal").changed() {
                                        self.select_thermal_mode(ExperimentThermalMode::NonIsothermal);
                                    }
                                });
                                if self.thermal_mode != ExperimentThermalMode::Isothermal {
                                    let label = ui.label("Heating Rate (K/min):");
                                    let input = ui.text_edit_singleline(&mut self.heating_rate_input);
                                    let _ = input.labelled_by(label.id);
                                }
                                if self.thermal_mode != ExperimentThermalMode::NonIsothermal {
                                    let label = ui.label("Isothermal Temperature (K):");
                                    let input = ui.text_edit_singleline(&mut self.isothermal_temp_input);
                                    let _ = input.labelled_by(label.id);
                                    ui.checkbox(&mut self.use_measured_temperature_mean, "Use mean of the bound measured-temperature column");
                                }

                                ui.add_space(5.0);
                                let label = ui.label("Comment:");
                                let comment = self.meta.comment.get_or_insert_with(String::new);
                                let input = ui.text_edit_multiline(comment);
                                let _ = input.labelled_by(label.id);
                            }

                            ui.add_space(8.0);
                            if !self.append_mode
                                && self.import_file_format == Some(ImportFileFormat::TypedDataset)
                            {
                                ui.label("Bindings come from the file; edit them later in Column Manager if needed.");
                            } else if !self.append_mode
                                && self.import_file_format == Some(ImportFileFormat::TypedSeries)
                            {
                                ui.label("Bindings are restored independently for every imported experiment.");
                            } else {
                                self.show_bindings_row(ui);
                            }
                                });

                            if let Some(err) = &self.last_error {
                                let mut preview = err.chars().take(80).collect::<String>();
                                if err.chars().count() > 80 {
                                    preview.push_str("...");
                                }
                                ui.add_sized(
                                    [ui.available_width(), 36.0],
                                    egui::Label::new(
                                        egui::RichText::new(preview).color(egui::Color32::RED),
                                    ),
                                )
                                .on_hover_text(err);
                                ui.add_space(4.0);
                            } else {
                                ui.add_space(12.0);
                            }
                            let action_label = if self.append_mode {
                                "Preview continuation"
                            } else if self.import_file_format == Some(ImportFileFormat::TypedSeries) {
                                if self.replace_series_on_import {
                                    "Open series"
                                } else {
                                    "Import series into current document"
                                }
                            } else {
                                "Import experiment"
                            };
                            let action = ui.add_sized(
                                [ui.available_width(), 42.0],
                                egui::Button::new(
                                    egui::RichText::new(action_label).size(18.0).strong(),
                                ),
                            );
                            let action_tooltip = if self.append_mode {
                                "Preview aligned rows, then choose whether to append or replace the experiment."
                            } else if self.import_file_format == Some(ImportFileFormat::TypedSeries) {
                                if self.replace_series_on_import {
                                    "Replace the current document with this typed series after the unsaved-changes prompt. Save will target this file."
                                } else {
                                    "Add these experiments to the current document without changing its Save path."
                                }
                            } else if self.import_file_format == Some(ImportFileFormat::TypedDataset) {
                                "Import this typed dataset as a new experiment; column bindings and units are restored from the file."
                            } else {
                                "Import this raw file as a new experiment. Set column bindings and experiment metadata above."
                            };
                            if action
                                .on_hover_text(action_tooltip)
                                .clicked()
                            {
                                should_load_experiment = true;
                            }

                        });
                });
            self.open = dialog_open;

            if choose_another_file {
                let title = if self.append_mode {
                    "Append to Experiment"
                } else {
                    "Import Data"
                };
                if let Some(path) = self.pick_open_file(title, false) {
                    self.set_import_file(path);
                    self.append_preview = None;
                    self.append_preview_rows = None;
                    self.append_preview_signature = None;
                }
            }

            // Resolve the selected format when a path was entered manually.
            if should_load_experiment && !self.append_mode {
                let path = if self.file_path_input.trim().is_empty() {
                    self.selected_file.clone()
                } else {
                    Some(PathBuf::from(self.file_path_input.trim()))
                };
                self.import_file_format = path
                    .as_deref()
                    .and_then(|path| detect_import_file_format(path).ok());
            }

            // Handle import or continuation preview outside the window closure.
            if should_load_experiment {
                if self.append_mode {
                    let path = if self.file_path_input.trim().is_empty() {
                        self.selected_file.clone()
                    } else {
                        Some(PathBuf::from(self.file_path_input.trim()))
                    };
                    let result = (|| -> Result<(), String> {
                        let path =
                            path.ok_or_else(|| "Select a continuation file first".to_string())?;
                        let target_id = self
                            .append_target_id
                            .as_deref()
                            .ok_or_else(|| "Select a target experiment first".to_string())?;
                        let target = model
                            .series
                            .get_experiment_by_id(target_id)
                            .map_err(|e| format!("{e:?}"))?;
                        target.dataset.schema.time.as_deref().ok_or_else(|| {
                            "Target experiment has no bound time column".to_string()
                        })?;
                        let mut continuation =
                            TGADataset::from_csv_raw(&path).map_err(|e| e.to_string())?;
                        let incoming_names = continuation
                            .frame
                            .clone()
                            .collect_schema()
                            .map_err(|e| e.to_string())?;
                        let target_names = target
                            .dataset
                            .frame
                            .clone()
                            .collect_schema()
                            .map_err(|e| e.to_string())?
                            .iter_names()
                            .map(|name| name.to_string())
                            .collect::<Vec<_>>();
                        let mappings = [
                            (
                                target.dataset.schema.time.as_deref(),
                                self.bind_time_input.trim(),
                            ),
                            (
                                target.dataset.schema.temperature.as_deref(),
                                self.bind_temperature_input.trim(),
                            ),
                            (
                                target.dataset.schema.mass.as_deref(),
                                self.bind_mass_input.trim(),
                            ),
                        ];
                        for (target_column, incoming_column) in mappings {
                            if let Some(target_column) = target_column {
                                let target_unit = target
                                    .dataset
                                    .schema
                                    .columns
                                    .get(target_column)
                                    .map(|meta| meta.unit);
                                let incoming_unit = match target_column {
                                    column
                                        if Some(column)
                                            == target.dataset.schema.time.as_deref() =>
                                    {
                                        self.bind_time_unit
                                    }
                                    column
                                        if Some(column)
                                            == target.dataset.schema.temperature.as_deref() =>
                                    {
                                        self.bind_temperature_unit
                                    }
                                    _ => self.bind_mass_unit,
                                };
                                if target_unit != Some(incoming_unit) {
                                    return Err(format!(
                                        "Incoming unit for '{target_column}' must match the target experiment unit {:?}",
                                        target_unit
                                    ));
                                }
                                if !incoming_names
                                    .iter_names()
                                    .any(|name| name.as_str() == incoming_column)
                                {
                                    return Err(format!(
                                        "Incoming file has no bound column '{incoming_column}'"
                                    ));
                                }
                            }
                        }
                        let mut selected_sources = Vec::with_capacity(target_names.len());
                        let mut expressions = Vec::with_capacity(target_names.len());
                        for target_column in &target_names {
                            let source_column = if Some(target_column.as_str())
                                == target.dataset.schema.time.as_deref()
                            {
                                self.bind_time_input.trim()
                            } else if Some(target_column.as_str())
                                == target.dataset.schema.temperature.as_deref()
                            {
                                self.bind_temperature_input.trim()
                            } else if Some(target_column.as_str())
                                == target.dataset.schema.mass.as_deref()
                            {
                                self.bind_mass_input.trim()
                            } else {
                                target_column.as_str()
                            };
                            if !incoming_names
                                .iter_names()
                                .any(|name| name.as_str() == source_column)
                            {
                                return Err(format!(
                                    "Incoming file is missing column '{source_column}' required for '{target_column}'"
                                ));
                            }
                            if selected_sources.iter().any(|name| name == source_column) {
                                return Err("Time, temperature, and mass bindings must use distinct incoming columns".into());
                            }
                            selected_sources.push(source_column.to_owned());
                            let target_dtype = target
                                .dataset
                                .frame
                                .clone()
                                .collect_schema()
                                .map_err(|error| error.to_string())?
                                .get(target_column)
                                .ok_or_else(|| {
                                    format!("Target schema is missing '{target_column}'")
                                })?
                                .clone();
                            let expression =
                                col(source_column).cast(target_dtype).alias(target_column);
                            expressions.push(expression);
                        }
                        if selected_sources.len() != incoming_names.len() {
                            return Err("The continuation file has extra columns; map or remove them explicitly before appending".into());
                        }
                        continuation.frame = continuation.frame.select(expressions);
                        continuation.schema = target.dataset.schema.clone();
                        let candidate = target
                            .dataset
                            .append_continuation(&continuation)
                            .map_err(|e| format!("{e:?}"))?;
                        self.append_preview_rows = Some(
                            candidate
                                .frame
                                .clone()
                                .collect()
                                .map_err(|e| e.to_string())?
                                .height(),
                        );
                        self.file_path_input = path.to_string_lossy().to_string();
                        self.append_preview = Some(candidate);
                        self.append_result_id = format!("{}_continued", target_id);
                        self.append_preview_signature = Some(self.append_signature(model));
                        Ok(())
                    })();
                    self.last_error = result.err();
                    if let Some(error) = &self.last_error {
                        let _ =
                            model.push_message(&format!("Continuation preview failed: {error}"));
                    }
                } else if self.import_file_format == Some(ImportFileFormat::TypedSeries) {
                    let path = if self.file_path_input.trim().is_empty() {
                        self.selected_file.clone()
                    } else {
                        Some(PathBuf::from(self.file_path_input.trim()))
                    };
                    match path {
                        Some(path) if path.is_file() => {
                            if self.replace_series_on_import {
                                self.request_open_typed_series(model, &path);
                            } else {
                                match self.import_typed_series_file(model, &path) {
                                    Ok(ids) => {
                                        info!(
                                            "Imported typed KiThe CSV with {} experiments",
                                            ids.len()
                                        );
                                        self.open = false;
                                        self.last_error = None;
                                    }
                                    Err(error) => {
                                        self.last_error = Some(error.to_string());
                                        let _ = model.push_message(&format!(
                                            "Series import failed: {error}"
                                        ));
                                    }
                                }
                            }
                        }
                        _ => {
                            self.last_error =
                                Some("Please select a valid typed CSV file".to_string());
                            let _ = model
                                .push_message("Series import failed: select a valid CSV file.");
                        }
                    }
                } else {
                    self.meta.thermal_mode = self.thermal_mode;
                    let parse_optional_number = |input: &str| -> Result<Option<f64>, ()> {
                        if input.trim().is_empty() {
                            Ok(None)
                        } else {
                            input.trim().parse::<f64>().map(Some).map_err(|_| ())
                        }
                    };
                    self.meta.heating_rate = if self.thermal_mode
                        == ExperimentThermalMode::Isothermal
                    {
                        None
                    } else {
                        match parse_optional_number(&self.heating_rate_input) {
                            Ok(value) => value,
                            Err(()) => {
                                self.last_error = Some("Heating rate must be a number".to_string());
                                let _ = model
                                    .push_message("Import failed: heating rate must be a number.");
                                return;
                            }
                        }
                    };
                    self.meta.isothermal_temperature = if self.thermal_mode
                        == ExperimentThermalMode::NonIsothermal
                        || self.use_measured_temperature_mean
                    {
                        None
                    } else {
                        match parse_optional_number(&self.isothermal_temp_input) {
                            Ok(value) => value,
                            Err(()) => {
                                self.last_error =
                                    Some("Isothermal temperature must be a number".to_string());
                                let _ = model.push_message(
                                    "Import failed: isothermal temperature must be a number.",
                                );
                                return;
                            }
                        }
                    };
                    if self.thermal_mode == ExperimentThermalMode::Unspecified
                        && self.meta.heating_rate.is_some()
                        && self.meta.isothermal_temperature.is_some()
                    {
                        self.last_error = Some("Choose Isothermal or Non-isothermal when both thermal metadata fields are provided".into());
                        let _ = model.push_message("Import failed: thermal regime is ambiguous.");
                        return;
                    }
                    if self
                        .meta
                        .heating_rate
                        .is_some_and(|rate| !rate.is_finite() || rate <= 0.0)
                    {
                        self.last_error =
                            Some("Heating rate must be finite and above 0 K/min".into());
                        let _ = model.push_message(
                            "Import failed: heating rate must be positive and finite.",
                        );
                        return;
                    }
                    if self
                        .meta
                        .isothermal_temperature
                        .is_some_and(|temperature| !temperature.is_finite() || temperature <= 0.0)
                    {
                        self.last_error =
                            Some("Isothermal temperature must be finite and above 0 K".into());
                        let _ = model.push_message(
                            "Import failed: isothermal temperature must be positive and finite.",
                        );
                        return;
                    }
                    if self.thermal_mode == ExperimentThermalMode::Isothermal
                        && !self.use_measured_temperature_mean
                    {
                        match self.meta.isothermal_temperature {
                            Some(temperature) if temperature.is_finite() && temperature > 0.0 => {}
                            _ => {
                                self.last_error = Some(
                                    "Enter a finite isothermal temperature above 0 K, or choose the measured-column mean".into(),
                                );
                                let _ = model.push_message(
                                    "Import failed: valid isothermal temperature is required.",
                                );
                                return;
                            }
                        }
                    }

                    let file_path_candidate = if self.file_path_input.trim().is_empty() {
                        self.selected_file.clone()
                    } else {
                        Some(PathBuf::from(self.file_path_input.trim()))
                    };

                    if let Some(file_path) = file_path_candidate {
                        if file_path.is_file() {
                            self.selected_file = Some(file_path.clone());
                            self.file_path_input = file_path.to_string_lossy().to_string();
                            if self.meta.id.trim().is_empty() {
                                self.meta.id = file_path
                                    .file_stem()
                                    .map(|stem| stem.to_string_lossy().into_owned())
                                    .filter(|stem| !stem.trim().is_empty())
                                    .unwrap_or_else(|| "experiment".to_string());
                            }
                            // Create a new experiment with the selected file and metadata
                            if let Err(e) =
                                model.push_from_file_with_meta(&file_path, self.meta.clone())
                            {
                                self.last_error = Some(format!("Error loading experiment: {e}"));
                                let _ = model.push_message(&format!("Import failed: {e}"));
                            } else {
                                let exp_id = self.meta.id.clone();
                                let apply_user_bindings =
                                    self.import_file_format != Some(ImportFileFormat::TypedDataset);
                                match self.finish_new_experiment_import_with_bindings(
                                    model,
                                    &exp_id,
                                    apply_user_bindings,
                                ) {
                                    Ok(()) => {
                                        info!("File parsed successfully");
                                        info!("Experiment ID set to: {}", exp_id);
                                        if let Ok(exp) = model.series.get_experiment_by_id(&exp_id)
                                        {
                                            info!(
                                                "found columns with names {:?}",
                                                exp.list_of_columns()
                                            );
                                        }
                                        self.last_error = None;
                                        self.open = false;
                                    }
                                    Err(e) => {
                                        if let TGAGUIError::PostCommitWarning(warning) = e {
                                            self.last_error = None;
                                            self.open = false;
                                            let _ =
                                                model.push_message(&format!("Warning: {warning}"));
                                        } else {
                                            self.last_error = Some(e.to_string());
                                            let _ =
                                                model.push_message(&format!("Import failed: {e}"));
                                            return;
                                        }
                                    }
                                }
                            }
                        } else {
                            self.last_error = Some("Selected path is not a file".to_string());
                            let _ =
                                model.push_message("Import failed: selected path is not a file.");
                        }
                    } else {
                        self.last_error = Some("Please select a file".to_string());
                        let _ = model.push_message("Import failed: select a file first.");
                    }
                }
            }

            if commit_append {
                let result = (|| -> Result<(String, Option<String>), String> {
                    let current_signature = self.append_signature(model);
                    if self.append_preview_signature.as_deref() != Some(current_signature.as_str())
                    {
                        self.append_preview = None;
                        self.append_preview_rows = None;
                        self.append_preview_signature = None;
                        return Err("Preview is stale. Rebuild it after changing the file, bindings, target, or data.".into());
                    }
                    let target_id = self
                        .append_target_id
                        .clone()
                        .ok_or_else(|| "Select a target experiment".to_string())?;
                    let candidate = self
                        .append_preview
                        .take()
                        .ok_or_else(|| "Preview the continuation before committing".to_string())?;
                    if self.append_replace_target {
                        model
                            .series
                            .try_mutate_by_id(&target_id, |experiment| {
                                experiment.dataset = candidate;
                                Ok(())
                            })
                            .map_err(|e| format!("{e:?}"))?;
                        let warning =
                            model
                                .rebuild_plots_for_experiment(&target_id)
                                .err()
                                .map(|error| match error {
                                    TGAGUIError::PostCommitWarning(message) => message,
                                    other => other.to_string(),
                                });
                        Ok((target_id, warning))
                    } else {
                        let id = self.append_result_id.trim();
                        if id.is_empty() {
                            return Err("Enter a new experiment ID".into());
                        }
                        let mut meta = model
                            .series
                            .get_experiment_by_id(&target_id)
                            .map_err(|e| format!("{e:?}"))?
                            .meta
                            .clone();
                        meta.id = id.to_string();
                        model
                            .series
                            .push(TGAExperiment {
                                dataset: candidate,
                                meta,
                            })
                            .map_err(|e| format!("{e:?}"))?;
                        let previous_view =
                            (model.interaction.view_range, model.interaction.view_y_range);
                        let warning = self
                            .create_default_plots_for_last_experiment(model)
                            .err()
                            .map(|error| {
                                model.plots.retain(|plot| plot.experiment_id != id);
                                model.interaction.view_range = previous_view.0;
                                model.interaction.view_y_range = previous_view.1;
                                format!(
                                    "Continuation experiment '{id}' was added, but its plots could not be created; data remain available ({error})"
                                )
                            });
                        Ok((id.to_string(), warning))
                    }
                })();
                match result {
                    Ok((_id, Some(warning))) => {
                        let _ = model.push_message(&format!("Warning: {warning}"));
                        self.open = false;
                        self.last_error = None;
                    }
                    Ok((id, None)) => {
                        let _ = model
                            .push_message(&format!("Continuation appended as experiment '{id}'"));
                        self.open = false;
                        self.last_error = None;
                    }
                    Err(error) => {
                        let _ = model.push_message(&format!("Continuation failed: {error}"));
                        self.last_error = Some(error);
                    }
                }
            }
        }
    }

    fn with_extension(filename: &str, ext: &str) -> String {
        let trimmed = filename.trim();
        let required_suffix = format!(".{}", ext);
        if trimmed.to_ascii_lowercase().ends_with(&required_suffix) {
            trimmed.to_string()
        } else {
            format!("{trimmed}{required_suffix}")
        }
    }

    fn open_save_document_dialog(&mut self, model: &PlotModel) {
        self.save_series_open = true;
        self.save_dialog_purpose = SaveDialogPurpose::Document;
        self.save_series_format = SaveSeriesFormat::Csv;
        self.save_scope = SaveScope::Series;
        self.save_series_last_status = None;
        self.last_error = None;
        if let Some(path) = &self.current_document_path {
            if let Some(parent) = path.parent() {
                self.current_dir = parent.to_path_buf();
            }
            self.save_series_filename_input = path
                .file_name()
                .and_then(|name| name.to_str())
                .unwrap_or("session.csv")
                .to_string();
        } else {
            let stem = model
                .list_of_experiments()
                .first()
                .filter(|id| !id.trim().is_empty())
                .cloned()
                .unwrap_or_else(|| "kinetics_session".to_string());
            self.save_series_filename_input = format!("{stem}.csv");
        }
    }

    fn open_export_dialog(&mut self, format: SaveSeriesFormat, model: &PlotModel) {
        self.save_series_open = true;
        self.save_dialog_purpose = SaveDialogPurpose::Export;
        self.save_series_format = format;
        self.save_scope = if matches!(format, SaveSeriesFormat::Csv | SaveSeriesFormat::Xlsx) {
            SaveScope::Series
        } else {
            SaveScope::Experiment
        };
        self.save_experiment_id = model.list_of_experiments().first().cloned();
        self.save_selected_experiments.clear();
        if let Some(id) = &self.save_experiment_id {
            self.save_selected_experiments.insert(id.clone());
        }
        self.save_selected_columns.clear();
        self.sync_export_columns(model);
        self.save_series_last_status = None;
        self.last_error = None;
        self.save_series_filename_input = match format {
            SaveSeriesFormat::Csv => "series_export.csv".to_string(),
            SaveSeriesFormat::Txt => self
                .save_experiment_id
                .as_deref()
                .filter(|id| !id.trim().is_empty())
                .map(|id| format!("{id}_export.txt"))
                .unwrap_or_else(|| "experiment_export.txt".to_string()),
            SaveSeriesFormat::Xlsx => "series_export.xlsx".to_string(),
        };
    }

    fn pick_open_file(&mut self, title: &str, typed_series_only: bool) -> Option<PathBuf> {
        #[cfg(test)]
        {
            let _ = (title, typed_series_only);
            return self.test_next_open_file.take();
        }
        #[cfg(not(test))]
        {
            let mut picker = rfd::FileDialog::new()
                .set_title(title)
                .set_directory(&self.current_dir);
            picker = if typed_series_only {
                picker.add_filter("KiThe series document", &["csv"])
            } else {
                picker.add_filter("CSV and text data", &["csv", "txt"])
            };
            picker.pick_file()
        }
    }

    fn pick_save_file(&mut self, title: &str, extension: &str) -> Option<PathBuf> {
        #[cfg(test)]
        {
            let _ = (title, extension);
            return self.test_next_save_file.take();
        }
        #[cfg(not(test))]
        {
            rfd::FileDialog::new()
                .set_title(title)
                .set_directory(&self.current_dir)
                .set_file_name(&self.save_series_filename_input)
                .add_filter(
                    &format!("{} file", extension.to_ascii_uppercase()),
                    &[extension],
                )
                .save_file()
        }
    }

    fn set_import_file(&mut self, path: PathBuf) {
        if let Some(parent) = path.parent() {
            self.current_dir = parent.to_path_buf();
        }
        self.file_path_input = path.to_string_lossy().to_string();
        self.update_import_file_format(&path);
        self.selected_file = Some(path.clone());
        if self.meta.id.trim().is_empty() {
            self.meta.id = path
                .file_stem()
                .map(|stem| stem.to_string_lossy().into_owned())
                .filter(|stem| !stem.trim().is_empty())
                .unwrap_or_else(|| "experiment".to_string());
        }
    }

    fn choose_import_file(&mut self) {
        if let Some(path) = self.pick_open_file("Import Data", false) {
            self.set_import_file(path);
            self.open_import_dialog();
        }
    }

    fn open_import_dialog(&mut self) {
        self.open = true;
        self.append_mode = false;
        self.import_file_format = self
            .selected_file
            .as_deref()
            .and_then(|path| detect_import_file_format(path).ok());
        self.replace_series_on_import = false;
        self.last_error = None;
    }

    fn open_document_dialog(&mut self, model: &mut PlotModel) {
        let Some(path) = self.pick_open_file("Open KiThe Document", true) else {
            return;
        };
        if detect_import_file_format(&path).ok() != Some(ImportFileFormat::TypedSeries) {
            let message = "Open Document accepts typed KiThe series CSV files. Use Import Data for raw files or typed datasets.";
            self.last_error = Some(message.to_string());
            let _ = model.push_message(message);
            return;
        }
        if let Some(parent) = path.parent() {
            self.current_dir = parent.to_path_buf();
        }
        self.request_open_typed_series(model, &path);
    }

    fn update_import_file_format(&mut self, path: &Path) {
        self.import_file_format = detect_import_file_format(path).ok();
    }

    fn append_signature(&self, model: &PlotModel) -> String {
        format!(
            "{}|{}|{}|{}|{}|{:?}|{:?}|{}|{}|{}",
            self.append_target_id.as_deref().unwrap_or_default(),
            self.selected_file
                .as_deref()
                .unwrap_or_else(|| Path::new(""))
                .display(),
            self.file_path_input.trim(),
            self.bind_time_input.trim(),
            self.bind_temperature_input.trim(),
            self.bind_time_unit,
            self.bind_temperature_unit,
            self.bind_mass_input.trim(),
            self.bind_mass_unit,
            model.series.data_revision(),
        )
    }

    fn choose_append_file(&mut self, model: &PlotModel) {
        let Some(path) = self.pick_open_file("Append to Experiment", false) else {
            return;
        };
        self.open_append_dialog(model);
        self.set_import_file(path);
    }

    fn open_append_dialog(&mut self, model: &PlotModel) {
        self.open = true;
        self.append_mode = true;
        self.last_error = None;
        self.append_preview = None;
        self.append_preview_rows = None;
        self.append_preview_signature = None;
        self.append_target_id = model
            .get_experiment_by_selected_curve()
            .ok()
            .or_else(|| model.list_of_experiments().into_iter().next());
        self.append_replace_target = false;
        self.append_result_id.clear();
        self.file_path_input.clear();
        self.selected_file = None;
        if let Some(id) = self.append_target_id.as_deref()
            && let Ok(experiment) = model.series.get_experiment_by_id(id)
        {
            self.bind_time_input = experiment.dataset.schema.time.clone().unwrap_or_default();
            self.bind_temperature_input = experiment
                .dataset
                .schema
                .temperature
                .clone()
                .unwrap_or_default();
            self.bind_mass_input = experiment.dataset.schema.mass.clone().unwrap_or_default();
            if let Some(unit) = experiment
                .dataset
                .schema
                .time
                .as_deref()
                .and_then(|name| experiment.dataset.schema.columns.get(name))
                .map(|meta| meta.unit)
            {
                self.bind_time_unit = unit;
            }
            if let Some(unit) = experiment
                .dataset
                .schema
                .temperature
                .as_deref()
                .and_then(|name| experiment.dataset.schema.columns.get(name))
                .map(|meta| meta.unit)
            {
                self.bind_temperature_unit = unit;
            }
            if let Some(unit) = experiment
                .dataset
                .schema
                .mass
                .as_deref()
                .and_then(|name| experiment.dataset.schema.columns.get(name))
                .map(|meta| meta.unit)
            {
                self.bind_mass_unit = unit;
            }
        }
    }

    fn export_experiment_ids(&self, model: &PlotModel) -> Result<Vec<String>, String> {
        match self.save_scope {
            SaveScope::Series
                if matches!(
                    self.save_series_format,
                    SaveSeriesFormat::Csv | SaveSeriesFormat::Xlsx
                ) =>
            {
                let ids = model.list_of_experiments();
                if ids.is_empty() {
                    Err("There are no experiments to export".to_string())
                } else {
                    Ok(ids)
                }
            }
            SaveScope::SelectedExperiments
                if matches!(
                    self.save_series_format,
                    SaveSeriesFormat::Csv | SaveSeriesFormat::Xlsx
                ) =>
            {
                let ids = model
                    .list_of_experiments()
                    .into_iter()
                    .filter(|id| self.save_selected_experiments.contains(id))
                    .collect::<Vec<_>>();
                if ids.is_empty() {
                    Err("Select at least one experiment to export".to_string())
                } else {
                    Ok(ids)
                }
            }
            SaveScope::SelectedExperiments => {
                Err("Selected-experiment scope requires CSV or Excel export".to_string())
            }
            SaveScope::Series => Err("Plain TXT export supports one experiment at a time".into()),
            SaveScope::Experiment => self
                .save_experiment_id
                .as_ref()
                .filter(|id| model.list_of_experiments().contains(id))
                .cloned()
                .map(|id| vec![id])
                .ok_or_else(|| "Select an experiment to export".to_string()),
        }
    }

    fn sync_export_columns(&mut self, model: &PlotModel) {
        let experiment_ids = model.list_of_experiments();
        self.save_selected_columns
            .retain(|id, _| experiment_ids.contains(id));
        for id in experiment_ids {
            let columns = model.list_of_columns(&id).unwrap_or_default();
            let selected = self
                .save_selected_columns
                .entry(id)
                .or_insert_with(|| columns.iter().cloned().collect());
            selected.retain(|name| columns.contains(name));
        }
    }

    fn selected_export_columns(
        &self,
        model: &PlotModel,
    ) -> Result<HashMap<String, Vec<String>>, String> {
        let mut result = HashMap::new();
        for id in self.export_experiment_ids(model)? {
            let all_columns = model
                .list_of_columns(&id)
                .map_err(|error| error.to_string())?;
            let selected = self
                .save_selected_columns
                .get(&id)
                .ok_or_else(|| format!("No columns are selected for experiment '{id}'"))?;
            let columns = all_columns
                .into_iter()
                .filter(|name| selected.contains(name))
                .collect::<Vec<_>>();
            if columns.is_empty() {
                return Err(format!("Select at least one column for experiment '{id}'"));
            }
            result.insert(id, columns);
        }
        Ok(result)
    }

    fn export_column_label(&self, model: &PlotModel, id: &str, column: &str) -> String {
        let metadata = model
            .series
            .get_experiment_by_id(id)
            .ok()
            .and_then(|experiment| experiment.dataset.schema.columns.get(column));
        let nature = metadata
            .map(|meta| column_nature_label(meta.nature))
            .unwrap_or("Unclassified");
        let unit = metadata
            .map(|meta| meta.unit.to_string())
            .unwrap_or_else(|| "Unknown".to_string());
        let binding = model
            .series
            .get_experiment_by_id(id)
            .ok()
            .and_then(|experiment| {
                let schema = &experiment.dataset.schema;
                [
                    ("time", schema.time.as_deref()),
                    ("temperature", schema.temperature.as_deref()),
                    ("mass", schema.mass.as_deref()),
                    ("relative mass", schema.alpha.as_deref()),
                    ("conversion", schema.eta.as_deref()),
                    ("mass rate", schema.dm_dt.as_deref()),
                    ("relative-mass rate", schema.dalpha_dt.as_deref()),
                    ("conversion rate", schema.deta_dt.as_deref()),
                    ("temperature rate", schema.dT_dt.as_deref()),
                    ("activation energy", schema.E.as_deref()),
                    ("R2", schema.R2.as_deref()),
                ]
                .into_iter()
                .find_map(|(role, bound_name)| (bound_name == Some(column)).then_some(role))
            });
        let binding = binding
            .map(|role| format!(" · active: {role}"))
            .unwrap_or_default();

        format!("{id}::{column} · {nature} · {unit}{binding}")
    }

    fn write_export(
        &self,
        model: &PlotModel,
        path: &Path,
        request: &ExportRequest,
    ) -> Result<HashMap<String, Vec<String>>, String> {
        let columns_by_experiment = &request.columns_by_experiment;
        match request.format {
            SaveSeriesFormat::Csv => {
                model
                    .to_csv_series_selected(path, &request.experiment_ids, columns_by_experiment)
                    .map_err(|error| error.to_string())?;
            }
            SaveSeriesFormat::Txt => {
                let id = request
                    .experiment_ids
                    .first()
                    .ok_or_else(|| "Select an experiment to export".to_string())?;
                model
                    .to_plain_txt_experiment(id, path, &columns_by_experiment[id])
                    .map_err(|error| error.to_string())?;
            }
            SaveSeriesFormat::Xlsx => {
                write_xlsx_export(model, path, &request.experiment_ids, columns_by_experiment)?;
            }
        }
        Ok(columns_by_experiment.clone())
    }

    fn prepare_export(&self, model: &PlotModel) -> Result<ExportRequest, String> {
        Ok(ExportRequest {
            format: self.save_series_format,
            experiment_ids: self.export_experiment_ids(model)?,
            columns_by_experiment: self.selected_export_columns(model)?,
        })
    }

    fn set_save_format(&mut self, format: SaveSeriesFormat) {
        if self.save_series_format == format {
            return;
        }
        self.save_series_format = format;
        if format == SaveSeriesFormat::Txt {
            self.save_scope = SaveScope::Experiment;
        }
        let stem = Path::new(self.save_series_filename_input.trim())
            .file_stem()
            .and_then(|name| name.to_str())
            .filter(|name| !name.is_empty())
            .unwrap_or("export");
        self.save_series_filename_input = format!("{stem}.{}", format.extension());
    }

    fn save_series_shortcut(&mut self, model: &mut PlotModel, save_as: bool) {
        if save_as || self.current_document_path.is_none() {
            self.open_save_document_dialog(model);
            return;
        }
        if let Some(path) = self.current_document_path.clone() {
            self.write_document(model, &path);
        }
    }

    fn is_document_dirty(&self, model: &PlotModel) -> bool {
        model.series.data_revision() != self.saved_data_revision
    }

    pub(crate) fn document_generation(&self) -> u64 {
        self.document_generation
    }

    pub(crate) fn reset_save_target_after_external_replacement(&mut self, model: &PlotModel) {
        self.current_document_path = None;
        self.saved_data_revision = model.series.data_revision();
        self.pending_file_write = None;
        self.pending_open_series = None;
        self.append_preview = None;
        self.append_preview_rows = None;
        self.append_preview_signature = None;
        self.document_generation = self.document_generation.wrapping_add(1);
    }

    fn write_document(&mut self, model: &mut PlotModel, path: &Path) -> bool {
        let result = write_atomically(path, |temporary_path| {
            model
                .to_csv_series(temporary_path)
                .map_err(|error| error.to_string())
        });
        match result {
            Ok(()) => {
                self.current_document_path = Some(path.to_path_buf());
                self.saved_data_revision = model.series.data_revision();
                self.last_error = None;
                let _ = model.push_message(&format!("Document saved to {}", path.display()));
                info!("Document saved to {}", path.display());
                if let Some(next_path) = self.pending_open_series.take() {
                    self.load_typed_series_document(model, &next_path);
                }
                true
            }
            Err(error) => {
                self.last_error = Some(format!("Failed to save document: {error}"));
                let _ = model.push_message(&format!("Document save failed: {error}"));
                false
            }
        }
    }

    fn commit_file_write(
        &mut self,
        model: &mut PlotModel,
        path: &Path,
        purpose: SaveDialogPurpose,
        export: Option<ExportRequest>,
    ) {
        if purpose == SaveDialogPurpose::Document {
            if self.write_document(model, path) {
                self.save_series_open = false;
            }
            return;
        }

        let result = export
            .ok_or_else(|| "Export request was not captured".to_string())
            .and_then(|request| {
                write_atomically(path, |temporary_path| {
                    self.write_export(model, temporary_path, &request)
                })
            });
        match result {
            Ok(_) => {
                self.last_error = None;
                self.save_series_last_status = Some(format!("Exported to {}", path.display()));
                let _ = model.push_message(&format!("Data exported to {}", path.display()));
                self.save_series_open = false;
            }
            Err(error) => {
                self.save_series_last_status = None;
                self.last_error = Some(format!("Export failed: {error}"));
                let _ = model.push_message(&format!("Data export failed: {error}"));
            }
        }
    }

    pub fn show_save_series_dialog(&mut self, ctx: &egui::Context, model: &mut PlotModel) {
        if !self.save_series_open {
            return;
        }

        self.sync_export_columns(model);
        let mut should_export = false;
        let mut dialog_open = self.save_series_open;

        let is_document = self.save_dialog_purpose == SaveDialogPurpose::Document;
        egui::Window::new(if is_document {
            "Save Document As"
        } else {
            self.save_series_format.dialog_title()
        })
            .open(&mut dialog_open)
            .resizable(true)
            .default_size([850.0, 570.0])
            .show(ctx, |ui| {
                ui.vertical(|ui| {
                        ui.heading(if is_document { "Document" } else { "Export Options" });
                        ui.separator();
                        if is_document {
                            ui.label("Save the complete series, all columns, metadata, bindings, provenance, and operation history.");
                            ui.label("Undo snapshots and plot display state are not part of the document.");
                        } else {
                        ui.horizontal(|ui| {
                            if ui
                                .selectable_label(
                                    self.save_series_format == SaveSeriesFormat::Csv,
                                    "Typed CSV",
                                )
                                .clicked()
                            {
                                self.set_save_format(SaveSeriesFormat::Csv);
                            }
                            if ui
                                .selectable_label(
                                    self.save_series_format == SaveSeriesFormat::Txt,
                                    "Plain TXT",
                                )
                                .clicked()
                            {
                                self.set_save_format(SaveSeriesFormat::Txt);
                            }
                            if ui
                                .selectable_label(
                                    self.save_series_format == SaveSeriesFormat::Xlsx,
                                    "Excel XLSX",
                                )
                                .clicked()
                            {
                                self.set_save_format(SaveSeriesFormat::Xlsx);
                            }
                        });
                        ui.label(match self.save_series_format {
                            SaveSeriesFormat::Csv => "KiThe metadata and bindings are preserved.",
                            SaveSeriesFormat::Txt => {
                                "Tab-delimited values only; no KiThe metadata."
                            }
                            SaveSeriesFormat::Xlsx => {
                                "Selected numeric columns are written to one worksheet per experiment; no KiThe metadata."
                            }
                        });

                        let experiment_ids = model.list_of_experiments();
                        if matches!(self.save_series_format, SaveSeriesFormat::Csv | SaveSeriesFormat::Xlsx) {
                            ui.horizontal(|ui| {
                                ui.selectable_value(
                                    &mut self.save_scope,
                                    SaveScope::Series,
                                    "Entire series",
                                );
                                ui.selectable_value(
                                    &mut self.save_scope,
                                    SaveScope::SelectedExperiments,
                                    "Selected experiments",
                                );
                                ui.selectable_value(
                                    &mut self.save_scope,
                                    SaveScope::Experiment,
                                    "One experiment",
                                );
                            });
                        } else {
                            self.save_scope = SaveScope::Experiment;
                        }

                        if self.save_scope == SaveScope::Experiment {
                            if let Some(first_id) = experiment_ids.first() {
                                if self
                                    .save_experiment_id
                                    .as_ref()
                                    .is_none_or(|id| !experiment_ids.contains(id))
                                {
                                    self.save_experiment_id = Some(first_id.clone());
                                }
                            }
                            egui::ComboBox::from_label("Experiment")
                                .selected_text(
                                    self.save_experiment_id
                                        .as_deref()
                                        .unwrap_or("Select experiment"),
                                )
                                .show_ui(ui, |ui| {
                                    for id in &experiment_ids {
                                        ui.selectable_value(
                                            &mut self.save_experiment_id,
                                            Some(id.clone()),
                                            id,
                                        );
                                    }
                                });
                        }

                        if self.save_scope == SaveScope::SelectedExperiments {
                            ui.horizontal(|ui| {
                                if ui.small_button("Select all experiments").clicked() {
                                    self.save_selected_experiments
                                        .extend(experiment_ids.iter().cloned());
                                }
                                if ui.small_button("Clear experiment selection").clicked() {
                                    self.save_selected_experiments.clear();
                                }
                            });
                            egui::ScrollArea::vertical()
                                .id_salt("export_experiment_selection")
                                .max_height(60.0)
                                .show(ui, |ui| {
                                    for id in &experiment_ids {
                                        let mut selected =
                                            self.save_selected_experiments.contains(id);
                                        let display_id = if id.trim().is_empty() {
                                            "(unnamed experiment)"
                                        } else {
                                            id.as_str()
                                        };
                                        if ui.checkbox(&mut selected, display_id).changed() {
                                            if selected {
                                                self.save_selected_experiments.insert(id.clone());
                                            } else {
                                                self.save_selected_experiments.remove(id);
                                            }
                                        }
                                    }
                                });
                        }

                        let visible_ids = match self.save_scope {
                            SaveScope::Series => experiment_ids.clone(),
                            SaveScope::SelectedExperiments => experiment_ids
                                .iter()
                                .filter(|id| self.save_selected_experiments.contains(*id))
                                .cloned()
                                .collect(),
                            SaveScope::Experiment => {
                                self.save_experiment_id.iter().cloned().collect()
                            }
                        };
                        ui.horizontal(|ui| {
                            if ui.small_button("Select all columns").clicked() {
                                for id in &visible_ids {
                                    if let (Some(selected), Ok(columns)) = (
                                        self.save_selected_columns.get_mut(id),
                                        model.list_of_columns(id),
                                    ) {
                                        selected.extend(columns);
                                    }
                                }
                            }
                            if ui.small_button("Clear column selection").clicked() {
                                for id in &visible_ids {
                                    if let Some(selected) = self.save_selected_columns.get_mut(id) {
                                        selected.clear();
                                    }
                                }
                            }
                        });
                        egui::ScrollArea::vertical()
                            .id_salt("export_column_selection")
                            .max_height(110.0)
                            .show(ui, |ui| {
                                for id in &visible_ids {
                                    let Ok(columns) = model.list_of_columns(id) else {
                                        continue;
                                    };
                                    ui.strong(id);
                                    ui.horizontal_wrapped(|ui| {
                                        for column in columns {
                                            let label =
                                                self.export_column_label(model, id, &column);
                                            if let Some(selected) =
                                                self.save_selected_columns.get_mut(id)
                                            {
                                                let mut checked = selected.contains(&column);
                                                if ui
                                                    .checkbox(
                                                        &mut checked,
                                                        label,
                                                    )
                                                    .changed()
                                                {
                                                    if checked {
                                                        selected.insert(column);
                                                    } else {
                                                        selected.remove(&column);
                                                    }
                                                }
                                            }
                                        }
                                    });
                                    ui.separator();
                                }
                            });
                        }

                        ui.label("The system save dialog will choose the destination folder.");
                        ui.add_space(6.0);
                        let label = ui.label("File name:");
                        let input = ui.text_edit_singleline(&mut self.save_series_filename_input);
                        let _ = input.labelled_by(label.id);
                        ui.label(format!("Extension: .{}", self.save_series_format.extension()));

                        if let Some(status) = &self.save_series_last_status {
                            ui.colored_label(egui::Color32::GREEN, status);
                        }
                        if let Some(err) = &self.last_error {
                            ui.colored_label(egui::Color32::RED, err);
                        }
                        ui.add_space(6.0);
                        if ui.button(if is_document { "Save Document" } else { "Export" }).clicked() {
                            should_export = true;
                        }
                });
            });

        let document_dialog_cancelled = self.save_dialog_purpose == SaveDialogPurpose::Document
            && self.save_series_open
            && !dialog_open
            && !should_export;
        self.save_series_open = dialog_open;
        if document_dialog_cancelled && self.pending_open_series.is_some() {
            self.pending_open_series = None;
            self.unsaved_changes_open = false;
        }
        if should_export {
            let export_request = if self.save_dialog_purpose == SaveDialogPurpose::Export {
                match self.prepare_export(model) {
                    Ok(request) => Some(request),
                    Err(error) => {
                        self.last_error = Some(error.clone());
                        let _ = model.push_message(&format!("Export failed: {error}"));
                        return;
                    }
                }
            } else {
                None
            };
            let filename = self.save_series_filename_input.trim();
            if filename.is_empty() {
                self.last_error = Some("File name cannot be empty".to_string());
                return;
            }
            let extension = if self.save_dialog_purpose == SaveDialogPurpose::Document {
                "csv"
            } else {
                self.save_series_format.extension()
            };
            let suggested_name = Self::with_extension(filename, extension);
            self.save_series_filename_input = suggested_name.clone();
            let title = if self.save_dialog_purpose == SaveDialogPurpose::Document {
                "Save KiThe Document"
            } else {
                self.save_series_format.dialog_title()
            };
            if let Some(mut save_path) = self.pick_save_file(title, extension) {
                if !save_path
                    .extension()
                    .and_then(|ext| ext.to_str())
                    .is_some_and(|ext| ext.eq_ignore_ascii_case(extension))
                {
                    save_path.set_extension(extension);
                }
                if let Some(parent) = save_path.parent() {
                    self.current_dir = parent.to_path_buf();
                }
                if save_path.exists() {
                    self.pending_file_write = Some(PendingFileWrite {
                        path: save_path,
                        purpose: self.save_dialog_purpose,
                        export: export_request,
                    });
                } else {
                    self.commit_file_write(
                        model,
                        &save_path,
                        self.save_dialog_purpose,
                        export_request,
                    );
                }
            }
        }

        let mut confirm_overwrite = false;
        let mut cancel_overwrite = false;
        if let Some(pending) = &self.pending_file_write {
            egui::Window::new("Confirm overwrite")
                .collapsible(false)
                .resizable(false)
                .show(ctx, |ui| {
                    ui.label(format!(
                        "'{}' already exists. Replace it?",
                        pending.path.display()
                    ));
                    ui.horizontal(|ui| {
                        confirm_overwrite = ui.button("Replace file").clicked();
                        cancel_overwrite = ui.button("Cancel").clicked();
                    });
                });
        }
        if cancel_overwrite {
            self.pending_file_write = None;
        } else if confirm_overwrite {
            if let Some(pending) = self.pending_file_write.take() {
                self.commit_file_write(model, &pending.path, pending.purpose, pending.export);
            }
        }
    }
}

impl Default for NewExperimentDialogState {
    fn default() -> Self {
        Self::new()
    }
}
//=========================================================================================================
//           MANAGE PLOT DIALOGUE
//========================================================================================================
impl NewExperimentDialogState {
    pub fn show_manage_plot_dialog(&mut self, ctx: &egui::Context, model: &mut PlotModel) {
        if !self.manage_plot_open {
            return;
        }

        let mut dialog_open = self.manage_plot_open;
        egui::Window::new("Manage Plots")
            .open(&mut dialog_open)
            .resizable(true)
            .default_size([450.0, 400.0])
            .show(ctx, |ui| {
                ui.heading("📊 Create Plots");
                ui.separator();

                let experiments = model.list_of_experiments();

                egui::ScrollArea::vertical().show(ui, |ui| {
                    let mut to_remove = Vec::new();

                    for (idx, config) in self.plot_configs.iter_mut().enumerate() {
                        ui.group(|ui| {
                            ui.horizontal(|ui| {
                                ui.label(format!("Plot {}", idx + 1));
                                if ui.button("❌").clicked() {
                                    to_remove.push(idx);
                                }
                            });

                            let label = ui.label("Experiment:");
                            let combo = egui::ComboBox::from_id_salt(format!("exp_{}", idx))
                                .selected_text(config.experiment.as_deref().unwrap_or("Select..."))
                                .show_ui(ui, |ui| {
                                    for exp_id in &experiments {
                                        ui.selectable_value(
                                            &mut config.experiment,
                                            Some(exp_id.clone()),
                                            exp_id,
                                        );
                                    }
                                });
                            let _ = combo.response.labelled_by(label.id);

                            if let Some(exp_id) = &config.experiment {
                                match model.list_of_columns(exp_id) {
                                    Ok(columns) => {
                                        let label = ui.label("X column:");
                                        let combo =
                                            egui::ComboBox::from_id_salt(format!("x_{}", idx))
                                                .selected_text(
                                                    config.x_column.as_deref().unwrap_or("None"),
                                                )
                                                .show_ui(ui, |ui| {
                                                    for col in &columns {
                                                        ui.selectable_value(
                                                            &mut config.x_column,
                                                            Some(col.clone()),
                                                            col,
                                                        );
                                                    }
                                                });
                                        let _ = combo.response.labelled_by(label.id);

                                        let label = ui.label("Y column:");
                                        let combo =
                                            egui::ComboBox::from_id_salt(format!("y_{}", idx))
                                                .selected_text(
                                                    config.y_column.as_deref().unwrap_or("None"),
                                                )
                                                .show_ui(ui, |ui| {
                                                    for col in &columns {
                                                        ui.selectable_value(
                                                            &mut config.y_column,
                                                            Some(col.clone()),
                                                            col,
                                                        );
                                                    }
                                                });
                                        let _ = combo.response.labelled_by(label.id);
                                    }
                                    Err(e) => {
                                        ui.colored_label(
                                            egui::Color32::RED,
                                            format!("Error: {:?} ", e),
                                        );
                                        ui.colored_label(
                                            egui::Color32::RED,
                                            format!("there is no experiment {:?} ", exp_id),
                                        );
                                    }
                                }
                            }
                        });
                        ui.add_space(5.0);
                    }

                    for idx in to_remove.iter().rev() {
                        self.plot_configs.remove(*idx);
                    }
                });

                ui.separator();

                if ui.button("➕ Add Plot").clicked() {
                    self.plot_configs.push(PlotConfig::default());
                }

                ui.add_space(10.0);

                if let Some(err) = &self.last_error {
                    ui.colored_label(egui::Color32::RED, err);
                }

                ui.add_space(10.0);

                if ui.button("✅ Create All").clicked() {
                    let mut success = true;
                    for config in &self.plot_configs {
                        if let (Some(exp_id), Some(x_col), Some(y_col)) =
                            (&config.experiment, &config.x_column, &config.y_column)
                        {
                            if let Err(e) = model.set_x(exp_id, x_col) {
                                self.last_error = Some(format!("Error setting X: {e}"));
                                success = false;
                                break;
                            } else if let Err(e) = model.set_y(exp_id, y_col) {
                                self.last_error = Some(format!("Error setting Y: {e}"));
                                success = false;
                                break;
                            } else if let Err(e) = model.create_points_for_curve(exp_id) {
                                self.last_error = Some(format!("Error creating curve: {e}"));
                                success = false;
                                break;
                            }
                        }
                    }

                    if success {
                        self.last_error = None;
                        self.manage_plot_open = false;
                        self.plot_configs = vec![PlotConfig::default()];
                        ctx.request_repaint();
                    }
                }
            });

        self.manage_plot_open = dialog_open;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Kinetics::experimental_kinetics::experiment_series_main::{
        CalibrationSnapshot, TGAExperiment,
    };
    use crate::Kinetics::experimental_kinetics::one_experiment_dataset::Unit;
    use crate::Kinetics::experimental_kinetics::testing_mod::VirtualTGA;
    use std::path::Path;

    fn typed_series() -> TGASeries {
        let mut series = TGASeries::new();
        let data = VirtualTGA {
            time: (0..40).map(|index| index as f64).collect(),
            temperature: vec![523.15; 40],
            mass: (0..40).map(|index| 100.0 - index as f64 * 0.5).collect(),
        };
        for id in ["typed_a", "typed_b"] {
            series
                .create_from_synthetic_data(
                    &data,
                    ExperimentMeta {
                        id: id.to_string(),
                        isothermal_temperature: Some(523.15),
                        thermal_mode: ExperimentThermalMode::Isothermal,
                        comment: Some(format!("restored {id}")),
                        ..ExperimentMeta::default()
                    },
                )
                .unwrap();
        }
        series
    }

    fn existing_experiment(id: &str) -> TGAExperiment {
        let data = VirtualTGA {
            time: (0..40).map(|index| index as f64).collect(),
            temperature: vec![500.0; 40],
            mass: (0..40).map(|index| 90.0 - index as f64 * 0.3).collect(),
        };
        let dataset = TGADataset::create_from_synthetic_data(&data).unwrap();
        TGAExperiment::new(dataset)
            .with_id(id)
            .with_isothermal_temperature(500.0)
    }

    fn write_typed_series() -> tempfile::NamedTempFile {
        let file = tempfile::NamedTempFile::new().unwrap();
        typed_series().to_csv_series(file.path()).unwrap();
        file
    }

    #[test]
    fn export_request_keeps_format_and_selection_snapshot() {
        let mut model = PlotModel::default();
        model.series = typed_series();
        let mut dialog = NewExperimentDialogState::new();
        dialog.save_series_format = SaveSeriesFormat::Csv;
        dialog.save_scope = SaveScope::SelectedExperiments;
        dialog.save_selected_experiments = ["typed_a".to_string()].into_iter().collect();
        dialog.save_selected_columns.insert(
            "typed_a".to_string(),
            ["time".to_string(), "mass".to_string()]
                .into_iter()
                .collect(),
        );

        let request = dialog.prepare_export(&model).unwrap();
        dialog.save_series_format = SaveSeriesFormat::Txt;
        dialog.save_selected_experiments.clear();
        dialog.save_selected_columns.clear();

        assert_eq!(request.format, SaveSeriesFormat::Csv);
        assert_eq!(request.experiment_ids, ["typed_a"]);
        let mut selected_columns = request.columns_by_experiment["typed_a"].clone();
        selected_columns.sort();
        assert_eq!(selected_columns, ["mass", "time"]);
    }

    #[test]
    fn append_preview_signature_changes_when_target_or_binding_changes() {
        let mut model = PlotModel::default();
        model.series = typed_series();
        let mut dialog = NewExperimentDialogState::new();
        dialog.append_target_id = Some("typed_a".into());
        dialog.selected_file = Some(PathBuf::from("continuation.csv"));
        dialog.file_path_input = "continuation.csv".into();
        dialog.bind_time_input = "t".into();
        dialog.bind_temperature_input = "T".into();
        dialog.bind_mass_input = "m".into();
        let initial = dialog.append_signature(&model);

        dialog.append_target_id = Some("typed_b".into());
        assert_ne!(initial, dialog.append_signature(&model));
        dialog.append_target_id = Some("typed_a".into());
        dialog.bind_mass_input = "signal".into();
        assert_ne!(initial, dialog.append_signature(&model));
    }

    #[test]
    fn import_format_detection_distinguishes_raw_and_both_typed_csv_formats() {
        let raw = tempfile::NamedTempFile::new().unwrap();
        std::fs::write(raw.path(), "t,T,m\n0,500,10\n1,501,9\n").unwrap();
        assert_eq!(
            detect_import_file_format(raw.path()).unwrap(),
            ImportFileFormat::Raw
        );

        let dataset = tempfile::NamedTempFile::new().unwrap();
        std::fs::write(
            dataset.path(),
            "# KiThe TGA Dataset\n# t [s], m [mg]\n# Bindings: {}\nt,m\n0,10\n",
        )
        .unwrap();
        assert_eq!(
            detect_import_file_format(dataset.path()).unwrap(),
            ImportFileFormat::TypedDataset
        );

        let series = write_typed_series();
        assert_eq!(
            detect_import_file_format(series.path()).unwrap(),
            ImportFileFormat::TypedSeries
        );

        let mut dialog = NewExperimentDialogState::new();
        dialog.replace_series_on_import = true;
        dialog.update_import_file_format(raw.path());
        assert!(dialog.replace_series_on_import);
        dialog.open_import_dialog();
        assert!(!dialog.replace_series_on_import);
    }

    #[test]
    fn open_document_loads_only_typed_series_without_opening_import_form() {
        let series = write_typed_series();
        let mut dialog = NewExperimentDialogState::new();
        let mut model = PlotModel::new();
        dialog.test_next_open_file = Some(series.path().to_path_buf());

        dialog.open_document_dialog(&mut model);

        assert!(!dialog.open);
        assert!(!dialog.append_mode);
        assert_eq!(model.series.len(), 2);
        assert_eq!(dialog.current_document_path.as_deref(), Some(series.path()));
        assert!(dialog.last_error.is_none());
    }

    #[test]
    fn open_document_rejects_raw_file_without_showing_binding_form() {
        let raw = tempfile::NamedTempFile::new().unwrap();
        std::fs::write(raw.path(), "t,T,m\n0,500,10\n1,501,9\n").unwrap();
        let mut dialog = NewExperimentDialogState::new();
        let mut model = PlotModel::new();
        dialog.test_next_open_file = Some(raw.path().to_path_buf());

        dialog.open_document_dialog(&mut model);

        assert!(!dialog.open);
        assert_eq!(model.series.len(), 0);
        assert!(
            dialog
                .last_error
                .as_deref()
                .unwrap()
                .contains("Import Data")
        );
    }

    #[test]
    fn failed_atomic_writer_leaves_previous_destination_untouched() {
        let destination = tempfile::NamedTempFile::new().unwrap();
        std::fs::write(destination.path(), "original contents").unwrap();

        let result = write_atomically(destination.path(), |_temporary_path| {
            Err::<(), _>("injected serialization failure".to_string())
        });

        assert!(result.is_err());
        assert_eq!(
            std::fs::read_to_string(destination.path()).unwrap(),
            "original contents"
        );
    }

    #[test]
    fn typed_series_import_adds_experiments_and_restores_metadata_atomically() {
        let file = write_typed_series();
        let mut model = PlotModel::new();
        model.series.push(existing_experiment("existing")).unwrap();
        let mut dialog = NewExperimentDialogState::new();

        let imported_ids = dialog
            .import_typed_series_file(&mut model, file.path())
            .unwrap();

        assert_eq!(imported_ids, vec!["typed_a", "typed_b"]);
        assert_eq!(
            model.list_of_experiments(),
            vec!["existing", "typed_a", "typed_b"]
        );
        for id in ["typed_a", "typed_b"] {
            let experiment = model.series.get_experiment_by_id(id).unwrap();
            assert_eq!(
                experiment.meta.comment.as_deref(),
                Some(format!("restored {id}").as_str())
            );
            assert_eq!(
                experiment.meta.thermal_mode,
                ExperimentThermalMode::Isothermal
            );
            assert_eq!(experiment.meta.isothermal_temperature, Some(523.15));
            assert_eq!(experiment.dataset.schema.time.as_deref(), Some("time"));
            assert_eq!(
                experiment.dataset.schema.temperature.as_deref(),
                Some("temperature")
            );
            assert_eq!(experiment.dataset.schema.mass.as_deref(), Some("mass"));
        }
        for id in ["typed_a", "typed_b"] {
            assert!(model.plots.iter().any(|curve| curve.experiment_id == id));
        }
    }

    #[test]
    fn typed_document_archives_calibration_without_overwriting_current_settings() {
        let path = tempfile::NamedTempFile::new().unwrap();
        let mut source_model = PlotModel::new();
        source_model.settings.set_calibration_coeffs(-2.5, 17.0);
        source_model.series = typed_series();
        source_model.to_csv_series(path.path()).unwrap();

        let mut model = PlotModel::new();
        model.settings.set_calibration_coeffs(3.25, 8.0);
        let mut dialog = NewExperimentDialogState::new();
        dialog.replace_series_on_import = true;
        dialog
            .import_typed_series_file(&mut model, path.path())
            .unwrap();

        let active_calibration = model.settings.calibration_line().unwrap();
        assert_eq!(active_calibration.k(), 3.25);
        assert_eq!(active_calibration.b(), 8.0);
        assert_eq!(
            model.series.calibration_snapshot(),
            Some(CalibrationSnapshot::new(-2.5, 17.0))
        );
        let contents = std::fs::read_to_string(path.path()).unwrap();
        assert!(contents.contains("#META\tCALIBRATION\t-2.5\t17"));
    }

    #[test]
    fn typed_series_id_conflict_does_not_partially_modify_current_series() {
        let file = write_typed_series();
        let mut model = PlotModel::new();
        model.series.push(existing_experiment("typed_a")).unwrap();
        let before_ids = model.list_of_experiments();
        let before_meta = model
            .series
            .get_experiment_by_id("typed_a")
            .unwrap()
            .meta
            .clone();
        let mut dialog = NewExperimentDialogState::new();

        assert!(
            dialog
                .import_typed_series_file(&mut model, file.path())
                .is_err()
        );

        assert_eq!(model.list_of_experiments(), before_ids);
        let after = &model.series.get_experiment_by_id("typed_a").unwrap().meta;
        assert_eq!(after.id, before_meta.id);
        assert_eq!(
            after.isothermal_temperature,
            before_meta.isothermal_temperature
        );
        assert!(model.series.get_experiment_by_id("typed_b").is_err());
    }

    #[test]
    fn typed_series_replace_discards_the_previous_series_only_when_requested() {
        let file = write_typed_series();
        let mut model = PlotModel::new();
        model.series.push(existing_experiment("old")).unwrap();
        let mut dialog = NewExperimentDialogState::new();
        dialog.replace_series_on_import = true;

        dialog
            .import_typed_series_file(&mut model, file.path())
            .unwrap();

        assert_eq!(model.list_of_experiments(), vec!["typed_a", "typed_b"]);
        assert!(model.series.get_experiment_by_id("old").is_err());
        assert!(model.plots.iter().all(|plot| plot.experiment_id != "old"));
    }

    #[test]
    fn typed_dataset_import_keeps_file_bindings_and_units_instead_of_applying_defaults() {
        let data = VirtualTGA {
            time: (0..40).map(|index| index as f64).collect(),
            temperature: vec![523.15; 40],
            mass: (0..40).map(|index| 100.0 - index as f64 * 0.5).collect(),
        };
        let dataset = TGADataset::create_from_synthetic_data(&data).unwrap();
        let file = tempfile::NamedTempFile::new().unwrap();
        dataset.to_csv_with_units(file.path()).unwrap();
        assert_eq!(
            detect_import_file_format(file.path()).unwrap(),
            ImportFileFormat::TypedDataset
        );

        let mut model = PlotModel::new();
        model
            .push_from_file_with_meta(
                file.path(),
                ExperimentMeta {
                    id: "typed_dataset".into(),
                    isothermal_temperature: Some(523.15),
                    thermal_mode: ExperimentThermalMode::Isothermal,
                    ..ExperimentMeta::default()
                },
            )
            .unwrap();
        let mut dialog = NewExperimentDialogState::new();

        dialog
            .finish_new_experiment_import_with_bindings(&mut model, "typed_dataset", false)
            .unwrap();

        let imported = &model
            .series
            .get_experiment_by_id("typed_dataset")
            .unwrap()
            .dataset;
        assert_eq!(imported.schema.time, dataset.schema.time);
        assert_eq!(imported.schema.temperature, dataset.schema.temperature);
        assert_eq!(imported.schema.mass, dataset.schema.mass);
        for name in ["time", "temperature", "mass"] {
            assert_eq!(
                imported.schema.columns[name].unit,
                dataset.schema.columns[name].unit
            );
            assert_eq!(
                imported.schema.columns[name].origin,
                dataset.schema.columns[name].origin
            );
        }
    }

    fn fixture_model() -> PlotModel {
        let mut model = PlotModel::new();
        model
            .push_from_file(&Path::new(env!("CARGO_MANIFEST_DIR")).join("TGAexample.txt"))
            .expect("failed to load fixture");
        model
    }

    #[test]
    fn new_experiment_autofills_common_bindings_from_column_names() {
        let model = fixture_model();
        let mut dialog = NewExperimentDialogState::new();

        dialog.autofill_bind_inputs_from_last_experiment(&model);

        assert_eq!(dialog.bind_time_input, "t");
        assert_eq!(dialog.bind_temperature_input, "T");
        assert_eq!(dialog.bind_mass_input, "m");
        assert_eq!(dialog.bind_time_unit, Unit::Second);
        assert_eq!(dialog.bind_temperature_unit, Unit::Celsius);
        assert_eq!(dialog.bind_mass_unit, Unit::MilliVolt);
    }

    #[test]
    fn new_experiment_select_creates_standard_plots_automatically() {
        let mut model = fixture_model();
        let mut dialog = NewExperimentDialogState::new();

        dialog.autofill_bind_inputs_from_last_experiment(&model);
        dialog
            .apply_bindings_to_last_experiment(&mut model)
            .expect("binding should succeed");
        dialog
            .create_default_plots_for_last_experiment(&mut model)
            .expect("plot creation should succeed");

        let exp_id = model.list_of_experiments()[0].clone();
        let curves = model
            .plots
            .iter()
            .filter(|curve| curve.experiment_id == exp_id)
            .collect::<Vec<_>>();

        assert_eq!(curves.len(), 2);
        assert!(
            curves
                .iter()
                .any(|curve| curve.x_name == "t" && curve.y_name == "T")
        );
        assert!(
            curves
                .iter()
                .any(|curve| curve.x_name == "t" && curve.y_name == "m")
        );
    }

    #[test]
    fn failed_new_experiment_binding_does_not_commit_partial_schema_or_history() {
        let mut model = fixture_model();
        let exp_id = model.list_of_experiments()[0].clone();
        let before = model.series.get_experiment_by_id(&exp_id).unwrap();
        let before_schema = before.dataset.schema.clone();
        let before_history = before
            .dataset
            .history_of_operations
            .vector_of_operations
            .len();
        let before_undo_depth = before.dataset.undo_depth();

        let mut dialog = NewExperimentDialogState::new();
        dialog.bind_mass_input = "m".to_owned();
        dialog.bind_temperature_input = "missing_temperature".to_owned();
        dialog.bind_time_input = "t".to_owned();

        assert!(
            dialog
                .apply_bindings_to_last_experiment(&mut model)
                .is_err()
        );

        let after = model.series.get_experiment_by_id(&exp_id).unwrap();
        assert_eq!(after.dataset.schema.time, before_schema.time);
        assert_eq!(after.dataset.schema.temperature, before_schema.temperature);
        assert_eq!(after.dataset.schema.mass, before_schema.mass);
        assert_eq!(
            after
                .dataset
                .history_of_operations
                .vector_of_operations
                .len(),
            before_history
        );
        assert_eq!(after.dataset.undo_depth(), before_undo_depth);
    }

    #[test]
    fn failed_new_experiment_setup_rolls_back_imported_experiment_and_curves() {
        let mut model = fixture_model();
        let source = model.series.experiments().last().unwrap().dataset.clone();
        model
            .series
            .push(TGAExperiment::new(source).with_id("rollback_candidate"))
            .unwrap();

        let mut dialog = NewExperimentDialogState::new();
        dialog.bind_mass_input = "m".to_owned();
        dialog.bind_temperature_input = "missing_temperature".to_owned();
        dialog.bind_time_input = "t".to_owned();
        let before_ids = model.list_of_experiments();

        assert!(
            dialog
                .finish_new_experiment_import(&mut model, "rollback_candidate")
                .is_err()
        );

        assert_eq!(
            model.list_of_experiments(),
            before_ids[..before_ids.len() - 1]
        );
        assert!(
            !model
                .plots
                .iter()
                .any(|curve| curve.experiment_id == "rollback_candidate")
        );
    }

    #[test]
    fn isothermal_t_m_import_binds_without_temperature_and_creates_mass_plot() {
        let file = tempfile::NamedTempFile::new().unwrap();
        std::fs::write(file.path(), "t,m\n0.0,10.0\n1.0,9.8\n2.0,9.5\n").unwrap();
        let mut model = PlotModel::new();
        let meta = ExperimentMeta {
            id: "isothermal_no_sensor".into(),
            isothermal_temperature: Some(523.15),
            thermal_mode: ExperimentThermalMode::Isothermal,
            ..ExperimentMeta::default()
        };
        model
            .push_from_file_with_meta(file.path(), meta)
            .expect("load t/m data");

        let mut dialog = NewExperimentDialogState::new();
        dialog.thermal_mode = ExperimentThermalMode::Isothermal;
        dialog.bind_time_input = "t".into();
        dialog.bind_mass_input = "m".into();
        dialog.bind_time_unit = Unit::Second;
        dialog.bind_mass_unit = Unit::Milligram;
        dialog
            .finish_new_experiment_import(&mut model, "isothermal_no_sensor")
            .expect("isothermal import should not require a measured T column");

        let experiment = model
            .series
            .get_experiment_by_id("isothermal_no_sensor")
            .unwrap();
        assert!(experiment.dataset.schema.temperature.is_none());
        assert_eq!(experiment.meta.isothermal_temperature, Some(523.15));
        assert_eq!(
            experiment.meta.resolved_thermal_mode(),
            ExperimentThermalMode::Isothermal
        );
        assert_eq!(
            model
                .plots
                .iter()
                .filter(|curve| curve.experiment_id == "isothermal_no_sensor")
                .count(),
            1
        );
    }

    #[test]
    fn explicit_sensor_mean_becomes_isothermal_temperature_in_kelvin() {
        let mut model = fixture_model();
        let exp_id = model.list_of_experiments()[0].clone();
        let mean_celsius = model
            .series
            .get_experiment_by_id(&exp_id)
            .unwrap()
            .dataset
            .mean_on_column("T")
            .unwrap();
        let valid_count = model
            .series
            .get_experiment_by_id(&exp_id)
            .unwrap()
            .dataset
            .preview_rows(&["T"], usize::MAX)
            .unwrap()
            .columns["T"]
            .iter()
            .flatten()
            .filter(|value| value.is_finite())
            .count();
        let mut dialog = NewExperimentDialogState::new();
        dialog.thermal_mode = ExperimentThermalMode::Isothermal;
        dialog.use_measured_temperature_mean = true;
        dialog.bind_time_input = "t".into();
        dialog.bind_temperature_input = "T".into();
        dialog.bind_mass_input = "m".into();
        dialog.bind_time_unit = Unit::Second;
        dialog.bind_temperature_unit = Unit::Celsius;
        dialog.bind_mass_unit = Unit::MilliVolt;

        dialog
            .finish_new_experiment_import(&mut model, &exp_id)
            .expect("sensor mean should be adopted after binding");
        let experiment = model.series.get_experiment_by_id(&exp_id).unwrap();
        let actual = experiment.meta.isothermal_temperature.unwrap();
        assert!((actual - (mean_celsius + 273.15)).abs() < 1e-9);
        assert_eq!(
            experiment.meta.resolved_thermal_mode(),
            ExperimentThermalMode::Isothermal
        );
        assert_eq!(experiment.dataset.schema.temperature.as_deref(), Some("T"));
        assert!(
            model
                .message
                .contains(&format!("{valid_count} valid measurements"))
        );
        assert!(model.message.contains("K"));
    }

    #[test]
    fn stale_sensor_mean_option_is_ignored_after_switching_to_nonisothermal() {
        let mut model = fixture_model();
        let exp_id = model.list_of_experiments()[0].clone();
        let mut dialog = NewExperimentDialogState::new();
        dialog.thermal_mode = ExperimentThermalMode::NonIsothermal;
        dialog.use_measured_temperature_mean = true;
        dialog.isothermal_temp_input = "523".into();
        dialog.select_thermal_mode(ExperimentThermalMode::Isothermal);
        dialog.use_measured_temperature_mean = true;
        dialog.isothermal_temp_input = "523".into();
        dialog.select_thermal_mode(ExperimentThermalMode::NonIsothermal);
        assert!(!dialog.use_measured_temperature_mean);
        assert!(dialog.isothermal_temp_input.is_empty());
        dialog.heating_rate_input = "5.0".into();
        dialog.bind_time_input = "t".into();
        dialog.bind_temperature_input = "T".into();
        dialog.bind_mass_input = "m".into();
        dialog.bind_time_unit = Unit::Second;
        dialog.bind_temperature_unit = Unit::Celsius;
        dialog.bind_mass_unit = Unit::MilliVolt;
        model
            .series
            .try_mutate_by_id(&exp_id, |experiment| {
                experiment.meta.thermal_mode = ExperimentThermalMode::NonIsothermal;
                experiment.meta.heating_rate = Some(5.0);
                experiment.meta.isothermal_temperature = None;
                Ok(())
            })
            .unwrap();

        dialog
            .finish_new_experiment_import(&mut model, &exp_id)
            .unwrap();

        let experiment = model.series.get_experiment_by_id(&exp_id).unwrap();
        assert_eq!(
            experiment.meta.resolved_thermal_mode(),
            ExperimentThermalMode::NonIsothermal
        );
        assert_eq!(experiment.meta.heating_rate, Some(5.0));
        assert!(experiment.meta.isothermal_temperature.is_none());
    }

    #[test]
    fn switching_from_nonisothermal_to_isothermal_import_uses_only_setpoint() {
        let mut model = fixture_model();
        let exp_id = model.list_of_experiments()[0].clone();
        let mut dialog = NewExperimentDialogState::new();
        dialog.thermal_mode = ExperimentThermalMode::NonIsothermal;
        dialog.heating_rate_input = "5.0".into();
        dialog.select_thermal_mode(ExperimentThermalMode::Isothermal);
        assert!(dialog.heating_rate_input.is_empty());
        dialog.isothermal_temp_input = "523.15".into();
        dialog.bind_time_input = "t".into();
        dialog.bind_temperature_input = "T".into();
        dialog.bind_mass_input = "m".into();
        dialog.bind_time_unit = Unit::Second;
        dialog.bind_temperature_unit = Unit::Celsius;
        dialog.bind_mass_unit = Unit::MilliVolt;
        model
            .series
            .try_mutate_by_id(&exp_id, |experiment| {
                experiment.meta.thermal_mode = ExperimentThermalMode::Isothermal;
                experiment.meta.isothermal_temperature = Some(523.15);
                experiment.meta.heating_rate = None;
                Ok(())
            })
            .unwrap();

        dialog
            .finish_new_experiment_import(&mut model, &exp_id)
            .unwrap();

        let experiment = model.series.get_experiment_by_id(&exp_id).unwrap();
        assert_eq!(
            experiment.meta.thermal_mode,
            ExperimentThermalMode::Isothermal
        );
        assert_eq!(experiment.meta.isothermal_temperature, Some(523.15));
        assert!(experiment.meta.heating_rate.is_none());
    }
}
