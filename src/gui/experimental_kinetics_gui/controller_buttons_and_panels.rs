//! Top menu bar: File Manager, Math, Kinetic Methods,
//! Direct Problem and Test Options with stub handlers.
#![allow(non_camel_case_types)]
#![allow(non_snake_case)]
use crate::Kinetics::experimental_kinetics::exp_engine_api::XY;
use crate::Kinetics::experimental_kinetics::experiment_series_main::{
    ExperimentMeta, TGAExperiment,
};
use crate::Kinetics::experimental_kinetics::one_experiment_dataset::{TGADataset, Unit};
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
use std::path::Path;
use std::path::PathBuf;

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
                new_experiment_dialog.open_import_dialog();
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
                if ui.button("Solve IVP").clicked() {
                    direct_problem_state.open();
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
            .button("Import Data")
            .on_hover_text(
                "Open a data file and import its columns as an experiment. Shortcut: Ctrl/Cmd+O.",
            )
            .clicked()
        {
            new_experiment_dialog.open_import_dialog();
        }
        if ui
            .button("New experiment")
            .on_hover_text("Open the dialog to create an experiment from imported or entered data.")
            .clicked()
        {
            new_experiment_dialog.open = true;
            new_experiment_dialog.last_error = None;
        }
        if ui.button("Append to experiment")
            .on_hover_text("Import a continuation, preview its time alignment, then append it or replace the selected experiment.")
            .clicked()
        {
            new_experiment_dialog.open_append_dialog(model);
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
        if ui
            .button("Save")
            .on_hover_text(
                "Save the current experiment series to its current file. Shortcut: Ctrl/Cmd+S.",
            )
            .clicked()
        {
            new_experiment_dialog.save_series_shortcut(model, false);
        }
        if ui
            .button("Save As CSV")
            .on_hover_text("Choose a path and export the series as CSV.")
            .clicked()
        {
            new_experiment_dialog.open_save_series_dialog(SaveSeriesFormat::Csv);
            ui.ctx().request_repaint();
        }
        if ui
            .button("Save As TXT")
            .on_hover_text("Choose a path and export the series as TXT.")
            .clicked()
        {
            new_experiment_dialog.open_save_series_dialog(SaveSeriesFormat::Txt);
            ui.ctx().request_repaint();
        }
        if ui
            .button("Plot Editor")
            .on_hover_text(
                "Open the print-ready plot editor for journal/report figures. This does not change the live interactive plots or experiment data.",
            )
            .clicked()
        {
            Self::save_as_image(model, kithe_plot_window);
        }
    }

    fn save_as_image(model: &mut PlotModel, kithe_plot_window: &mut KiThePlotWindowState) {
        match kithe_plot_window.open_from_model(model) {
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
                "Choose X or Y, then select a plot range. Removes rows outside that range in place; no new column."
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
                                    ui.selectable_value(&mut temperature_source, col.clone(), col);
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
                        model
                            .set_active_temperature_source_for_experiment(&id, &temperature_source);
                        ui.ctx().request_repaint();
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
        .on_hover_text("When enabled, arithmetic changes only the selected curve's X/Y values inside its plotted range; otherwise it changes the whole selected column.");
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
                info!("averaging column {}", column);
                let r = if state.apply_only_to_selected_chart {
                    model.mean_on_interval_on_own_range_for_selected()?
                } else {
                    model.mean_on_column_for_selected()?
                };
                info!("AVERAGE ON COLUMN {} =  {}", column, r);
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
}

impl SaveSeriesFormat {
    fn extension(self) -> &'static str {
        match self {
            Self::Csv => "csv",
            Self::Txt => "txt",
        }
    }

    fn dialog_title(self) -> &'static str {
        match self {
            Self::Csv => "Save Series As CSV",
            Self::Txt => "Save Series As TXT",
        }
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
    pub current_dir: PathBuf,
    pub selected_file: Option<PathBuf>,
    file_path_input: String,
    browse_disks_level: bool,
    pub last_error: Option<String>,
    save_series_format: SaveSeriesFormat,
    save_series_filename_input: String,
    save_series_last_status: Option<String>,
    last_saved_series_path: Option<PathBuf>,
    heating_rate_input: String,
    isothermal_temp_input: String,
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
            ("New Experiment", &mut self.open),
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
            && ["Save Series As CSV", "Save Series As TXT"]
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
            current_dir: std::env::current_dir().unwrap_or_else(|_| PathBuf::from(".")),
            selected_file: None,
            file_path_input: String::new(),
            browse_disks_level: false,
            last_error: None,
            save_series_format: SaveSeriesFormat::Csv,
            save_series_filename_input: "series_export.csv".to_string(),
            save_series_last_status: None,
            last_saved_series_path: None,
            heating_rate_input: String::new(),
            isothermal_temp_input: String::new(),
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
        let result = (|| {
            model.series.try_mutate_by_id(exp_id, |exp| {
                exp.check_nulls_for_operation_borrowed("parsing file")?;
                exp.sample_all_columns_even_layers_table(50, Vec::new())
            })?;
            self.autofill_bind_inputs_from_last_experiment(model);
            self.apply_bindings_to_last_experiment(model)?;
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

        let time_col = model.series.get_time_col(&exp_id).or_else(|_| {
            model
                .list_of_columns(&exp_id)
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
            model.series.get_temperature_col(&exp_id).ok(),
            model.series.get_mass_col(&exp_id).ok(),
        ]
        .into_iter()
        .flatten()
        {
            if y_col == time_col {
                continue;
            }
            model.set_x(&exp_id, &time_col)?;
            model.set_y(&exp_id, &y_col)?;
            model.create_points_for_curve(&exp_id)?;
            created_any = true;
        }

        if created_any {
            model.reset_view();
        }

        Ok(())
    }

    /// Show the new experiment dialog
    pub fn show_new_experiment_dialogue(&mut self, ctx: &egui::Context, model: &mut PlotModel) {
        let mut should_load_experiment = false;
        let mut commit_append = false;

        if self.open {
            let mut dialog_open = self.open;
            egui::Window::new(if self.append_mode { "Append to Experiment" } else { "New Experiment" })
                .open(&mut dialog_open)
                .resizable(true)
                .default_size(if self.append_mode { [760.0, 640.0] } else { [600.0, 400.0] })
                .show(ctx, |ui| {
                    ui.horizontal(|ui| {
                        ui.label("📁 Current directory:");
                        if ui.button("Up").clicked() {
                            self.move_up();
                        }
                        if self.browse_disks_level {
                            ui.label("<Disks>");
                        } else {
                            ui.label(self.current_dir.to_string_lossy().as_ref());
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

                    // File browser panel (left side)
                    egui::Panel::left("file_browser")
                        .resizable(true)
                        .default_size(300.0)
                        .show(ui, |ui| {
                            self.show_file_browser(ui);
                        });

                    // Metadata panel (right side)
                    egui::Panel::right("metadata_panel")
                        .resizable(true)
                        .default_size(250.0)
                        .show(ui, |ui| {
                            ui.heading("📝 Metadata");
                            if ui.button(if self.append_mode { "Preview continuation" } else { "Import experiment" }).clicked() {
                                should_load_experiment = true;
                            }
                            ui.separator();

                            if self.append_mode {
                                let experiments = model.list_of_experiments();
                                if self.append_target_id.as_ref().is_none_or(|id| !experiments.contains(id)) {
                                    self.append_target_id = experiments.first().cloned();
                                }
                                egui::ComboBox::from_label("Append to")
                                    .selected_text(self.append_target_id.as_deref().unwrap_or("Select experiment"))
                                    .show_ui(ui, |ui| {
                                        for id in &experiments {
                                            ui.selectable_value(&mut self.append_target_id, Some(id.clone()), id);
                                        }
                                    });
                                ui.checkbox(&mut self.append_replace_target, "Replace the original experiment");
                                if !self.append_replace_target {
                                    let label = ui.label("New experiment ID:");
                                    let input = ui.text_edit_singleline(&mut self.append_result_id);
                                    let _ = input.labelled_by(label.id);
                                }
                                ui.label("Incoming columns must match the selected experiment's column names and order. Time is preserved if already continuous; otherwise the segment is shifted to the next typical time step.");
                            }

                            let label = ui.label("File path:");
                            let input = ui.text_edit_singleline(&mut self.file_path_input);
                            let _ = input.labelled_by(label.id);

                            if !self.append_mode {
                                ui.add_space(5.0);
                                let label = ui.label("Experiment ID:");
                                let input = ui.text_edit_singleline(&mut self.meta.id);
                                let _ = input.labelled_by(label.id);

                                ui.add_space(5.0);
                                let label = ui.label("Heating Rate (K/min):");
                                let input = ui.text_edit_singleline(&mut self.heating_rate_input);
                                let _ = input.labelled_by(label.id);

                                ui.add_space(5.0);
                                let label = ui.label("Isothermal Temperature (K):");
                                let input = ui.text_edit_singleline(&mut self.isothermal_temp_input);
                                let _ = input.labelled_by(label.id);

                                ui.add_space(5.0);
                                let label = ui.label("Comment:");
                                let comment = self.meta.comment.get_or_insert_with(String::new);
                                let input = ui.text_edit_multiline(comment);
                                let _ = input.labelled_by(label.id);
                            }

                            ui.add_space(8.0);
                            self.show_bindings_row(ui);

                            ui.add_space(10.0);

                            if let Some(err) = &self.last_error {
                                ui.colored_label(egui::Color32::RED, err);
                                ui.add_space(5.0);
                            }
                        });
                });
            self.open = dialog_open;

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
                        self.append_preview = Some(candidate);
                        self.append_result_id = format!("{}_continued", target_id);
                        self.file_path_input = path.to_string_lossy().to_string();
                        Ok(())
                    })();
                    self.last_error = result.err();
                    if let Some(error) = &self.last_error {
                        let _ =
                            model.push_message(&format!("Continuation preview failed: {error}"));
                    }
                } else {
                    self.meta.heating_rate = if self.heating_rate_input.trim().is_empty() {
                        None
                    } else {
                        match self.heating_rate_input.trim().parse::<f64>() {
                            Ok(v) => Some(v),
                            Err(_) => {
                                self.last_error = Some("Heating rate must be a number".to_string());
                                let _ = model
                                    .push_message("Import failed: heating rate must be a number.");
                                return;
                            }
                        }
                    };
                    self.meta.isothermal_temperature =
                        if self.isothermal_temp_input.trim().is_empty() {
                            None
                        } else {
                            match self.isothermal_temp_input.trim().parse::<f64>() {
                                Ok(v) => Some(v),
                                Err(_) => {
                                    self.last_error =
                                        Some("Isothermal temperature must be a number".to_string());
                                    let _ = model.push_message(
                                        "Import failed: isothermal temperature must be a number.",
                                    );
                                    return;
                                }
                            }
                        };

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
                                match self.finish_new_experiment_import(model, &exp_id) {
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
                                            let _ = model.push_message(&format!("Warning: {warning}"));
                                        } else {
                                            self.last_error = Some(e.to_string());
                                            let _ = model.push_message(&format!("Import failed: {e}"));
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
                        let warning = model
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
                        let previous_view = (
                            model.interaction.view_range,
                            model.interaction.view_y_range,
                        );
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
                        let _ = model.push_message(&format!("Continuation appended as experiment '{id}'"));
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

    fn move_up(&mut self) {
        if self.browse_disks_level {
            return;
        }
        if let Some(parent) = self.current_dir.parent() {
            self.current_dir = parent.to_path_buf();
        } else {
            self.browse_disks_level = true;
        }
        self.selected_file = None;
    }

    fn list_windows_disks() -> Vec<PathBuf> {
        let mut disks = Vec::new();
        for letter in b'A'..=b'Z' {
            let drive = format!("{}:\\", letter as char);
            let drive_path = PathBuf::from(drive);
            if Path::new(&drive_path).exists() {
                disks.push(drive_path);
            }
        }
        disks
    }

    fn show_file_browser(&mut self, ui: &mut egui::Ui) {
        ui.heading("Files");
        ui.separator();

        if self.browse_disks_level {
            ui.label("Available disks:");
            egui::ScrollArea::vertical().show(ui, |ui| {
                for disk in Self::list_windows_disks() {
                    let disk_label = disk.to_string_lossy().to_string();
                    if ui.button(disk_label).clicked() {
                        self.current_dir = disk;
                        self.browse_disks_level = false;
                        self.selected_file = None;
                    }
                }
            });
            return;
        }

        if let Ok(entries) = std::fs::read_dir(&self.current_dir) {
            egui::ScrollArea::vertical().show(ui, |ui| {
                for entry in entries.filter_map(Result::ok) {
                    let path = entry.path();
                    let file_name_str = entry.file_name().to_string_lossy().to_string();

                    if path.is_dir() {
                        if ui.button(format!("[D] {}", file_name_str)).clicked() {
                            self.current_dir = path;
                            self.selected_file = None;
                        }
                    } else {
                        let is_selected = self.selected_file.as_ref() == Some(&path);
                        if ui
                            .selectable_label(is_selected, format!("[F] {}", file_name_str))
                            .clicked()
                        {
                            self.file_path_input = path.to_string_lossy().to_string();
                            self.selected_file = Some(path);
                        }
                    }
                }
            });
        } else {
            ui.label("Unable to read directory");
        }
    }

    fn show_directory_browser(&mut self, ui: &mut egui::Ui) {
        ui.heading("Directories");
        ui.separator();

        if self.browse_disks_level {
            ui.label("Available disks:");
            egui::ScrollArea::vertical().show(ui, |ui| {
                for disk in Self::list_windows_disks() {
                    let disk_label = disk.to_string_lossy().to_string();
                    if ui.button(disk_label).clicked() {
                        self.current_dir = disk;
                        self.browse_disks_level = false;
                    }
                }
            });
            return;
        }

        if let Ok(entries) = std::fs::read_dir(&self.current_dir) {
            egui::ScrollArea::vertical().show(ui, |ui| {
                for entry in entries.filter_map(Result::ok) {
                    let path = entry.path();
                    if !path.is_dir() {
                        continue;
                    }
                    let file_name_str = entry.file_name().to_string_lossy().to_string();
                    if ui.button(format!("[D] {}", file_name_str)).clicked() {
                        self.current_dir = path;
                    }
                }
            });
        } else {
            ui.label("Unable to read directory");
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

    fn open_save_series_dialog(&mut self, format: SaveSeriesFormat) {
        self.save_series_open = true;
        self.save_series_format = format;
        self.save_series_last_status = None;
        self.last_error = None;
        self.save_series_filename_input = match format {
            SaveSeriesFormat::Csv => "series_export.csv".to_string(),
            SaveSeriesFormat::Txt => "series_export.txt".to_string(),
        };
    }

    fn open_import_dialog(&mut self) {
        self.open = true;
        self.append_mode = false;
        self.last_error = None;
    }

    fn open_append_dialog(&mut self, model: &PlotModel) {
        self.open = true;
        self.append_mode = true;
        self.last_error = None;
        self.append_preview = None;
        self.append_preview_rows = None;
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

    fn save_series_shortcut(&mut self, model: &mut PlotModel, save_as: bool) {
        if save_as {
            self.open_save_series_dialog(SaveSeriesFormat::Csv);
            return;
        }

        let Some(path) = self.last_saved_series_path.clone() else {
            self.open_save_series_dialog(SaveSeriesFormat::Csv);
            return;
        };

        match model.to_csv_series(&path) {
            Ok(()) => {
                self.last_error = None;
                let _ = model.push_message(&format!("Series saved to {}", path.display()));
                info!("Series saved to {}", path.display());
            }
            Err(err) => {
                self.last_error = Some(format!("Failed to save series: {err}"));
                let _ = model.push_message(&format!("Saving series failed: {err}"));
                self.save_series_open = true;
            }
        }
    }

    pub fn show_save_series_dialog(&mut self, ctx: &egui::Context, model: &mut PlotModel) {
        if !self.save_series_open {
            return;
        }

        let mut should_save_series = false;
        let mut dialog_open = self.save_series_open;

        egui::Window::new(self.save_series_format.dialog_title())
            .open(&mut dialog_open)
            .resizable(true)
            .default_size([700.0, 420.0])
            .show(ctx, |ui| {
                ui.horizontal(|ui| {
                    ui.label("Current directory:");
                    if ui.button("Up").clicked() {
                        self.move_up();
                    }
                    if self.browse_disks_level {
                        ui.label("<Disks>");
                    } else {
                        ui.label(self.current_dir.to_string_lossy().as_ref());
                    }
                });

                ui.separator();

                egui::Panel::left("save_series_dir_browser")
                    .resizable(true)
                    .default_size(360.0)
                    .show(ui, |ui| {
                        self.show_directory_browser(ui);
                    });

                egui::Panel::right("save_series_meta_panel")
                    .resizable(true)
                    .default_size(280.0)
                    .show(ui, |ui| {
                        ui.heading("Save Options");
                        ui.separator();
                        ui.label("Directory:");
                        ui.label(self.current_dir.to_string_lossy().as_ref());
                        ui.add_space(8.0);

                        let label = ui.label("File name:");
                        let input = ui.text_edit_singleline(&mut self.save_series_filename_input);
                        let _ = input.labelled_by(label.id);
                        ui.label(format!(
                            "Extension will be: .{}",
                            self.save_series_format.extension()
                        ));
                        ui.add_space(8.0);

                        if let Some(status) = &self.save_series_last_status {
                            ui.colored_label(egui::Color32::GREEN, status);
                        }
                        if let Some(err) = &self.last_error {
                            ui.colored_label(egui::Color32::RED, err);
                        }
                        ui.add_space(6.0);

                        if ui.button("Save").clicked() {
                            should_save_series = true;
                        }
                    });
            });

        self.save_series_open = dialog_open;

        if should_save_series {
            let filename = self.save_series_filename_input.trim();
            if filename.is_empty() {
                self.last_error = Some("File name cannot be empty".to_string());
                return;
            }

            let output_file_name =
                Self::with_extension(filename, self.save_series_format.extension());
            let mut save_path = self.current_dir.clone();
            save_path.push(output_file_name);

            match model.to_csv_series(&save_path) {
                Ok(()) => {
                    self.last_error = None;
                    self.last_saved_series_path = Some(save_path.clone());
                    self.save_series_last_status =
                        Some(format!("Saved to {}", save_path.to_string_lossy()));
                    self.save_series_open = false;
                }
                Err(err) => {
                    self.save_series_last_status = None;
                    self.last_error = Some(format!("Failed to save series: {err}"));
                    let _ = model.push_message(&format!("Saving series failed: {err}"));
                }
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
    use crate::Kinetics::experimental_kinetics::experiment_series_main::TGAExperiment;
    use crate::Kinetics::experimental_kinetics::one_experiment_dataset::Unit;
    use std::path::Path;

    fn fixture_model() -> PlotModel {
        let mut model = PlotModel::new();
        model
            .push_from_file(Path::new("src/assets/TGAexample.txt"))
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
}
