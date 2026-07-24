//! Top menu bar: File Manager, Math, Kinetic Methods,
//! Direct Problem and Test Options with stub handlers.
#![allow(non_camel_case_types)]
#![allow(non_snake_case)]
use crate::Kinetics::experimental_kinetics::exp_engine_api::XY;
use crate::Kinetics::experimental_kinetics::experiment_series_main::ExperimentMeta;
use crate::Kinetics::experimental_kinetics::one_experiment_dataset::Unit;
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
use log::info;
use std::path::Path;
use std::path::PathBuf;

//===================================================================================
//  DROP DOWN MENUE AT THE TOP OF WINDOW
//====================================================================================
/// High-level entrypoint for top menus and panels.
pub struct TopDropDownMenues;

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
        ui.horizontal_wrapped(|ui| {
            ui.menu_button("File Manager", |ui| {
                FileManager::show(ui, model, new_experiment_dialog, kithe_plot_window);
            });

            ui.menu_button("Math", |ui| {
                mathematics.show_menu(ui);
            });

            ui.menu_button("Kinetic Methods", |ui| {
                KineticMethods::show_menu(ui, kinetic_methods);
            });

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
        if ui.button("Import Data").clicked() {
            Self::import_data(model);
        }
        if ui.button("New experiment").clicked() {
            new_experiment_dialog.open = true;
            new_experiment_dialog.last_error = None;
        }

        if ui.button("Append data to experiment").clicked() {}

        if ui.button("Export Data").clicked() {
            Self::export_data(model);
        }

        if ui.button("Manage Data").clicked() {
            Self::manage_data(model);
        }

        if ui.button("Manage Plot").clicked() {
            new_experiment_dialog.manage_plot_open = true;
            ui.ctx().request_repaint();
        }
        ui.separator();
        if ui.button("Save As CSV").clicked() {
            new_experiment_dialog.open_save_series_dialog(SaveSeriesFormat::Csv);
            ui.ctx().request_repaint();
        }
        if ui.button("Save As TXT").clicked() {
            new_experiment_dialog.open_save_series_dialog(SaveSeriesFormat::Txt);
            ui.ctx().request_repaint();
        }
        if ui.button("Save As Excel").clicked() {}
        if ui.button("Save As Image").clicked() {
            Self::save_as_image(model, kithe_plot_window);
        }
    }

    fn import_data(_model: &mut PlotModel) {
        info!("Stub: Import Data clicked");
    }

    fn export_data(_model: &mut PlotModel) {
        info!("Stub: Export Data clicked");
    }

    fn manage_data(_model: &mut PlotModel) {
        info!("Stub: Manage Data clicked");
    }

    fn save_as_image(model: &mut PlotModel, kithe_plot_window: &mut KiThePlotWindowState) {
        match kithe_plot_window.open_from_model(model) {
            Ok(()) => {
                let _ = model.push_message("KiThe Plot Redactor window opened");
            }
            Err(err) => {
                let message = match err {
                    TGAGUIError::TGADomainError(domain) => format!("{:?}", domain),
                    TGAGUIError::SettingsErrors(msg) => msg,
                    TGAGUIError::BindingError(msg) => msg,
                };
                let _ = model.push_message(&format!("Save As Image failed: {message}"));
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
            "golden pipeline",
            "Clear Selected",
            "Reset View",
            "Refresh",
            "cut before time",
            "cut selected",
            "Zoom To Selection",
            "move time to zero",
            "Delete Plot",
            "Reverse",
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
                    if ui
                        .add_sized(button_size, egui::Button::new(*label))
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
            ui.label("Input value:");
            ui.add(egui::DragValue::new(&mut state.input_value).speed(0.1));

            ui.separator();

            // X/Y selector: sets state.x_or_y
            if ui.button("X").clicked() {
                state.x_or_y = Some(XY::X);
            }
            if ui.button("Y").clicked() {
                state.x_or_y = Some(XY::Y);
            }

            ui.separator();

            // Optional output column name for transforms that materialize a derived column.
            let mut tmp = state.input_string.clone().unwrap_or_default();
            ui.label("output column:");
            if ui.text_edit_singleline(&mut tmp).changed() {
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
                    ui.label("Mass source:");
                    let mut mass_source = model
                        .active_mass_source_for_experiment(&id)
                        .unwrap_or_else(|_| String::new());
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

                    ui.label("Temperature source:");
                    let mut temperature_source = model
                        .active_temperature_source_for_experiment(&id)
                        .unwrap_or_else(|_| String::new());
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
                    if ui
                        .add_sized(button_size_small, egui::Button::new(*label))
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
        );
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
            "Reverse" => {
                model.reverse_selected()?;
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
                info!(
                    "Stub: {action} clicked (input_value={}, selected_only={})",
                    state.input_value, state.apply_only_to_selected_chart
                );
            }
        }
        Ok(())
    }

    fn set_action_error(model: &mut PlotModel, err: TGAGUIError) {
        let text = match err {
            TGAGUIError::TGADomainError(domain) => format!("{:?}", domain),
            TGAGUIError::SettingsErrors(msg) => msg,
            TGAGUIError::BindingError(msg) => msg,
        };
        let _ = model.push_message(&text);
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
}

impl NewExperimentDialogState {
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
            ui.label(label);
            ui.text_edit_singleline(input);
            egui::ComboBox::from_id_salt(combo_id)
                .selected_text(unit.to_string())
                .show_ui(ui, |ui| {
                    for option in options {
                        ui.selectable_value(unit, *option, option.to_string());
                    }
                });
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
        let Some(last_exp) = model.series.experiments.last() else {
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
        let Some(last_exp) = model.series.experiments.last() else {
            return Ok(());
        };
        let exp_id = last_exp.meta.id.clone();

        let bind_mass = self.bind_mass_input.trim();
        if !bind_mass.is_empty() {
            model
                .series
                .bind_mass(&exp_id, bind_mass, self.bind_mass_unit)
                .map_err(|e| TGAGUIError::BindingError(format!("bind m failed: {:?}", e)))?;
            match model.series.get_mass_col(&exp_id) {
                Ok(mass_name) => {
                    info!("binded mass column {}", mass_name);
                }
                Err(e) => {
                    return Err(TGAGUIError::BindingError(format!(
                        "failed to read bound mass column for '{}': {:?}",
                        exp_id, e
                    )));
                }
            }
        }

        let bind_temperature = self.bind_temperature_input.trim();
        if !bind_temperature.is_empty() {
            model
                .series
                .bind_temperature(&exp_id, bind_temperature, self.bind_temperature_unit)
                .map_err(|e| TGAGUIError::BindingError(format!("bind T failed: {:?}", e)))?;
            match model.series.get_temperature_col(&exp_id) {
                Ok(t_name) => {
                    info!("binded temperature column {}", t_name);
                }
                Err(e) => {
                    return Err(TGAGUIError::BindingError(format!(
                        "failed to read bound temperature column for '{}': {:?}",
                        exp_id, e
                    )));
                }
            }
        }

        let bind_time = self.bind_time_input.trim();
        if !bind_time.is_empty() {
            model
                .series
                .bind_time(&exp_id, bind_time, self.bind_time_unit)
                .map_err(|e| TGAGUIError::BindingError(format!("bind t failed: {:?}", e)))?;
            match model.series.get_time_col(&exp_id) {
                Ok(t_name) => {
                    info!("binded time column {}", t_name);
                }
                Err(e) => {
                    return Err(TGAGUIError::BindingError(format!(
                        "failed to read bound time column for '{}': {:?}",
                        exp_id, e
                    )));
                }
            }
        }

        Ok(())
    }

    fn create_default_plots_for_last_experiment(
        &self,
        model: &mut PlotModel,
    ) -> Result<(), TGAGUIError> {
        let Some(last_exp) = model.series.experiments.last() else {
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

        if self.open {
            let mut dialog_open = self.open;
            egui::Window::new("New Experiment")
                .open(&mut dialog_open)
                .resizable(true)
                .default_size([600.0, 400.0])
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
                            ui.separator();

                            ui.label("File path:");
                            ui.text_edit_singleline(&mut self.file_path_input);

                            ui.add_space(5.0);

                            ui.label("Experiment ID:");
                            ui.text_edit_singleline(&mut self.meta.id);

                            ui.add_space(5.0);

                            ui.label("Heating Rate (K/min):");
                            ui.text_edit_singleline(&mut self.heating_rate_input);

                            ui.add_space(5.0);

                            ui.label("Isothermal Temperature (K):");
                            ui.text_edit_singleline(&mut self.isothermal_temp_input);

                            ui.add_space(5.0);

                            ui.label("Comment:");
                            let comment = self.meta.comment.get_or_insert_with(String::new);
                            ui.text_edit_multiline(comment);

                            ui.add_space(8.0);
                            self.show_bindings_row(ui);

                            ui.add_space(10.0);

                            if let Some(err) = &self.last_error {
                                ui.colored_label(egui::Color32::RED, err);
                                ui.add_space(5.0);
                            }

                            // Select button
                            if ui.button("✅ Select").clicked() {
                                should_load_experiment = true;
                            }
                        });
                });
            self.open = dialog_open;

            // Handle experiment loading outside the window closure
            if should_load_experiment {
                self.meta.heating_rate = if self.heating_rate_input.trim().is_empty() {
                    None
                } else {
                    match self.heating_rate_input.trim().parse::<f64>() {
                        Ok(v) => Some(v),
                        Err(_) => {
                            self.last_error = Some("Heating rate must be a number".to_string());
                            return;
                        }
                    }
                };
                self.meta.isothermal_temperature = if self.isothermal_temp_input.trim().is_empty() {
                    None
                } else {
                    match self.isothermal_temp_input.trim().parse::<f64>() {
                        Ok(v) => Some(v),
                        Err(_) => {
                            self.last_error =
                                Some("Isothermal temperature must be a number".to_string());
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
                        // Create a new experiment with the selected file and metadata
                        if let Err(e) = model.push_from_file(&file_path) {
                            self.last_error = Some(format!("Error loading experiment: {:?}", e));
                        } else {
                            if let Some(exp) = model.series.experiments.last_mut() {
                                exp.meta.id = self.meta.id.clone();
                                exp.meta.heating_rate = self.meta.heating_rate;
                                exp.meta.isothermal_temperature = self.meta.isothermal_temperature;
                                exp.meta.comment = self.meta.comment.clone();
                                info!("File parsed successfully");
                                info!("Experiment ID set to: {}", exp.meta.id);
                                info!("found columns with names {:?}", &exp.list_of_columns());
                                exp.check_nulls_for_operation_borrowed("parsing file");
                                if let Err(e) =
                                    exp.sample_all_columns_even_layers_table(50, Vec::new())
                                {
                                    self.last_error = Some(format!("can't make a table: {:?}", e));
                                };
                            }
                            model.series.exp_map.clear();
                            for (idx, exp) in model.series.experiments.iter().enumerate() {
                                model.series.exp_map.insert(exp.meta.id.clone(), idx);
                            }
                            self.autofill_bind_inputs_from_last_experiment(model);
                            if let Err(e) = self.apply_bindings_to_last_experiment(model) {
                                self.last_error = Some(format!("{:?}", e));
                                return;
                            }
                            if let Err(e) = self.create_default_plots_for_last_experiment(model) {
                                self.last_error = Some(format!("{:?}", e));
                                return;
                            }
                            self.last_error = None;
                            self.open = false;
                        }
                    } else {
                        self.last_error = Some("Selected path is not a file".to_string());
                    }
                } else {
                    self.last_error = Some("Please select a file".to_string());
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

                        ui.label("File name:");
                        ui.text_edit_singleline(&mut self.save_series_filename_input);
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
                    self.save_series_last_status =
                        Some(format!("Saved to {}", save_path.to_string_lossy()));
                    self.save_series_open = false;
                }
                Err(err) => {
                    self.save_series_last_status = None;
                    self.last_error = Some(format!("Failed to save series: {:?}", err));
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

                            ui.label("Experiment:");
                            egui::ComboBox::from_id_salt(format!("exp_{}", idx))
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

                            if let Some(exp_id) = &config.experiment {
                                match model.list_of_columns(exp_id) {
                                    Ok(columns) => {
                                        ui.label("X column:");
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

                                        ui.label("Y column:");
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
                                self.last_error = Some(format!("Error setting X: {:?}", e));
                                success = false;
                                break;
                            } else if let Err(e) = model.set_y(exp_id, y_col) {
                                self.last_error = Some(format!("Error setting Y: {:?}", e));
                                success = false;
                                break;
                            } else if let Err(e) = model.create_points_for_curve(exp_id) {
                                self.last_error = Some(format!("Error creating curve: {:?}", e));
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
}
