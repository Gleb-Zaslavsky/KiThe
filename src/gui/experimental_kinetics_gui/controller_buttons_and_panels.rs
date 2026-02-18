//! Top menu bar: File Manager, Math, Kinetic Methods,
//! Direct Problem and Test Options with stub handlers.
#![allow(non_camel_case_types)]
#![allow(non_snake_case)]
use crate::Kinetics::experimental_kinetics::experiment_series_main::ExperimentMeta;
use crate::gui::experimental_kinetics_gui::model::PlotModel;
use crate::gui::experimental_kinetics_gui::model::TGAGUIError;
use eframe::egui;
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
        new_experiment_dialog: &mut NewExperimentDialogState,
    ) {
        ui.horizontal_wrapped(|ui| {
            ui.menu_button("File Manager", |ui| {
                FileManager::show(ui, model, new_experiment_dialog);
            });

            ui.menu_button("Math", |ui| {
                Mathematics::show(ui, model);
            });

            ui.menu_button("Kinetic Methods", |ui| {
                KineticMethods::show(ui, model);
            });

            ui.menu_button("Direct Problem", |ui| {
                DirectProblem::show(ui, model);
            });

            ui.menu_button("Test Options", |ui| {
                TestOptions::show(ui, model);
            });
        });
        ui.separator();
    }
}

// Group: File Manager
struct FileManager;

impl FileManager {
    fn show(
        ui: &mut egui::Ui,
        model: &mut PlotModel,
        new_experiment_dialog: &mut NewExperimentDialogState,
    ) {
        if ui.button("Import Data").clicked() {
            Self::import_data(model);
        }
        if ui.button("New experiment").clicked() {
            new_experiment_dialog.open = true;
            new_experiment_dialog.last_error = None;
        }

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

        if ui.button("Save As Image").clicked() {
            Self::save_as_image(model);
        }
    }

    fn import_data(_model: &mut PlotModel) {
        println!("Stub: Import Data clicked");
    }

    fn export_data(_model: &mut PlotModel) {
        println!("Stub: Export Data clicked");
    }

    fn manage_data(_model: &mut PlotModel) {
        println!("Stub: Manage Data clicked");
    }

    fn save_as_image(_model: &mut PlotModel) {
        println!("Stub: Save As Image clicked");
    }
}

// Group: Mathematics
struct Mathematics;

impl Mathematics {
    fn show(ui: &mut egui::Ui, model: &mut PlotModel) {
        if ui.button("Filter Data").clicked() {
            Self::filter_data(model);
        }

        if ui.button("Differentiate").clicked() {
            Self::differentiate(model);
        }

        if ui.button("Smooth").clicked() {
            Self::smooth(model);
        }

        if ui.button("Average Data").clicked() {
            Self::average_data(model);
        }

        if ui.button("Splines").clicked() {
            Self::splines(model);
        }
    }

    fn filter_data(_model: &mut PlotModel) {
        println!("Stub: Filter Data clicked");
    }

    fn differentiate(_model: &mut PlotModel) {
        println!("Stub: Differentiate clicked");
    }

    fn smooth(_model: &mut PlotModel) {
        println!("Stub: Smooth clicked");
    }

    fn average_data(_model: &mut PlotModel) {
        println!("Stub: Average Data clicked");
    }

    fn splines(_model: &mut PlotModel) {
        println!("Stub: Splines clicked");
    }
}

// Group: Kinetic Methods (stubs)
struct KineticMethods;

impl KineticMethods {
    fn show(ui: &mut egui::Ui, _model: &mut PlotModel) {
        if ui.button("Estimate Rates").clicked() {
            println!("Stub: Estimate Rates");
        }

        if ui.button("Fit Mechanism").clicked() {
            println!("Stub: Fit Mechanism");
        }
    }
}

// Group: Direct Problem (stubs)
struct DirectProblem;

impl DirectProblem {
    fn show(ui: &mut egui::Ui, _model: &mut PlotModel) {
        if ui.button("Solve IVP").clicked() {
            println!("Stub: Solve IVP");
        }

        if ui.button("Solve BVP").clicked() {
            println!("Stub: Solve BVP");
        }
    }
}

// Group: Test Options (stubs)
struct TestOptions;

impl TestOptions {
    fn show(ui: &mut egui::Ui, _model: &mut PlotModel) {
        if ui.button("Run Tests").clicked() {
            println!("Stub: Run Tests");
        }

        if ui.button("Show Logs").clicked() {
            println!("Stub: Show Logs");
        }
    }
}
//=================================================================================================

#[derive(Default)]
pub struct QuickActionPanelState {
    pub input_value: f64,
    pub apply_only_to_selected_chart: bool,
}

pub struct WrightPanelControllers;

impl WrightPanelControllers {
    pub fn quick_action_panel(
        ui: &mut egui::Ui,
        model: &mut PlotModel,
        state: &mut QuickActionPanelState,
        manage_plot_dialog: &mut NewExperimentDialogState,
    ) -> Result<(), TGAGUIError> {
        let labels = [
            "Manage Plots",
            "Clear Selection",
            "Reset View",
            "Refresh",
            "Zoom To Selection",
            "from C° to K",
            "move time to zero",
            "Delete",
            "Divide",
            "Sub",
            "Mul",
            "Add",
            "calc relative mass",
            "calc conversion",
            "from mV to mg",
            "from s to h",
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
                        let _ = Self::handle_quick_action(
                            label,
                            model,
                            state,
                            manage_plot_dialog,
                            ui.ctx(),
                        );
                    }
                    if (index + 1) % 4 == 0 {
                        ui.end_row();
                    }
                }
            });

        ui.separator();
        ui.horizontal(|ui| {
            ui.label("Input value:");
            ui.add(egui::DragValue::new(&mut state.input_value).speed(0.1));
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
        state: &QuickActionPanelState,
        manage_plot_dialog: &mut NewExperimentDialogState,
        ctx: &egui::Context,
    ) -> Result<(), TGAGUIError> {
        match action {
            "Manage Plots" => {
                manage_plot_dialog.manage_plot_open = true;
                ctx.request_repaint();
            }
            "Clear Selection" => {
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

            "from C° to K" => {
                model.from_C_to_K_of_selected()?;
                ctx.request_repaint();
            }
            "from s to h" => {
                model.from_s_to_h_of_selected()?;
                ctx.request_repaint();
            }

            _ => {
                println!(
                    "Stub: {action} clicked (input_value={}, selected_only={})",
                    state.input_value, state.apply_only_to_selected_chart
                );
            }
        }
        Ok(())
    }

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
        manage_plot_dialog: &mut NewExperimentDialogState,
    ) {
        if ui.button("Reset View").clicked() {
            model.reset_view();
        }

        if ui.button("Clear Selection").clicked() {
            model.clear_selection_rect();
            model.clear_selection();
        }

        if ui.button("📊 Manage Plots").clicked() {
            manage_plot_dialog.manage_plot_open = true;
            ui.ctx().request_repaint();
        }

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
    pub current_dir: PathBuf,
    pub selected_file: Option<PathBuf>,
    file_path_input: String,
    browse_disks_level: bool,
    pub last_error: Option<String>,
    heating_rate_input: String,
    isothermal_temp_input: String,
    meta: ExperimentMeta,
    plot_configs: Vec<PlotConfig>,
}

impl NewExperimentDialogState {
    /// Create a new dialog state with the current directory
    pub fn new() -> Self {
        Self {
            open: false,
            manage_plot_open: false,
            current_dir: std::env::current_dir().unwrap_or_else(|_| PathBuf::from(".")),
            selected_file: None,
            file_path_input: String::new(),
            browse_disks_level: false,
            last_error: None,
            heating_rate_input: String::new(),
            isothermal_temp_input: String::new(),
            meta: ExperimentMeta::new(),
            plot_configs: vec![PlotConfig::default()],
        }
    }

    /// Show the new experiment dialog
    pub fn show_dialogue(&mut self, ctx: &egui::Context, model: &mut PlotModel) {
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
                    egui::SidePanel::left("file_browser")
                        .resizable(true)
                        .default_width(300.0)
                        .show_inside(ui, |ui| {
                            self.show_file_browser(ui);
                        });

                    // Metadata panel (right side)
                    egui::SidePanel::right("metadata_panel")
                        .resizable(true)
                        .default_width(250.0)
                        .show_inside(ui, |ui| {
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
                                println!("File parsed successfully");
                                println!("Experiment ID set to: {}", exp.meta.id);
                                println!("found columns with names {:?}", &exp.list_of_columns());
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
}

impl Default for NewExperimentDialogState {
    fn default() -> Self {
        Self::new()
    }
}

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
