// Group: Test Options
use crate::Kinetics::experimental_kinetics::testing_mod::{
    AdvancedTGAConfig, ExperimentMode, KineticModel, NoiseConfig, NoiseKind, SpikeModel,
};
use crate::gui::experimental_kinetics_gui::gui_test::load_model_with_one_curve;
use crate::gui::experimental_kinetics_gui::model::PlotModel;
use crate::gui::experimental_kinetics_gui::settings::{CalibrationLine, Settings};

#[derive(Clone, Debug, Default)]
struct HelpDoc {
    sections: Vec<HelpSection>,
}
#[derive(Clone, Debug, Default)]
struct HelpSection {
    title: String,
    paragraphs: Vec<HelpParagraph>,
}
#[derive(Clone, Debug, Default)]
struct HelpParagraph {
    title: String,
    body: String,
}
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
enum HelpLanguage {
    Eng,
    Rus,
}
impl Default for HelpLanguage {
    fn default() -> Self {
        HelpLanguage::Eng
    }
}

pub struct TestOptions {
    /// Flag to show/hide the settings window
    show_settings_window: bool,
    /// Temporary settings being edited in the UI
    temp_k: String,
    temp_b: String,
    temp_n_points: String,
    temp_log_level: String,
    show_help_window: bool,
    help_lang: HelpLanguage,
    parsed_help: HelpDoc,
    help_error: Option<String>,
    // Synthetic data window
    show_synthetic_data: bool,

    // Synthetic data: general
    syn_n_points: String,
    syn_dt: String,
    syn_seed: String,

    // Experiment mode selection: 0 = Isothermal, 1 = NonIsothermal
    syn_mode: usize,
    // Isothermal
    syn_iso_temperatures: String, // csv list of K
    // NonIsothermal
    syn_ni_t0: String,
    syn_ni_heating_rates: String, // csv list of K/min

    // Kinetic model selection: 0 = ArrheniusSingle, 1 = ArrheniusTwoComponent
    syn_kmodel: usize,
    // Single
    syn_m0: String,
    syn_k0: String,
    syn_e: String,
    syn_r: String,
    // Two-component
    syn_m01: String,
    syn_k01: String,
    syn_e1: String,
    syn_m02: String,
    syn_k02: String,
    syn_e2: String,
    syn_r2: String,

    // Noise (mass)
    syn_mass_noise_enabled: bool,
    syn_mass_noise_kind: usize, // 0=Gaussian,1=Uniform,2=Drift
    syn_mass_sigma: String,
    syn_mass_ampl: String,
    syn_mass_slope: String,

    // Noise (temperature)
    syn_temp_noise_enabled: bool,
    syn_temp_noise_kind: usize,
    syn_temp_sigma: String,
    syn_temp_ampl: String,
    syn_temp_slope: String,

    // Spikes
    syn_spikes_enabled: bool,
    syn_spike_prob: String,
    syn_spike_ampl: String,

    syn_status: Option<String>,
}

impl Default for TestOptions {
    fn default() -> Self {
        Self {
            show_settings_window: false,
            temp_k: String::new(),
            temp_b: String::new(),
            temp_n_points: String::new(),
            temp_log_level: String::new(),
            show_help_window: false,
            help_lang: HelpLanguage::Eng,
            parsed_help: HelpDoc::default(),
            help_error: None,
            show_synthetic_data: false,

            syn_n_points: "70000".to_string(),
            syn_dt: "0.1".to_string(),
            syn_seed: "42".to_string(),
            syn_mode: 1,
            syn_iso_temperatures: "400.0".to_string(),
            syn_ni_t0: "400.0".to_string(),
            syn_ni_heating_rates: "1.0, 2.0, 3.0".to_string(),
            syn_kmodel: 0,
            syn_m0: "100.0".to_string(),
            syn_k0: "1e6".to_string(),
            syn_e: "80000.0".to_string(),
            syn_r: "8.314".to_string(),
            syn_m01: "50.0".to_string(),
            syn_k01: "1e6".to_string(),
            syn_e1: "90000.0".to_string(),
            syn_m02: "50.0".to_string(),
            syn_k02: "1e3".to_string(),
            syn_e2: "60000.0".to_string(),
            syn_r2: "8.314".to_string(),
            syn_mass_noise_enabled: false,
            syn_mass_noise_kind: 0,
            syn_mass_sigma: "5.0".to_string(),
            syn_mass_ampl: "1.0".to_string(),
            syn_mass_slope: "0.0".to_string(),
            syn_temp_noise_enabled: false,
            syn_temp_noise_kind: 0,
            syn_temp_sigma: "0.1".to_string(),
            syn_temp_ampl: "0.0".to_string(),
            syn_temp_slope: "0.0".to_string(),
            syn_spikes_enabled: false,
            syn_spike_prob: "0.01".to_string(),
            syn_spike_ampl: "1.0".to_string(),
            syn_status: None,
        }
    }
}

impl TestOptions {
    pub fn new() -> Self {
        Self::default()
    }

    /// Main UI display function
    pub fn show(ui: &mut egui::Ui, model: &mut PlotModel, test_options: &mut TestOptions) {
        if ui.button("Run Tests").clicked() {
            println!("Stub: Run Tests");
        }
        if ui.button("Load testing data").clicked() {
            println!("Stub: Load testing data");
            model.push_testing_data();
        }
        if ui.button("Create synthetic data").clicked() {
            test_options.show_synthetic_data = true;
        }
        if ui.button("Show and redact settings").clicked() {
            test_options.show_settings_window = true;
            // Initialize temp values from current settings
            test_options.load_current_settings(&model.settings);
        }

        if ui.button("Show Logs").clicked() {
            println!("Stub: Show Logs");
        }

        ui.separator();
        if ui.button("Help").clicked() {
            println!("Opening help window");
            test_options.show_help_window = true;
            let _ = test_options.reload_help();
            ui.ctx().request_repaint();
        }
    }

    /// Displays the settings editor window
    pub fn show_settings_window(&mut self, ctx: &egui::Context, model: &mut PlotModel) {
        let mut is_open = self.show_settings_window;
        egui::Window::new("TGA Settings")
            .open(&mut is_open)
            .resizable(true)
            .default_width(400.0)
            .show(ctx, |ui| {
                ui.heading("Current Settings");
                ui.separator();

                // Display current settings
                ui.group(|ui| {
                    ui.label("Active Configuration:");
                    if let Some(cal) = model.settings.calibration_line() {
                        ui.label(format!("Calibration k: {}", cal.k()));
                        ui.label(format!("Calibration b: {}", cal.b()));
                    }
                    if let Some(n) = model.settings.n_points() {
                        ui.label(format!("Number of points: {}", n));
                    }
                    if let Some(log) = model.settings.log_level() {
                        ui.label(format!("Log level: {}", log));
                    }
                });

                ui.add_space(10.0);
                ui.separator();
                ui.heading("Edit Settings");

                // Edit calibration line
                ui.group(|ui| {
                    ui.label("Calibration Line (mass = k * voltage + b):");
                    ui.horizontal(|ui| {
                        ui.label("k:");
                        ui.text_edit_singleline(&mut self.temp_k);
                    });
                    ui.horizontal(|ui| {
                        ui.label("b:");
                        ui.text_edit_singleline(&mut self.temp_b);
                    });
                });

                ui.add_space(5.0);

                // Edit number of points
                ui.horizontal(|ui| {
                    ui.label("Number of plot points:");
                    ui.text_edit_singleline(&mut self.temp_n_points);
                });

                ui.add_space(5.0);

                // Edit log level
                ui.horizontal(|ui| {
                    ui.label("Log level:");
                    egui::ComboBox::from_id_salt("log_level_combo")
                        .selected_text(&self.temp_log_level)
                        .show_ui(ui, |ui| {
                            ui.selectable_value(&mut self.temp_log_level, "off".to_string(), "off");
                            ui.selectable_value(
                                &mut self.temp_log_level,
                                "error".to_string(),
                                "error",
                            );
                            ui.selectable_value(
                                &mut self.temp_log_level,
                                "warn".to_string(),
                                "warn",
                            );
                            ui.selectable_value(
                                &mut self.temp_log_level,
                                "info".to_string(),
                                "info",
                            );
                            ui.selectable_value(
                                &mut self.temp_log_level,
                                "debug".to_string(),
                                "debug",
                            );
                            ui.selectable_value(
                                &mut self.temp_log_level,
                                "trace".to_string(),
                                "trace",
                            );
                        });
                });

                ui.add_space(10.0);
                ui.separator();

                // Action buttons
                ui.horizontal(|ui| {
                    if ui.button("Save Settings").clicked() {
                        if let Err(e) = self.apply_settings(model) {
                            eprintln!("Error applying settings: {:?}", e);
                            model.push_message(&format!("Error: {:?}", e)).ok();
                        } else {
                            model.push_message("Settings saved successfully").ok();
                        }
                    }

                    if ui.button("Reset to Defaults").clicked() {
                        model.settings = Settings::default();
                        self.load_current_settings(&model.settings);
                        model.push_message("Settings reset to defaults").ok();
                    }

                    if ui.button("Close").clicked() {
                        self.show_settings_window = false;
                    }
                });

                // Show any messages
                if !model.message.is_empty() {
                    ui.add_space(5.0);
                    ui.colored_label(egui::Color32::GREEN, &model.message);
                }
            });
        self.show_settings_window = is_open;
    }

    /// Synthetic data window
    pub fn show_synthetic_data_window(&mut self, ctx: &egui::Context, model: &mut PlotModel) {
        if !self.show_synthetic_data {
            return;
        }
        let mut open = self.show_synthetic_data;
        egui::Window::new("Create Synthetic TGA Data")
            .open(&mut open)
            .resizable(true)
            .vscroll(true)
            .default_width(720.0)
            .default_height(620.0)
            .show(ctx, |ui| {
                ui.heading("AdvancedTGAConfig");
                ui.separator();
                // General
                ui.group(|ui| {
                    ui.heading("General");
                    ui.horizontal(|ui| {
                        ui.label("n_points:");
                        ui.text_edit_singleline(&mut self.syn_n_points);
                        ui.label("dt (s):");
                        ui.text_edit_singleline(&mut self.syn_dt);
                        ui.label("seed:");
                        ui.text_edit_singleline(&mut self.syn_seed);
                    });
                });
                ui.add_space(6.0);
                // Experiment Mode
                ui.group(|ui| {
                    ui.heading("Experiment Mode");
                    egui::ComboBox::from_id_salt("syn_mode")
                        .selected_text(match self.syn_mode {
                            0 => "Isothermal",
                            _ => "NonIsothermal",
                        })
                        .show_ui(ui, |ui| {
                            ui.selectable_value(&mut self.syn_mode, 0, "Isothermal");
                            ui.selectable_value(&mut self.syn_mode, 1, "NonIsothermal");
                        });
                    match self.syn_mode {
                        0 => {
                            ui.horizontal(|ui| {
                                ui.label("temperatures (K, csv):");
                                ui.text_edit_singleline(&mut self.syn_iso_temperatures);
                            });
                        }
                        _ => {
                            ui.horizontal(|ui| {
                                ui.label("t0 (K):");
                                ui.text_edit_singleline(&mut self.syn_ni_t0);
                                ui.label("heating_rates (K/min, csv):");
                                ui.text_edit_singleline(&mut self.syn_ni_heating_rates);
                            });
                        }
                    }
                });
                ui.add_space(6.0);
                // Kinetic Model
                ui.group(|ui| {
                    ui.heading("Kinetic Model");
                    egui::ComboBox::from_id_salt("syn_kmodel")
                        .selected_text(match self.syn_kmodel {
                            0 => "ArrheniusSingle",
                            _ => "ArrheniusTwoComponent",
                        })
                        .show_ui(ui, |ui| {
                            ui.selectable_value(&mut self.syn_kmodel, 0, "ArrheniusSingle");
                            ui.selectable_value(&mut self.syn_kmodel, 1, "ArrheniusTwoComponent");
                        });
                    match self.syn_kmodel {
                        0 => {
                            egui::Grid::new("km_single").num_columns(8).show(ui, |ui| {
                                ui.label("m0");
                                ui.text_edit_singleline(&mut self.syn_m0);
                                ui.label("k0");
                                ui.text_edit_singleline(&mut self.syn_k0);
                                ui.label("E");
                                ui.text_edit_singleline(&mut self.syn_e);
                                ui.label("R");
                                ui.text_edit_singleline(&mut self.syn_r);
                            });
                        }
                        _ => {
                            egui::Grid::new("km_two").num_columns(14).show(ui, |ui| {
                                ui.label("m01");
                                ui.text_edit_singleline(&mut self.syn_m01);
                                ui.label("k01");
                                ui.text_edit_singleline(&mut self.syn_k01);
                                ui.label("E1");
                                ui.text_edit_singleline(&mut self.syn_e1);
                                ui.label("m02");
                                ui.text_edit_singleline(&mut self.syn_m02);
                                ui.label("k02");
                                ui.text_edit_singleline(&mut self.syn_k02);
                                ui.label("E2");
                                ui.text_edit_singleline(&mut self.syn_e2);
                                ui.label("R");
                                ui.text_edit_singleline(&mut self.syn_r2);
                            });
                        }
                    }
                });
                ui.add_space(6.0);
                // Noise and spikes
                ui.group(|ui| {
                    ui.heading("Mass noise");
                    ui.checkbox(&mut self.syn_mass_noise_enabled, "enabled");
                    egui::ComboBox::from_id_salt("mass_noise_kind")
                        .selected_text(match self.syn_mass_noise_kind {
                            0 => "Gaussian",
                            1 => "Uniform",
                            _ => "Drift",
                        })
                        .show_ui(ui, |ui| {
                            ui.selectable_value(&mut self.syn_mass_noise_kind, 0, "Gaussian");
                            ui.selectable_value(&mut self.syn_mass_noise_kind, 1, "Uniform");
                            ui.selectable_value(&mut self.syn_mass_noise_kind, 2, "Drift");
                        });
                    match self.syn_mass_noise_kind {
                        0 => {
                            ui.horizontal(|ui| {
                                ui.label("sigma");
                                ui.text_edit_singleline(&mut self.syn_mass_sigma);
                            });
                        }
                        1 => {
                            ui.horizontal(|ui| {
                                ui.label("amplitude");
                                ui.text_edit_singleline(&mut self.syn_mass_ampl);
                            });
                        }
                        _ => {
                            ui.horizontal(|ui| {
                                ui.label("slope");
                                ui.text_edit_singleline(&mut self.syn_mass_slope);
                            });
                        }
                    }
                    ui.separator();
                    ui.heading("Temperature noise");
                    ui.checkbox(&mut self.syn_temp_noise_enabled, "enabled");
                    egui::ComboBox::from_id_salt("temp_noise_kind")
                        .selected_text(match self.syn_temp_noise_kind {
                            0 => "Gaussian",
                            1 => "Uniform",
                            _ => "Drift",
                        })
                        .show_ui(ui, |ui| {
                            ui.selectable_value(&mut self.syn_temp_noise_kind, 0, "Gaussian");
                            ui.selectable_value(&mut self.syn_temp_noise_kind, 1, "Uniform");
                            ui.selectable_value(&mut self.syn_temp_noise_kind, 2, "Drift");
                        });
                    match self.syn_temp_noise_kind {
                        0 => {
                            ui.horizontal(|ui| {
                                ui.label("sigma");
                                ui.text_edit_singleline(&mut self.syn_temp_sigma);
                            });
                        }
                        1 => {
                            ui.horizontal(|ui| {
                                ui.label("amplitude");
                                ui.text_edit_singleline(&mut self.syn_temp_ampl);
                            });
                        }
                        _ => {
                            ui.horizontal(|ui| {
                                ui.label("slope");
                                ui.text_edit_singleline(&mut self.syn_temp_slope);
                            });
                        }
                    }
                    ui.separator();
                    ui.heading("Spikes");
                    ui.checkbox(&mut self.syn_spikes_enabled, "enabled");
                    ui.horizontal(|ui| {
                        ui.label("probability");
                        ui.text_edit_singleline(&mut self.syn_spike_prob);
                        ui.label("amplitude");
                        ui.text_edit_singleline(&mut self.syn_spike_ampl);
                    });
                });

                if let Some(status) = &self.syn_status {
                    ui.colored_label(egui::Color32::YELLOW, status);
                }

                ui.separator();
                ui.horizontal(|ui| {
                    if ui.button("Generate synthetic data").clicked() {
                        match self.try_build_config() {
                            Ok(cfg) => {
                                match model.generate_synthetic_data_from_config(&cfg) {
                                    Ok(()) => {
                                        self.syn_status =
                                            Some("Synthetic data generated".to_string());
                                        // keep window open to allow further runs
                                    }
                                    Err(e) => {
                                        self.syn_status =
                                            Some(format!("Error generating data: {:?}", e));
                                    }
                                }
                            }
                            Err(err) => {
                                self.syn_status = Some(err);
                            }
                        }
                    }
                    if ui.button("Close").clicked() {
                        self.show_synthetic_data = false;
                    }
                });
            });
        self.show_synthetic_data = open;
    }

    fn parse_csv_f64(s: &str) -> Result<Vec<f64>, String> {
        let t = s.trim();
        if t.is_empty() {
            return Ok(Vec::new());
        }
        t.split([',', ';', ' '])
            .filter(|p| !p.trim().is_empty())
            .map(|p| {
                p.trim()
                    .parse::<f64>()
                    .map_err(|_| format!("Invalid float: {}", p))
            })
            .collect()
    }

    fn parse_bool_opt_noise(
        kind: usize,
        sigma: &str,
        ampl: &str,
        slope: &str,
    ) -> Result<NoiseConfig, String> {
        let k = match kind {
            0 => {
                let s = sigma
                    .trim()
                    .parse::<f64>()
                    .map_err(|_| format!("Invalid sigma: {}", sigma))?;
                NoiseKind::Gaussian { sigma: s }
            }
            1 => {
                let a = ampl
                    .trim()
                    .parse::<f64>()
                    .map_err(|_| format!("Invalid amplitude: {}", ampl))?;
                NoiseKind::Uniform { amplitude: a }
            }
            2 => {
                let sl = slope
                    .trim()
                    .parse::<f64>()
                    .map_err(|_| format!("Invalid slope: {}", slope))?;
                NoiseKind::Drift { slope: sl }
            }
            _ => return Err("Unknown noise kind".to_string()),
        };
        Ok(NoiseConfig { kind: k })
    }

    fn try_build_config(&self) -> Result<AdvancedTGAConfig, String> {
        let n_points = self
            .syn_n_points
            .trim()
            .parse::<usize>()
            .map_err(|_| format!("Invalid n_points: {}", self.syn_n_points))?;
        let dt = self
            .syn_dt
            .trim()
            .parse::<f64>()
            .map_err(|_| format!("Invalid dt: {}", self.syn_dt))?;
        let seed = self
            .syn_seed
            .trim()
            .parse::<u64>()
            .map_err(|_| format!("Invalid seed: {}", self.syn_seed))?;

        let experiment_mode = match self.syn_mode {
            0 => {
                let temps = Self::parse_csv_f64(&self.syn_iso_temperatures)?;
                if temps.is_empty() {
                    return Err("Isothermal mode requires non-empty temperatures".to_string());
                }
                ExperimentMode::Isothermal {
                    temperatures: temps,
                }
            }
            _ => {
                let t0 = self
                    .syn_ni_t0
                    .trim()
                    .parse::<f64>()
                    .map_err(|_| format!("Invalid t0: {}", self.syn_ni_t0))?;
                let betas = Self::parse_csv_f64(&self.syn_ni_heating_rates)?;
                if betas.is_empty() {
                    return Err("NonIsothermal mode requires non-empty heating_rates".to_string());
                }
                ExperimentMode::NonIsothermal {
                    t0,
                    heating_rates: betas,
                }
            }
        };

        let kinetic_model = match self.syn_kmodel {
            0 => {
                let m0 = self
                    .syn_m0
                    .trim()
                    .parse::<f64>()
                    .map_err(|_| format!("Invalid m0: {}", self.syn_m0))?;
                let k0 = self
                    .syn_k0
                    .trim()
                    .parse::<f64>()
                    .map_err(|_| format!("Invalid k0: {}", self.syn_k0))?;
                let e = self
                    .syn_e
                    .trim()
                    .parse::<f64>()
                    .map_err(|_| format!("Invalid E: {}", self.syn_e))?;
                let r = self
                    .syn_r
                    .trim()
                    .parse::<f64>()
                    .map_err(|_| format!("Invalid R: {}", self.syn_r))?;
                KineticModel::ArrheniusSingle { m0, k0, e, r }
            }
            _ => {
                let m01 = self
                    .syn_m01
                    .trim()
                    .parse::<f64>()
                    .map_err(|_| format!("Invalid m01: {}", self.syn_m01))?;
                let k01 = self
                    .syn_k01
                    .trim()
                    .parse::<f64>()
                    .map_err(|_| format!("Invalid k01: {}", self.syn_k01))?;
                let e1 = self
                    .syn_e1
                    .trim()
                    .parse::<f64>()
                    .map_err(|_| format!("Invalid e1: {}", self.syn_e1))?;
                let m02 = self
                    .syn_m02
                    .trim()
                    .parse::<f64>()
                    .map_err(|_| format!("Invalid m02: {}", self.syn_m02))?;
                let k02 = self
                    .syn_k02
                    .trim()
                    .parse::<f64>()
                    .map_err(|_| format!("Invalid k02: {}", self.syn_k02))?;
                let e2 = self
                    .syn_e2
                    .trim()
                    .parse::<f64>()
                    .map_err(|_| format!("Invalid e2: {}", self.syn_e2))?;
                let r = self
                    .syn_r2
                    .trim()
                    .parse::<f64>()
                    .map_err(|_| format!("Invalid R: {}", self.syn_r2))?;
                KineticModel::ArrheniusTwoComponent {
                    m01,
                    k01,
                    e1,
                    m02,
                    k02,
                    e2,
                    r,
                }
            }
        };

        let mass_noise = if self.syn_mass_noise_enabled {
            Some(Self::parse_bool_opt_noise(
                self.syn_mass_noise_kind,
                &self.syn_mass_sigma,
                &self.syn_mass_ampl,
                &self.syn_mass_slope,
            )?)
        } else {
            None
        };

        let temp_noise = if self.syn_temp_noise_enabled {
            Some(Self::parse_bool_opt_noise(
                self.syn_temp_noise_kind,
                &self.syn_temp_sigma,
                &self.syn_temp_ampl,
                &self.syn_temp_slope,
            )?)
        } else {
            None
        };

        let spikes = if self.syn_spikes_enabled {
            let p = self
                .syn_spike_prob
                .trim()
                .parse::<f64>()
                .map_err(|_| format!("Invalid spike probability: {}", self.syn_spike_prob))?;
            let a = self
                .syn_spike_ampl
                .trim()
                .parse::<f64>()
                .map_err(|_| format!("Invalid spike amplitude: {}", self.syn_spike_ampl))?;
            Some(SpikeModel {
                probability: p,
                amplitude: a,
            })
        } else {
            None
        };

        Ok(AdvancedTGAConfig {
            n_points,
            dt,
            experiment_mode,
            kinetic_model,
            mass_noise,
            temp_noise,
            spikes,
            seed,
        })
    }

    /// Loads current settings into temporary edit fields
    fn load_current_settings(&mut self, settings: &Settings) {
        if let Some(cal) = settings.calibration_line() {
            self.temp_k = cal.k().to_string();
            self.temp_b = cal.b().to_string();
        } else {
            self.temp_k = "40.0".to_string();
            self.temp_b = "-1.0".to_string();
        }

        self.temp_n_points = settings
            .n_points()
            .map(|n| n.to_string())
            .unwrap_or_else(|| "1000".to_string());

        self.temp_log_level = settings.log_level().unwrap_or("info").to_string();
    }

    /// Applies edited settings to the model and saves to file
    fn apply_settings(&self, model: &mut PlotModel) -> Result<(), String> {
        // Parse calibration coefficients
        let k = self
            .temp_k
            .parse::<f64>()
            .map_err(|_| format!("Invalid k value: '{}'", self.temp_k))?;
        let b = self
            .temp_b
            .parse::<f64>()
            .map_err(|_| format!("Invalid b value: '{}'", self.temp_b))?;

        // Parse number of points
        let n_points = self
            .temp_n_points
            .parse::<usize>()
            .map_err(|_| format!("Invalid n_points value: '{}'", self.temp_n_points))?;

        // Update model settings
        model
            .settings
            .set_calibration_line(Some(CalibrationLine::new(k, b)));
        model.settings.set_n_points(Some(n_points));
        model
            .settings
            .set_log_level(Some(self.temp_log_level.clone()));

        // Save to file
        model
            .settings
            .create_or_recreate_config()
            .map_err(|e| format!("Failed to save settings: {:?}", e))?;

        Ok(())
    }

    /// Returns whether the settings window is currently shown
    pub fn is_settings_window_open(&self) -> bool {
        self.show_settings_window
    }

    pub fn show_help_window_ui(&mut self, ctx: &egui::Context) {
        if !self.show_help_window {
            return;
        }
        let mut open = self.show_help_window;
        egui::Window::new("Help")
            .open(&mut open)
            .resizable(true)
            .vscroll(true)
            .default_width(640.0)
            .default_height(520.0)
            .show(ctx, |ui| {
                ui.horizontal(|ui| {
                    ui.label("Language:");
                    let mut lang = self.help_lang;
                    if ui
                        .radio_value(&mut lang, HelpLanguage::Eng, "ENG")
                        .clicked()
                    {
                        self.help_lang = HelpLanguage::Eng;
                        let _ = self.reload_help();
                    }
                    if ui
                        .radio_value(&mut lang, HelpLanguage::Rus, "RUS")
                        .clicked()
                    {
                        self.help_lang = HelpLanguage::Rus;
                        let _ = self.reload_help();
                    }
                    ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                        if ui.button("Close").clicked() {
                            self.show_help_window = false;
                        }
                    });
                });
                if let Some(err) = &self.help_error {
                    ui.colored_label(egui::Color32::RED, err);
                }
                ui.separator();
                egui::ScrollArea::vertical().show(ui, |ui| {
                    for (si, sec) in self.parsed_help.sections.iter().enumerate() {
                        let id = egui::Id::new(("help_sec", si));
                        egui::CollapsingHeader::new(&sec.title)
                            .id_salt(id)
                            .show(ui, |ui| {
                                for (pi, par) in sec.paragraphs.iter().enumerate() {
                                    let pid = egui::Id::new(("help_par", si, pi));
                                    egui::CollapsingHeader::new(&par.title).id_salt(pid).show(
                                        ui,
                                        |ui| {
                                            ui.label(
                                                egui::RichText::new(par.body.clone()).monospace(),
                                            );
                                        },
                                    );
                                }
                            });
                    }
                });
            });
        self.show_help_window = open;
    }

    fn reload_help(&mut self) -> Result<(), String> {
        let path = match self.help_lang {
            HelpLanguage::Eng => "assets/help_eng.txt",
            HelpLanguage::Rus => "assets/help_rus.txt",
        };
        match std::fs::read_to_string(path) {
            Ok(contents) => {
                self.parsed_help = parse_help(&contents);
                self.help_error = None;
                Ok(())
            }
            Err(e) => {
                self.parsed_help = HelpDoc::default();
                self.help_error = Some(format!("Failed to read {}: {}", path, e));
                Err(format!("{}", e))
            }
        }
    }
}

fn parse_help(src: &str) -> HelpDoc {
    let mut doc = HelpDoc {
        sections: Vec::new(),
    };
    let mut current_section: Option<HelpSection> = None;
    let mut current_paragraph: Option<HelpParagraph> = None;
    let finalize_paragraph = |sec: &mut Option<HelpSection>, par: &mut Option<HelpParagraph>| {
        if let (Some(s), Some(p)) = (sec.as_mut(), par.take()) {
            s.paragraphs.push(p);
        }
    };
    let finalize_section =
        |doc: &mut HelpDoc, sec: &mut Option<HelpSection>, par: &mut Option<HelpParagraph>| {
            finalize_paragraph(sec, par);
            if let Some(s) = sec.take() {
                doc.sections.push(s);
            }
        };
    for line in src.lines() {
        let trimmed = line.trim_end_matches(['\r']);
        if trimmed.starts_with("## ") {
            finalize_paragraph(&mut current_section, &mut current_paragraph);
            current_paragraph = Some(HelpParagraph {
                title: trimmed[3..].to_string(),
                body: String::new(),
            });
        } else if trimmed.starts_with("# ") {
            finalize_section(&mut doc, &mut current_section, &mut current_paragraph);
            current_section = Some(HelpSection {
                title: trimmed[2..].to_string(),
                paragraphs: Vec::new(),
            });
        } else {
            if let Some(p) = current_paragraph.as_mut() {
                if !p.body.is_empty() {
                    p.body.push('\n');
                }
                p.body.push_str(trimmed);
            }
        }
    }
    finalize_section(&mut doc, &mut current_section, &mut current_paragraph);
    doc
}

//============================================================================================
//
//============================================================================================
impl PlotModel {
    fn push_testing_data(&mut self) {
        *self = load_model_with_one_curve();
    }
}
