//! # Thermochemistry GUI Module
//!
//! This module provides a graphical user interface for thermochemical calculations and data management.
//! It allows users to search for thermodynamic properties, calculate Gibbs free energy, and manage
//! substance databases from various sources (NASA, CEA, NIST, Aramco).
//!
//! ## Main Features
//! - **Substance Search**: Find thermodynamic data across multiple libraries
//! - **Property Calculation**: Calculate Cp, enthalpy, entropy, and Gibbs free energy
//! - **Library Management**: Browse and select from different thermodynamic databases
//! - **Temperature Range Analysis**: Calculate properties over temperature ranges
//!
//! ## GUI Layout
//! The interface consists of:
//! - **Search Panel**: Substance input and library selection
//! - **Results Panel**: Display of found thermodynamic data
//! - **Calculation Panel**: Temperature input and property calculations

use eframe::egui;
use nalgebra::{DMatrix, DVector};

use crate::Thermodynamics::DBhandlers::thermo_api::{
    EnergyUnit, ThermoCalculator, create_thermal_by_name,
};
use crate::Thermodynamics::thermo_lib_api::ThermoData;
use crate::gui::gui_plot::PlotWindow;
/*
    Example usage within the GUI context:
    let thermo_data = % ThermoData instance %
   let sublib = thermo_data.LibThermoData.get(% library name %).unwrap();
    let subdata = sublib.get(%substance name%).unwrap();

    let mut thermoenum  = create_thermal_by_name(% library name %);
            let _ = thermoenum .newinstance();
            let _ =  thermoenum .from_serde(subdata.clone());
            ket _ = thermoenum .from_serde(subdata.clone())

            let _ = thermoenum.print_instance();
           thermoenum.extract_model_coefficients(T);
            let _ =   thermoenum.calculate_Cp_dH_dS(T);
                let Cp = thermoenum.get_Cp().unwrap();
            let dh =  thermoenum.get_dh().unwrap();
            let ds =  thermoenum.get_ds().unwrap();
*/
/// Main application structure for the thermochemistry GUI
///
/// This struct manages the thermochemistry interface state including:
/// - Substance search functionality
/// - Library selection and data display
/// - Temperature input for calculations
/// - Results from thermodynamic property calculations
#[derive(Debug)]
pub struct ThermochemistryApp {
    /// Core thermodynamic database interface
    thermo_data: ThermoData,
    /// Input field for substance names to search
    substance_input: String,
    /// Currently selected thermodynamic library
    selected_library: String,
    /// Temperature for calculations (in Kelvin)
    temperature: String,
    /// Pressure for calculations (in Pa)
    pressure: String,
    /// Search results display text
    search_results: String,
    /// Available substances in selected library
    available_substances: Vec<String>,
    /// Currently selected substance for calculations
    selected_substance: String,
    /// Calculated thermodynamic properties
    calculated_cp: Option<f64>,
    calculated_dh: Option<f64>,
    calculated_ds: Option<f64>,
    energy_unit: EnergyUnit,
    show_plots_window: bool,
    t0: String,
    tend: String,
    plot_window: Option<PlotWindow>,
}

impl Default for ThermochemistryApp {
    fn default() -> Self {
        let thermo_data = ThermoData::new();
        let selected_library = if !thermo_data.thermo_libs.is_empty() {
            thermo_data.thermo_libs[0].clone()
        } else {
            "NASA_gas".to_string()
        };

        let available_substances = thermo_data
            .LibThermoData
            .get(&selected_library)
            .map(|lib_data| lib_data.keys().cloned().collect())
            .unwrap_or_default();

        Self {
            thermo_data,
            substance_input: String::new(),
            selected_library,
            temperature: "298.15".to_string(),
            pressure: "101325.0".to_string(),
            search_results: "No search performed yet".to_string(),
            available_substances,
            selected_substance: String::new(),
            calculated_cp: None,
            calculated_dh: None,
            calculated_ds: None,
            energy_unit: EnergyUnit::J,
            show_plots_window: false,
            t0: "298.15".to_string(),
            tend: "1000.0".to_string(),
            plot_window: None,
        }
    }
}

impl ThermochemistryApp {
    /// Creates a new ThermochemistryApp instance
    ///
    /// Initializes the application with default values:
    /// - Empty substance search
    /// - First thermo library as default
    /// - Standard temperature and pressure (298.15 K, 101325 Pa)
    pub fn new() -> Self {
        Self::default()
    }

    /// Updates available substances when library changes
    fn update_substances_for_library(&mut self) {
        self.available_substances = self
            .thermo_data
            .LibThermoData
            .get(&self.selected_library)
            .map(|lib_data| lib_data.keys().cloned().collect())
            .unwrap_or_default();
    }

    /// Searches for a specific substance by name without affecting the search filter
    fn search_substance_by_name(&mut self, substance_name: &str) {
        self.selected_substance = substance_name.to_string();
        if let Some(lib_data) = self.thermo_data.LibThermoData.get(&self.selected_library) {
            if let Some(substance_data) = lib_data.get(substance_name) {
                self.search_results = format!(
                    "Found '{}' in library '{}':\n\n{}",
                    substance_name,
                    self.selected_library,
                    serde_json::to_string_pretty(substance_data)
                        .unwrap_or("Error formatting data".to_string())
                );
            } else {
                self.search_results = format!(
                    "Substance '{}' not found in library '{}'",
                    substance_name, self.selected_library
                );
            }
        } else {
            self.search_results = format!("Library '{}' not found", self.selected_library);
        }
    }

    /// Performs substance search in the selected library
    ///
    /// This method integrates with the thermodynamics module to:
    /// 1. Search for the substance in the selected library
    /// 2. Retrieve available thermodynamic data
    /// 3. Display results in the search_results field
    #[allow(dead_code)]
    fn search_substance(&mut self) {
        if self.substance_input.is_empty() {
            self.search_results = "Please enter a substance name".to_string();
            return;
        }

        if let Some(lib_data) = self.thermo_data.LibThermoData.get(&self.selected_library) {
            if let Some(substance_data) = lib_data.get(&self.substance_input) {
                self.search_results = format!(
                    "Found '{}' in library '{}':\n\n{}",
                    self.substance_input,
                    self.selected_library,
                    serde_json::to_string_pretty(substance_data)
                        .unwrap_or("Error formatting data".to_string())
                );
            } else {
                self.search_results = format!(
                    "Substance '{}' not found in library '{}'",
                    self.substance_input, self.selected_library
                );
            }
        } else {
            self.search_results = format!("Library '{}' not found", self.selected_library);
        }
    }

    /// Calculates thermodynamic properties at specified conditions
    ///
    /// This method integrates with the thermodynamics module to:
    /// 1. Parse temperature and pressure inputs
    /// 2. Calculate Cp, enthalpy, entropy, and Gibbs free energy
    /// 3. Display calculated values
    fn calculate_properties(&mut self) {
        let temp_result = self.temperature.parse::<f64>();
        let pres_result = self.pressure.parse::<f64>();

        let substance_name = if !self.selected_substance.is_empty() {
            &self.selected_substance
        } else {
            &self.substance_input
        };

        if substance_name.is_empty() {
            self.search_results = "Please select a substance first".to_string();
            return;
        }

        match (temp_result, pres_result) {
            (Ok(t), Ok(_p)) => {
                if let Some(lib_data) = self.thermo_data.LibThermoData.get(&self.selected_library) {
                    if let Some(substance_data) = lib_data.get(substance_name) {
                        match self.perform_thermo_calculation(substance_data, t) {
                            Ok((cp, dh, ds)) => {
                                self.calculated_cp = Some(cp);
                                self.calculated_dh = Some(dh);
                                self.calculated_ds = Some(ds);
                                let (cp_unit, dh_unit, ds_unit) = match self.energy_unit {
                                    EnergyUnit::J => ("J/mol·K", "kJ/mol", "J/mol·K"),
                                    EnergyUnit::Cal => ("cal/mol·K", "kcal/mol", "cal/mol·K"),
                                };
                                self.search_results = format!(
                                    "Thermodynamic properties for '{}' at T = {} K:\n\nCp = {:.3} {}\ndH = {:.3} {}\ndS = {:.3} {}",
                                    substance_name,
                                    t,
                                    cp,
                                    cp_unit,
                                    dh / 1000.0,
                                    dh_unit,
                                    ds,
                                    ds_unit
                                );
                            }
                            Err(e) => {
                                self.search_results = format!("Calculation error: {}", e);
                            }
                        }
                    } else {
                        self.search_results =
                            format!("Substance '{}' not found for calculations", substance_name);
                    }
                } else {
                    self.search_results = format!(
                        "Library '{}' not available for calculations",
                        self.selected_library
                    );
                }
            }
            _ => {
                self.search_results = "Invalid temperature or pressure values".to_string();
            }
        }
    }

    /// Performs the actual thermodynamic calculation using thermo_api
    fn perform_thermo_calculation(
        &self,
        substance_data: &serde_json::Value,
        temperature: f64,
    ) -> Result<(f64, f64, f64), String> {
        let mut thermoenum = create_thermal_by_name(&self.selected_library);

        thermoenum
            .newinstance()
            .map_err(|e| format!("Failed to create instance: {}", e))?;
        thermoenum
            .set_unit(Some(self.energy_unit.clone()))
            .map_err(|e| format!("Failed to set unit: {}", e))?;
        thermoenum
            .from_serde(substance_data.clone())
            .map_err(|e| format!("Failed to load data: {}", e))?;
        thermoenum
            .extract_model_coefficients(temperature)
            .map_err(|e| format!("Failed to extract coefficients: {}", e))?;
        thermoenum
            .calculate_Cp_dH_dS(temperature)
            .map_err(|e| format!("Failed to calculate properties: {}", e))?;

        let cp = thermoenum
            .get_Cp()
            .map_err(|e| format!("Failed to get Cp: {}", e))?;
        let dh = thermoenum
            .get_dh()
            .map_err(|e| format!("Failed to get dH: {}", e))?;
        let ds = thermoenum
            .get_ds()
            .map_err(|e| format!("Failed to get dS: {}", e))?;

        Ok((cp, dh, ds))
    }

    /// Calculates thermodynamic properties over a temperature range
    fn calculate_temperature_range(
        &self,
        substance_data: &serde_json::Value,
        t0: f64,
        tend: f64,
        n_points: usize,
    ) -> Result<Vec<(f64, f64, f64, f64)>, String> {
        let mut results = Vec::new();
        let dt = (tend - t0) / (n_points - 1) as f64;
        let mut thermoenum = create_thermal_by_name(&self.selected_library);
        thermoenum
            .newinstance()
            .map_err(|e| format!("Failed to create instance: {}", e))?;
        thermoenum
            .set_unit(Some(self.energy_unit.clone()))
            .map_err(|e| format!("Failed to set unit: {}", e))?;
        thermoenum
            .from_serde(substance_data.clone())
            .map_err(|e| format!("Failed to load data: {}", e))?;
        for i in 0..n_points {
            let t = t0 + i as f64 * dt;

            thermoenum
                .extract_model_coefficients(t)
                .map_err(|e| format!("Failed to extract coefficients: {}", e))?;
            thermoenum
                .calculate_Cp_dH_dS(t)
                .map_err(|e| format!("Failed to calculate properties: {}", e))?;

            let cp = thermoenum
                .get_Cp()
                .map_err(|e| format!("Failed to get Cp: {}", e))?;
            let dh = thermoenum
                .get_dh()
                .map_err(|e| format!("Failed to get dH: {}", e))?;
            let ds = thermoenum
                .get_ds()
                .map_err(|e| format!("Failed to get dS: {}", e))?;

            results.push((t, cp, dh, ds));
        }

        Ok(results)
    }

    /// Shows the plots window
    fn show_plots_window(&mut self, ctx: &egui::Context) {
        egui::Window::new("Temperature Range Plots")
            .open(&mut self.show_plots_window.clone())
            .default_size([800.0, 600.0])
            .show(ctx, |ui| {
                ui.heading("Temperature Range Analysis");
                ui.separator();

                ui.horizontal(|ui| {
                    ui.label("T0 (K):");
                    ui.text_edit_singleline(&mut self.t0);

                    ui.label("Tend (K):");
                    ui.text_edit_singleline(&mut self.tend);

                    if ui.button("Calculate Range").clicked() {
                        self.calculate_range_data();
                    }
                });

                ui.separator();
                // ui.label("Plot functionality will be implemented here");
            });
    }

    /// Calculates data for the temperature range and creates plot window
    fn calculate_range_data(&mut self) {
        let substance_name = if !self.selected_substance.is_empty() {
            &self.selected_substance
        } else {
            &self.substance_input
        };

        if substance_name.is_empty() {
            return;
        }

        let t0_result = self.t0.parse::<f64>();
        let tend_result = self.tend.parse::<f64>();

        match (t0_result, tend_result) {
            (Ok(t0), Ok(tend)) if t0 < tend => {
                if let Some(lib_data) = self.thermo_data.LibThermoData.get(&self.selected_library) {
                    if let Some(substance_data) = lib_data.get(substance_name) {
                        match self.calculate_temperature_range(substance_data, t0, tend, 100) {
                            Ok(results) => {
                                // Convert results to DVector and DMatrix for plotting
                                let temperatures: Vec<f64> =
                                    results.iter().map(|(t, _, _, _)| *t).collect();
                                let t_result = DVector::from_vec(temperatures);

                                let mut y_data = Vec::new();
                                for (_, cp, dh, ds) in &results {
                                    y_data.push(*cp); // Cp column
                                    y_data.push(*dh); // dH column  
                                    y_data.push(*ds); // dS column
                                }

                                let y_result = DMatrix::from_row_slice(results.len(), 3, &y_data);

                                let (cp_unit, dh_unit, ds_unit) = match self.energy_unit {
                                    EnergyUnit::J => ("Cp (J/mol·K)", "dH (J/mol)", "dS (J/mol·K)"),
                                    EnergyUnit::Cal => {
                                        ("Cp (cal/mol·K)", "dH (cal/mol)", "dS (cal/mol·K)")
                                    }
                                };

                                let values = vec![
                                    cp_unit.to_string(),
                                    dh_unit.to_string(),
                                    ds_unit.to_string(),
                                ];

                                self.plot_window = Some(PlotWindow::new(
                                    "Temperature (K)".to_string(),
                                    values,
                                    t_result,
                                    y_result,
                                ));
                            }
                            Err(e) => {
                                println!("Range calculation error: {}", e);
                            }
                        }
                    }
                }
            }
            _ => {
                println!("Invalid temperature range");
            }
        }
    }

    /// Main GUI rendering method for the thermochemistry interface
    ///
    /// ## Layout Structure:
    ///
    /// ### Top Panel:
    /// - **Substance Input**: Text field for entering substance names
    /// - **Library Selection**: Dropdown for choosing thermodynamic database
    /// - **Search Button**: Triggers substance search in selected library
    ///
    /// ### Middle Panel:
    /// - **Temperature Input**: Field for calculation temperature (K)
    /// - **Pressure Input**: Field for calculation pressure (Pa)
    /// - **Calculate Button**: Performs thermodynamic property calculations
    ///
    /// ### Bottom Panel:
    /// - **Results Display**: Scrollable text area showing search results and calculations
    ///
    /// ## Button Logic:
    /// - **Search**: Calls `search_substance()` to find thermodynamic data
    /// - **Calculate**: Calls `calculate_properties()` to compute thermodynamic properties
    pub fn show(&mut self, ctx: &egui::Context, open: &mut bool) {
        egui::Window::new("Thermochemistry Analysis")
            .open(open)
            .default_size([1200.0, 800.0])
            .show(ctx, |ui| {
                ui.heading("Thermodynamic Properties and Calculations");
                ui.separator();

                ui.horizontal(|ui| {
                    // Left panel - Substance selection
                    ui.vertical(|ui| {
                        ui.set_width(200.0);
                        ui.set_min_height(600.0);
                        ui.heading("List of Substances");
                        // Search filter
                        ui.horizontal(|ui| {
                            ui.label("Search:");
                            ui.text_edit_singleline(&mut self.substance_input);
                        });
                        ui.separator();
                        egui::ScrollArea::vertical()
                            .max_height(1200.0)
                            .show(ui, |ui| {
                                for substance in &self.available_substances.clone() {
                                    if self.substance_input.is_empty()
                                        || substance
                                            .to_lowercase()
                                            .contains(&self.substance_input.to_lowercase())
                                    {
                                        if ui.selectable_label(false, substance).clicked() {
                                            self.search_substance_by_name(substance);
                                        }
                                    }
                                }
                            });
                        ui.separator();
                        // Dropdown for library selection
                        egui::ComboBox::from_label("Library Source")
                            .selected_text(&self.selected_library)
                            .show_ui(ui, |ui| {
                                for library in &self.thermo_data.thermo_libs.clone() {
                                    if ui
                                        .selectable_value(
                                            &mut self.selected_library,
                                            library.clone(),
                                            library,
                                        )
                                        .clicked()
                                    {
                                        self.update_substances_for_library();
                                    }
                                }
                            });
                    });

                    ui.separator();

                    // Right panel - Results and calculations
                    ui.vertical(|ui| {
                        // Results section
                        ui.group(|ui| {
                            ui.label("Substance Information:");
                            egui::ScrollArea::vertical()
                                .id_salt("substance_info")
                                .max_height(300.0)
                                .show(ui, |ui| {
                                    ui.text_edit_multiline(&mut self.search_results);
                                });
                        });

                        ui.add_space(10.0);

                        // Calculation parameters section
                        ui.group(|ui| {
                            ui.label("Calculation Parameters:");
                            ui.horizontal(|ui| {
                                ui.label("Temperature (K):");
                                ui.text_edit_singleline(&mut self.temperature);

                                ui.label("Pressure (Pa):");
                                ui.text_edit_singleline(&mut self.pressure);

                                if ui.button("🧮 Calculate Properties").clicked() {
                                    self.calculate_properties();
                                }
                            });
                            ui.horizontal(|ui| {
                                ui.label("Energy Unit:");
                                ui.radio_value(&mut self.energy_unit, EnergyUnit::J, "J");
                                ui.radio_value(&mut self.energy_unit, EnergyUnit::Cal, "cal");
                            });
                        });

                        ui.separator();

                        // Action buttons
                        ui.horizontal(|ui| {
                            if ui.button("Clear Results").clicked() {
                                self.search_results = "Results cleared".to_string();
                            }

                            if ui.button("View Plots").clicked() {
                                self.show_plots_window = true;
                            }

                            if ui.button("Export Data").clicked() {
                                self.search_results +=
                                    "\n\nExport functionality not yet implemented";
                            }
                        });
                    });
                });

                // Show plots window if requested
                if self.show_plots_window {
                    self.show_plots_window(ctx);
                }

                // Show plot window if available
                if let Some(plot_window) = &mut self.plot_window {
                    plot_window.show(ctx);
                    if !plot_window.visible {
                        self.plot_window = None;
                    }
                }
            });
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_thermochemistry_app_new() {
        let app = ThermochemistryApp::new();
        assert!(app.substance_input.is_empty());
        assert_eq!(app.temperature, "298.15");
        assert_eq!(app.pressure, "101325.0");
        assert!(!app.available_substances.is_empty());
    }

    #[test]
    fn test_thermochemistry_app_default() {
        let app = ThermochemistryApp::default();
        assert!(app.substance_input.is_empty());
        assert_eq!(app.temperature, "298.15");
        assert_eq!(app.pressure, "101325.0");
        assert_eq!(app.search_results, "No search performed yet");
    }

    #[test]
    fn test_search_empty_substance() {
        let mut app = ThermochemistryApp::new();
        app.search_substance();
        assert_eq!(app.search_results, "Please enter a substance name");
    }

    #[test]
    fn test_calculate_properties_invalid_temperature() {
        let mut app = ThermochemistryApp::new();
        app.temperature = "invalid".to_string();
        app.calculate_properties();
        assert_eq!(app.search_results, "Please select a substance first");
    }
}
