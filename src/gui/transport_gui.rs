//! # Transport Properties GUI Module
//!
//! This module provides a graphical user interface for transport property calculations.
//! It allows users to search for transport data and calculate viscosity and thermal conductivity
//! from various sources (CEA, Aramco transport).
//!
//! ## Main Features
//! - **Substance Search**: Find transport data across multiple libraries
//! - **Property Calculation**: Calculate viscosity and thermal conductivity
//! - **Library Management**: Browse and select from different transport databases
//! - **Unit Management**: Support for different unit systems
//!
//! ## GUI Layout
//! The interface consists of:
//! - **Search Panel**: Substance input and library selection
//! - **Results Panel**: Display of found transport data
//! - **Calculation Panel**: Temperature input and transport property calculations

use eframe::egui;

use crate::Thermodynamics::DBhandlers::transport_api::{
    LambdaUnit, TransportCalculator, ViscosityUnit, create_transport_calculator_by_name,
};
use crate::Thermodynamics::thermo_lib_api::ThermoData;
/*
    Example usage within the GUI context:
    let transport_data = ThermoData::new();
    let sublib = transport_data.LibThermoData.get(library_name).unwrap();
    let subdata = sublib.get(substance_name).unwrap();

    let mut transport_calc = create_transport_calculator_by_name(library_name);
    transport_calc.from_serde(subdata.clone())?;
    transport_calc.set_lambda_unit(Some(LambdaUnit::MWPerMK))?;
    transport_calc.set_viscosity_unit(Some(ViscosityUnit::MKPaS))?;
    transport_calc.extract_coefficients(T)?;
    let lambda = transport_calc.calculate_lambda(Some(cp), Some(density), T)?;
    let viscosity = transport_calc.calculate_viscosity(T)?;
*/
/// Main application structure for the transport properties GUI
///
/// This struct manages the transport interface state including:
/// - Substance search functionality
/// - Library selection and data display
/// - Temperature input for calculations
/// - Results from transport property calculations
#[derive(Debug)]
pub struct TransportApp {
    /// Core transport database interface
    transport_data: ThermoData,
    /// Input field for substance names to search
    substance_input: String,
    /// Currently selected transport library
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
    /// Calculated transport properties
    calculated_lambda: Option<f64>,
    calculated_viscosity: Option<f64>,
    lambda_unit: LambdaUnit,
    viscosity_unit: ViscosityUnit,
}

impl Default for TransportApp {
    fn default() -> Self {
        let transport_data = ThermoData::new();
        let selected_library = if !transport_data.transport_libs.is_empty() {
            transport_data.transport_libs[0].clone()
        } else {
            "Aramco_transport".to_string()
        };

        let available_substances = transport_data
            .LibThermoData
            .get(&selected_library)
            .map(|lib_data| lib_data.keys().cloned().collect())
            .unwrap_or_default();

        Self {
            transport_data,
            substance_input: String::new(),
            selected_library,
            temperature: "298.15".to_string(),
            pressure: "101325.0".to_string(),
            search_results: "No search performed yet".to_string(),
            available_substances,
            selected_substance: String::new(),
            calculated_lambda: None,
            calculated_viscosity: None,
            lambda_unit: LambdaUnit::MKWPerMK,
            viscosity_unit: ViscosityUnit::MKPaS,
        }
    }
}

impl TransportApp {
    /// Creates a new TransportApp instance
    ///
    /// Initializes the application with default values:
    /// - Empty substance search
    /// - First transport library as default
    /// - Standard temperature and pressure (298.15 K, 101325 Pa)
    pub fn new() -> Self {
        Self::default()
    }

    /// Updates available substances when library changes
    fn update_substances_for_library(&mut self) {
        self.available_substances = self
            .transport_data
            .LibThermoData
            .get(&self.selected_library)
            .map(|lib_data| lib_data.keys().cloned().collect())
            .unwrap_or_default();
    }

    /// Searches for a specific substance by name without affecting the search filter
    fn search_substance_by_name(&mut self, substance_name: &str) {
        self.selected_substance = substance_name.to_string();
        if let Some(lib_data) = self
            .transport_data
            .LibThermoData
            .get(&self.selected_library)
        {
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
    /// This method integrates with the transport module to:
    /// 1. Search for the substance in the selected library
    /// 2. Retrieve available transport data
    /// 3. Display results in the search_results field
    #[allow(dead_code)]
    fn search_substance(&mut self) {
        if self.substance_input.is_empty() {
            self.search_results = "Please enter a substance name".to_string();
            return;
        }

        if let Some(lib_data) = self
            .transport_data
            .LibThermoData
            .get(&self.selected_library)
        {
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

    /// Calculates transport properties at specified conditions
    ///
    /// This method integrates with the transport module to:
    /// 1. Parse temperature and pressure inputs
    /// 2. Calculate viscosity and thermal conductivity
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
                if let Some(lib_data) = self
                    .transport_data
                    .LibThermoData
                    .get(&self.selected_library)
                {
                    if let Some(substance_data) = lib_data.get(substance_name) {
                        match self.perform_transport_calculation(substance_data, t) {
                            Ok((lambda, viscosity)) => {
                                self.calculated_lambda = Some(lambda);
                                self.calculated_viscosity = Some(viscosity);
                                let lambda_unit_str = match self.lambda_unit {
                                    LambdaUnit::WPerMK => "W/m·K",
                                    LambdaUnit::MWPerMK => "mW/m·K",
                                    LambdaUnit::MKWPerMK => "μW/m·K",
                                    LambdaUnit::MKWPerSMK => "μW/s·m·K",
                                };
                                let viscosity_unit_str = match self.viscosity_unit {
                                    ViscosityUnit::KgPerMS => "kg/m·s",
                                    ViscosityUnit::PaS => "Pa·s",
                                    ViscosityUnit::MKPaS => "μPa·s",
                                };
                                self.search_results = format!(
                                    "Transport properties for '{}' at T = {} K:\n\nThermal conductivity = {:.5} {}\nViscosity = {:.8} {}",
                                    substance_name,
                                    t,
                                    lambda,
                                    lambda_unit_str,
                                    viscosity,
                                    viscosity_unit_str
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

    /// Performs the actual transport calculation using transport_api
    fn perform_transport_calculation(
        &self,
        substance_data: &serde_json::Value,
        temperature: f64,
    ) -> Result<(f64, f64), String> {
        let mut transport_calc = create_transport_calculator_by_name(&self.selected_library);

        transport_calc
            .from_serde(substance_data.clone())
            .map_err(|e| format!("Failed to load data: {}", e))?;
        transport_calc
            .set_lambda_unit(Some(self.lambda_unit))
            .map_err(|e| format!("Failed to set lambda unit: {}", e))?;
        transport_calc
            .set_viscosity_unit(Some(self.viscosity_unit))
            .map_err(|e| format!("Failed to set viscosity unit: {}", e))?;
        transport_calc
            .extract_coefficients(temperature)
            .map_err(|e| format!("Failed to extract coefficients: {}", e))?;
        let mut lambda: f64 = 0.0;
        if self.selected_library == "CEA" {
            lambda = transport_calc
                .calculate_lambda(None, None, temperature)
                .map_err(|e| format!("Failed to calculate thermal conductivity: {}", e))?;
        } else {
        };
        let viscosity = transport_calc
            .calculate_viscosity(temperature)
            .map_err(|e| format!("Failed to calculate viscosity: {}", e))?;

        Ok((lambda, viscosity))
    }

    /// Main GUI rendering method for the transport properties interface
    ///
    /// ## Layout Structure:
    ///
    /// ### Top Panel:
    /// - **Substance Input**: Text field for entering substance names
    /// - **Library Selection**: Dropdown for choosing transport database
    /// - **Search Button**: Triggers substance search in selected library
    ///
    /// ### Middle Panel:
    /// - **Temperature Input**: Field for calculation temperature (K)
    /// - **Pressure Input**: Field for calculation pressure (Pa)
    /// - **Calculate Button**: Performs transport property calculations
    ///
    /// ### Bottom Panel:
    /// - **Results Display**: Scrollable text area showing search results and calculations
    ///
    /// ## Button Logic:
    /// - **Search**: Calls `search_substance()` to find transport data
    /// - **Calculate**: Calls `calculate_properties()` to compute transport properties
    pub fn show(&mut self, ctx: &egui::Context, open: &mut bool) {
        egui::Window::new("Transport Properties Analysis")
            .open(open)
            .default_size([1200.0, 800.0])
            .show(ctx, |ui| {
                ui.heading("Transport Properties and Calculations");
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
                                for library in &self.transport_data.transport_libs.clone() {
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
                                ui.label("Lambda Unit:");
                                ui.radio_value(&mut self.lambda_unit, LambdaUnit::WPerMK, "W/m·K");
                                ui.radio_value(
                                    &mut self.lambda_unit,
                                    LambdaUnit::MWPerMK,
                                    "mW/m·K",
                                );
                                ui.radio_value(
                                    &mut self.lambda_unit,
                                    LambdaUnit::MKWPerMK,
                                    "μW/m·K",
                                );
                            });
                            ui.horizontal(|ui| {
                                ui.label("Viscosity Unit:");
                                ui.radio_value(
                                    &mut self.viscosity_unit,
                                    ViscosityUnit::PaS,
                                    "Pa·s",
                                );
                                ui.radio_value(
                                    &mut self.viscosity_unit,
                                    ViscosityUnit::MKPaS,
                                    "μPa·s",
                                );
                            });
                        });

                        ui.separator();

                        // Action buttons
                        ui.horizontal(|ui| {
                            if ui.button("Clear Results").clicked() {
                                self.search_results = "Results cleared".to_string();
                            }

                            if ui.button("Export Data").clicked() {
                                self.search_results +=
                                    "\n\nExport functionality not yet implemented";
                            }

                            if ui.button("Load from File").clicked() {
                                self.search_results +=
                                    "\n\nFile loading functionality not yet implemented";
                            }
                        });
                    });
                });
            });
    }
}
