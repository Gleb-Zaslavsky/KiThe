//! # Thermochemistry GUI Module
//!
//! This module provides a graphical user interface for thermochemical calculations and data management.
//! It allows users to search for thermodynamic properties, calculate Cp/enthalpy/entropy values,
//! and manage substance databases from thermochemistry sources such as NASA and NIST.
//!
//! ## Main Features
//! - **Substance Selection**: Inspect thermodynamic data across multiple libraries
//! - **Property Calculation**: Calculate Cp, enthalpy, and entropy
//! - **Library Management**: Browse and select from different thermodynamic databases
//! - **Temperature Range Analysis**: Calculate properties over temperature ranges
//!
//! ## GUI Layout
//! The interface consists of:
//! - **Selection Panel**: Substance input and library selection
//! - **Results Panel**: Display of found thermodynamic data
//! - **Calculation Panel**: Temperature input and property calculations

use eframe::egui;
use log::warn;
use nalgebra::{DMatrix, DVector};

use crate::Thermodynamics::DBhandlers::thermo_api::EnergyUnit;
use crate::Thermodynamics::User_substances::{LibraryPriority, SubsData};
use crate::Thermodynamics::thermo_lib_api::{LibraryId, ThermoData, ThermoLibraryError};
use crate::gui::NIST_gui::NISTApp;
use crate::gui::condition_parser::{parse_temperature, parse_temperature_range};
use crate::gui::gui_plot::PlotWindow;
use crate::gui::read_only_snapshot::multiline as read_only_multiline;
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
    pub(crate) thermo_data: ThermoData,
    /// Visible startup error if the catalog could not be loaded.
    pub(crate) startup_error: Option<String>,
    /// Input field for substance names to search
    pub(crate) substance_input: String,
    /// Currently selected thermodynamic library
    pub(crate) selected_library: String,
    /// Temperature for calculations (in Kelvin)
    pub(crate) temperature: String,
    /// Legacy pressure field kept only for compatibility with older state snapshots.
    pub(crate) pressure: String,
    /// Search results display text
    pub(crate) search_results: String,
    /// Available substances in selected library
    pub(crate) available_substances: Vec<String>,
    /// Currently selected substance for calculations
    pub(crate) selected_substance: String,
    /// Calculated thermodynamic properties
    pub(crate) calculated_cp: Option<f64>,
    pub(crate) calculated_dh: Option<f64>,
    pub(crate) calculated_ds: Option<f64>,
    pub(crate) energy_unit: EnergyUnit,
    pub(crate) show_plots_window: bool,
    pub(crate) t0: String,
    pub(crate) tend: String,
    pub(crate) plot_status: String,
    pub(crate) plot_window: Option<PlotWindow>,
    pub(crate) show_nist_window: bool,
    pub(crate) nist_app: Option<NISTApp>,
}

impl Default for ThermochemistryApp {
    fn default() -> Self {
        Self::from_catalog_result(ThermoData::try_new())
    }
}

impl ThermochemistryApp {
    pub(crate) fn from_catalog_result(
        catalog_result: Result<ThermoData, ThermoLibraryError>,
    ) -> Self {
        let (thermo_data, startup_error) = match catalog_result {
            Ok(thermo_data) => (thermo_data, None),
            Err(err) => (
                ThermoData::empty(),
                Some(format!("Thermo catalog failed to load: {}", err)),
            ),
        };

        let selected_library = thermo_data
            .thermo_libs
            .as_ref()
            .first()
            .cloned()
            .unwrap_or_else(|| "NASA_gas".to_string());

        let available_substances = thermo_data
            .LibThermoData
            .get(&selected_library)
            .map(|lib_data| lib_data.keys().cloned().collect())
            .unwrap_or_default();

        Self {
            thermo_data,
            startup_error,
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
            plot_status: "No temperature range calculated yet".to_string(),
            plot_window: None,
            show_nist_window: false,
            nist_app: None,
        }
    }

    /// Creates a new ThermochemistryApp instance
    ///
    /// Initializes the application with default values:
    /// - Empty substance search
    /// - First thermo library as default
    /// - Standard temperature snapshot (298.15 K)
    pub fn new() -> Self {
        Self::default()
    }

    fn startup_error_message(&self) -> Option<&str> {
        self.startup_error.as_deref()
    }

    fn clear_calculated_results(&mut self) {
        self.calculated_cp = None;
        self.calculated_dh = None;
        self.calculated_ds = None;
    }

    fn clear_plot_snapshot(&mut self) {
        self.plot_window = None;
    }

    pub(crate) fn invalidate_calculation_snapshot(&mut self) {
        self.clear_calculated_results();
        self.search_results = "No search performed yet".to_string();
        self.clear_plot_snapshot();
        self.show_plots_window = false;
        self.plot_status = "No temperature range calculated yet".to_string();
    }

    pub(crate) fn set_energy_unit(&mut self, unit: EnergyUnit) {
        if self.energy_unit != unit {
            self.energy_unit = unit;
            self.invalidate_calculation_snapshot();
        }
    }

    pub(crate) fn selected_thermochemistry_library_id(&self) -> Result<LibraryId, String> {
        LibraryId::resolve(&self.selected_library).ok_or_else(|| {
            format!(
                "Selected thermochemistry library '{}' is not recognized by the typed solver boundary",
                self.selected_library
            )
        })
    }

    pub(crate) fn build_calculation_subs_data(
        &self,
        substance_name: &str,
    ) -> Result<SubsData, String> {
        let library_id = self.selected_thermochemistry_library_id()?;
        let mut subs_data = SubsData::empty();
        subs_data.thermo_data = self.thermo_data.clone();
        subs_data.set_substances(vec![substance_name.to_string()]);
        subs_data.set_library_priority_id(library_id, LibraryPriority::Priority);
        subs_data
            .search_substances()
            .map_err(|e| format!("Failed to resolve substance routing: {}", e))?;
        Ok(subs_data)
    }

    /// Updates available substances when library changes
    pub(crate) fn update_substances_for_library(&mut self) {
        self.selected_substance.clear();
        self.invalidate_calculation_snapshot();
        self.available_substances = self
            .thermo_data
            .LibThermoData
            .get(&self.selected_library)
            .map(|lib_data| {
                let mut substances: Vec<String> = lib_data.keys().cloned().collect();
                substances.sort();
                substances
            })
            .unwrap_or_default();
    }

    /// Rejects transport-only libraries before we enter the thermochemistry factories.
    fn ensure_thermochemistry_backend(&self) -> Result<(), String> {
        if self
            .thermo_data
            .thermo_libs
            .as_ref()
            .iter()
            .any(|library| library == &self.selected_library)
        {
            Ok(())
        } else {
            Err(format!(
                "Library '{}' is not a supported thermochemistry backend",
                self.selected_library
            ))
        }
    }

    /// Updates the selected substance and shows the canonical raw entry snapshot.
    pub(crate) fn search_substance_by_name(&mut self, substance_name: &str) {
        if let Some(error) = self.startup_error_message() {
            self.search_results = error.to_string();
            return;
        }
        self.invalidate_calculation_snapshot();
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

    /// Calculates thermodynamic properties at specified conditions.
    ///
    /// This method integrates with the thermodynamics module to:
    /// 1. Parse the temperature input
    /// 2. Calculate Cp, enthalpy, and entropy
    /// 3. Display the canonical thermochemistry snapshot
    pub(crate) fn calculate_properties(&mut self) {
        if let Some(error) = self.startup_error_message() {
            self.search_results = error.to_string();
            self.clear_calculated_results();
            return;
        }

        self.clear_calculated_results();

        let temp_result = parse_temperature(&self.temperature);

        let substance_name = &self.selected_substance;

        if substance_name.is_empty() {
            self.search_results = "Please select a substance first".to_string();
            return;
        }

        if let Err(message) = self.ensure_thermochemistry_backend() {
            self.search_results = message;
            return;
        }

        match temp_result {
            Ok(t) => match self.build_calculation_subs_data(substance_name) {
                Ok(mut subs_data) => {
                    match self.perform_thermo_calculation(&mut subs_data, substance_name, t) {
                        Ok((cp, dh, ds)) => {
                            self.calculated_cp = Some(cp);
                            self.calculated_dh = Some(dh);
                            self.calculated_ds = Some(ds);
                            let (cp_unit, dh_unit, ds_unit) = match self.energy_unit {
                                EnergyUnit::J => ("J/(mol*K)", "kJ/mol", "J/(mol*K)"),
                                EnergyUnit::Cal => ("cal/(mol*K)", "kcal/mol", "cal/(mol*K)"),
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
                }
                Err(e) => {
                    self.search_results = format!("Calculation error: {}", e);
                }
            },
            Err(err) => {
                self.search_results = err.to_string();
            }
        }
    }

    /// Performs the actual thermodynamic calculation through the canonical `SubsData` boundary.
    fn perform_thermo_calculation(
        &self,
        subs_data: &mut SubsData,
        substance_name: &str,
        temperature: f64,
    ) -> Result<(f64, f64, f64), String> {
        self.ensure_thermochemistry_backend()?;
        subs_data
            .extract_thermal_coeffs(substance_name, temperature)
            .map_err(|e| format!("Failed to extract thermochemistry coefficients: {}", e))?;
        let (cp, dh, ds) = subs_data
            .calculate_thermo_properties(substance_name, temperature)
            .map_err(|e| {
                format!(
                    "Failed to calculate thermochemistry through SubsData: {}",
                    e
                )
            })?;

        Ok((cp, dh, ds))
    }

    /// Calculates thermodynamic properties over a temperature range
    fn calculate_temperature_range(
        &self,
        subs_data: &mut SubsData,
        substance_name: &str,
        t0: f64,
        tend: f64,
        n_points: usize,
    ) -> Result<Vec<(f64, f64, f64, f64)>, String> {
        let mut results = Vec::new();
        let dt = (tend - t0) / (n_points - 1) as f64;
        for i in 0..n_points {
            let t = t0 + i as f64 * dt;
            subs_data
                .extract_thermal_coeffs(substance_name, t)
                .map_err(|e| format!("Failed to extract thermochemistry coefficients: {}", e))?;
            let (cp, dh, ds) = subs_data
                .calculate_thermo_properties(substance_name, t)
                .map_err(|e| format!("Failed to calculate properties through SubsData: {}", e))?;

            results.push((t, cp, dh, ds));
        }

        Ok(results)
    }

    /// Shows the temperature-range controls window.
    fn render_temperature_range_window(&mut self, ctx: &egui::Context) {
        let mut window_open = self.show_plots_window;
        egui::Window::new("Temperature Range Plots")
            .open(&mut window_open)
            .default_size([800.0, 600.0])
            .show(ctx, |ui| {
                ui.heading("Temperature Range Analysis");
                ui.separator();
                ui.label(&self.plot_status);
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
        self.show_plots_window = window_open;
    }

    /// Calculates data for the temperature range and creates plot window
    pub(crate) fn calculate_range_data(&mut self) {
        if let Some(error) = self.startup_error_message() {
            let error_message = error.to_string();
            self.plot_status = error_message.clone();
            self.clear_plot_snapshot();
            warn!("{}", error_message);
            return;
        }

        let substance_name = &self.selected_substance;

        if substance_name.is_empty() {
            self.plot_status = "Please select a substance first".to_string();
            self.clear_plot_snapshot();
            warn!("{}", self.plot_status);
            return;
        }

        if self.ensure_thermochemistry_backend().is_err() {
            self.plot_status = format!(
                "Skipping temperature range calculation for unsupported thermochemistry backend '{}'",
                self.selected_library
            );
            self.clear_plot_snapshot();
            warn!("{}", self.plot_status);
            return;
        }

        match parse_temperature_range(&self.t0, &self.tend) {
            Ok(range) => {
                match self.build_calculation_subs_data(substance_name) {
                    Ok(mut subs_data) => match self.calculate_temperature_range(
                        &mut subs_data,
                        substance_name,
                        range.t0,
                        range.tend,
                        100,
                    ) {
                        Ok(results) => {
                            self.plot_status = format!(
                                "Temperature range ready for '{}' ({} points)",
                                substance_name,
                                results.len()
                            );
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
                                EnergyUnit::J => ("Cp (J/(mol*K))", "dH (J/mol)", "dS (J/(mol*K))"),
                                EnergyUnit::Cal => {
                                    ("Cp (cal/(mol*K))", "dH (cal/mol)", "dS (cal/(mol*K))")
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
                            self.plot_status = format!("Range calculation error: {}", e);
                            warn!("Range calculation error: {}", e);
                        }
                    },
                    Err(e) => {
                        self.plot_status = format!("Range preparation error: {}", e);
                        warn!("Range preparation error: {}", e);
                    }
                }
            }
            Err(err) => {
                self.plot_status = format!("Invalid temperature range: {}", err);
                self.clear_plot_snapshot();
                warn!("{}", self.plot_status);
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
    ///
    /// ### Middle Panel:
    /// - **Temperature Input**: Field for calculation temperature (K)
    /// - **Calculate Button**: Performs thermodynamic property calculations
    ///
    /// ### Bottom Panel:
    /// - **Results Display**: Scrollable text area showing search results and calculations
    ///
    /// ## Button Logic:
    /// - **Selection**: clicking a substance row updates the raw snapshot
    /// - **Calculate**: Calls `calculate_properties()` to compute thermodynamic properties
    pub fn show(&mut self, ctx: &egui::Context, open: &mut bool) {
        egui::Window::new("Thermochemistry Analysis")
            .open(open)
            .default_size([1200.0, 800.0])
            .show(ctx, |ui| {
                ui.heading("Thermodynamic Properties and Calculations");
                ui.separator();

                if let Some(error) = self.startup_error_message() {
                    ui.colored_label(egui::Color32::from_rgb(180, 40, 40), error);
                    ui.separator();
                }

                ui.add_enabled_ui(self.startup_error.is_none(), |ui| {
                    ui.horizontal(|ui| {
                        ui.vertical(|ui| {
                            ui.set_width(200.0);
                            ui.set_min_height(600.0);
                            ui.heading("List of Substances");
                            ui.horizontal(|ui| {
                                ui.label("Search:");
                                ui.text_edit_singleline(&mut self.substance_input);
                            });
                            ui.separator();
                            egui::ScrollArea::vertical()
                                .max_height(1200.0)
                                .show(ui, |ui| {
                                    let available_substances = self.available_substances.clone();
                                    for substance in &available_substances {
                                        if self.substance_input.is_empty()
                                            || substance
                                                .to_lowercase()
                                                .contains(&self.substance_input.to_lowercase())
                                        {
                                            if ui
                                                .selectable_label(
                                                    self.selected_substance == *substance,
                                                    substance,
                                                )
                                                .clicked()
                                            {
                                                self.search_substance_by_name(substance);
                                            }
                                        }
                                    }
                                });
                            ui.separator();
                            egui::ComboBox::from_label("Library Source")
                                .selected_text(&self.selected_library)
                                .show_ui(ui, |ui| {
                                    let thermolibs = self.thermo_data.thermo_libs.as_ref().clone();
                                    for library in &thermolibs {
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

                        ui.vertical(|ui| {
                            ui.group(|ui| {
                                ui.label("Substance Information (read-only):");
                                egui::ScrollArea::vertical()
                                    .id_salt("substance_info")
                                    .max_height(300.0)
                                    .show(ui, |ui| {
                                        read_only_multiline(ui, &mut self.search_results, 12);
                                    });
                            });

                            ui.add_space(10.0);

                            ui.group(|ui| {
                                ui.label("Calculation Parameters:");
                                ui.horizontal(|ui| {
                                    ui.label("Temperature (K):");
                                    if ui.text_edit_singleline(&mut self.temperature).changed() {
                                        self.invalidate_calculation_snapshot();
                                    }

                                    if ui.button("Calculate Properties").clicked() {
                                        self.calculate_properties();
                                    }
                                });
                                ui.horizontal(|ui| {
                                    ui.label("Energy Unit:");
                                    if ui
                                        .radio_value(&mut self.energy_unit, EnergyUnit::J, "J")
                                        .changed()
                                    {
                                        self.set_energy_unit(EnergyUnit::J);
                                    }
                                    if ui
                                        .radio_value(&mut self.energy_unit, EnergyUnit::Cal, "cal")
                                        .changed()
                                    {
                                        self.set_energy_unit(EnergyUnit::Cal);
                                    }
                                });
                            });

                            ui.separator();

                            ui.horizontal(|ui| {
                                if ui.button("Clear Results").clicked() {
                                    self.invalidate_calculation_snapshot();
                                }

                                if ui.button("View Plots").clicked() {
                                    self.show_plots_window = true;
                                }

                                if ui.button("Export Data").clicked() {
                                    self.search_results +=
                                        "\n\nExport functionality not yet implemented";
                                }

                                if ui.button("Search in NIST").clicked() {
                                    self.show_nist_window = true;
                                    if self.nist_app.is_none() {
                                        self.nist_app = Some(NISTApp::new());
                                    }
                                }
                            });
                        });
                    });
                });

                // Show plots window if requested
                if self.show_plots_window {
                    self.render_temperature_range_window(ctx);
                }

                // Show plot window if available
                if let Some(plot_window) = &mut self.plot_window {
                    plot_window.show(ctx);
                    if !plot_window.visible {
                        self.plot_window = None;
                    }
                }

                // Show NIST window if requested
                if self.show_nist_window {
                    if let Some(nist_app) = &mut self.nist_app {
                        nist_app.show(ctx, &mut self.show_nist_window);
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
    fn test_calculate_properties_invalid_temperature() {
        let mut app = ThermochemistryApp::new();
        app.temperature = "invalid".to_string();
        app.calculate_properties();
        assert_eq!(app.search_results, "Please select a substance first");
    }
}
