//! # Transport Properties GUI Module
//!
//! This module provides a graphical user interface for transport property calculations.
//! It allows users to search for transport data and calculate viscosity and thermal conductivity
//! from various sources (CEA, Aramco transport).
//!
//! ## Main Features
//! - **Substance Selection**: Inspect transport data across multiple libraries
//! - **Property Calculation**: Calculate viscosity and thermal conductivity
//! - **Library Management**: Browse and select from different transport databases
//! - **Unit Management**: Support for different unit systems
//!
//! ## GUI Layout
//! The interface consists of:
//! - **Selection Panel**: Substance input and library selection
//! - **Results Panel**: Display of found transport data
//! - **Calculation Panel**: Temperature input and transport property calculations

use eframe::egui;

use crate::Thermodynamics::DBhandlers::transport_api::{LambdaUnit, ViscosityUnit};
use crate::Thermodynamics::User_substances::{LibraryPriority, SubsData};
use crate::Thermodynamics::thermo_lib_api::{LibraryId, ThermoData, ThermoLibraryError};
use crate::gui::condition_parser::{parse_pressure, parse_temperature};
use crate::gui::read_only_snapshot::multiline as read_only_multiline;
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
    pub(crate) transport_data: ThermoData,
    /// Visible startup error if the catalog could not be loaded.
    pub(crate) startup_error: Option<String>,
    /// Input field for substance names to search
    pub(crate) substance_input: String,
    /// Currently selected transport library
    pub(crate) selected_library: String,
    /// Temperature for calculations (in Kelvin)
    pub(crate) temperature: String,
    /// Pressure for calculations (in Pa)
    pub(crate) pressure: String,
    /// Search results display text
    pub(crate) search_results: String,
    /// Available substances in selected library
    pub(crate) available_substances: Vec<String>,
    /// Currently selected substance for calculations
    pub(crate) selected_substance: String,
    /// Calculated transport properties
    pub(crate) calculated_lambda: Option<f64>,
    pub(crate) calculated_viscosity: Option<f64>,
    /// Temperature used for the latest successful transport report.
    pub(crate) last_calculated_temperature: Option<f64>,
    pub(crate) lambda_unit: LambdaUnit,
    pub(crate) viscosity_unit: ViscosityUnit,
}

impl Default for TransportApp {
    fn default() -> Self {
        Self::from_catalog_result(ThermoData::try_new())
    }
}

impl TransportApp {
    pub(crate) fn from_catalog_result(
        catalog_result: Result<ThermoData, ThermoLibraryError>,
    ) -> Self {
        let (transport_data, startup_error) = match catalog_result {
            Ok(transport_data) => (transport_data, None),
            Err(err) => (
                ThermoData::empty(),
                Some(format!("Thermo catalog failed to load: {}", err)),
            ),
        };

        let selected_library = transport_data
            .transport_libs
            .as_ref()
            .first()
            .cloned()
            .unwrap_or_else(|| "Aramco_transport".to_string());

        let available_substances = transport_data
            .LibThermoData
            .get(&selected_library)
            .map(|lib_data| lib_data.keys().cloned().collect())
            .unwrap_or_default();

        Self {
            transport_data,
            startup_error,
            substance_input: String::new(),
            selected_library,
            temperature: "298.15".to_string(),
            pressure: "101325.0".to_string(),
            search_results: "No search performed yet".to_string(),
            available_substances,
            selected_substance: String::new(),
            calculated_lambda: None,
            calculated_viscosity: None,
            last_calculated_temperature: None,
            lambda_unit: LambdaUnit::MKWPerMK,
            viscosity_unit: ViscosityUnit::MKPaS,
        }
    }

    /// Creates a new TransportApp instance
    ///
    /// Initializes the application with default values:
    /// - Empty substance search
    /// - First transport library as default
    /// - Standard temperature and pressure (298.15 K, 101325 Pa)
    pub fn new() -> Self {
        Self::default()
    }

    fn startup_error_message(&self) -> Option<&str> {
        self.startup_error.as_deref()
    }

    fn clear_calculated_results(&mut self) {
        self.calculated_lambda = None;
        self.calculated_viscosity = None;
        self.last_calculated_temperature = None;
    }

    pub(crate) fn invalidate_calculation_snapshot(&mut self) {
        self.clear_calculated_results();
        self.search_results = "No search performed yet".to_string();
    }

    pub(crate) fn set_lambda_unit(&mut self, unit: LambdaUnit) {
        if self.lambda_unit != unit {
            self.lambda_unit = unit;
            if self.has_calculation_snapshot() {
                self.refresh_calculation_snapshot();
            }
        }
    }

    pub(crate) fn set_viscosity_unit(&mut self, unit: ViscosityUnit) {
        if self.viscosity_unit != unit {
            self.viscosity_unit = unit;
            if self.has_calculation_snapshot() {
                self.refresh_calculation_snapshot();
            }
        }
    }

    fn has_calculation_snapshot(&self) -> bool {
        self.calculated_lambda.is_some()
            && self.calculated_viscosity.is_some()
            && self.last_calculated_temperature.is_some()
            && !self.selected_substance.is_empty()
    }

    fn lambda_display_multiplier(&self) -> f64 {
        match self.lambda_unit {
            LambdaUnit::WPerMK => 1.0,
            LambdaUnit::MWPerMK => 1.0e3,
            LambdaUnit::MKWPerMK => 1.0e6,
            LambdaUnit::MKWPerSMK => 1.0e4,
        }
    }

    fn viscosity_display_multiplier(&self) -> f64 {
        match self.viscosity_unit {
            ViscosityUnit::KgPerMS => 1.0,
            ViscosityUnit::PaS => 1.0,
            ViscosityUnit::MKPaS => 1.0e6,
        }
    }

    fn lambda_unit_label(&self) -> &'static str {
        match self.lambda_unit {
            LambdaUnit::WPerMK => "W/m*K",
            LambdaUnit::MWPerMK => "mW/m*K",
            LambdaUnit::MKWPerMK => "uW/m*K",
            LambdaUnit::MKWPerSMK => "uW/s*m*K",
        }
    }

    fn viscosity_unit_label(&self) -> &'static str {
        match self.viscosity_unit {
            ViscosityUnit::KgPerMS => "kg/m*s",
            ViscosityUnit::PaS => "Pa*s",
            ViscosityUnit::MKPaS => "uPa*s",
        }
    }

    fn format_calculation_report(
        &self,
        substance_name: &str,
        temperature: f64,
        lambda: f64,
        viscosity: f64,
    ) -> String {
        let lambda_display = lambda * self.lambda_display_multiplier();
        let viscosity_display = viscosity * self.viscosity_display_multiplier();

        format!(
            "Transport properties for '{}' at T = {} K:\n\nThermal conductivity = {:.5} {}\nViscosity = {:.8} {}",
            substance_name,
            temperature,
            lambda_display,
            self.lambda_unit_label(),
            viscosity_display,
            self.viscosity_unit_label()
        )
    }

    fn refresh_calculation_snapshot(&mut self) {
        if let (Some(lambda), Some(viscosity), Some(temperature)) = (
            self.calculated_lambda,
            self.calculated_viscosity,
            self.last_calculated_temperature,
        ) {
            self.search_results = self.format_calculation_report(
                &self.selected_substance,
                temperature,
                lambda,
                viscosity,
            );
        }
    }

    pub(crate) fn selected_transport_library_id(&self) -> Result<LibraryId, String> {
        LibraryId::resolve(&self.selected_library).ok_or_else(|| {
            format!(
                "Selected transport library '{}' is not recognized by the typed solver boundary",
                self.selected_library
            )
        })
    }

    pub(crate) fn build_calculation_subs_data(
        &self,
        substance_name: &str,
        pressure: f64,
    ) -> Result<SubsData, String> {
        let transport_library = self.selected_transport_library_id()?;
        let mut subs_data = SubsData::empty();
        subs_data.thermo_data = self.transport_data.clone();
        subs_data.set_substances(vec![substance_name.to_string()]);
        subs_data.set_library_priority_id(transport_library, LibraryPriority::Priority);

        for library in self.transport_data.thermo_libs.as_ref().iter() {
            if let Some(library_id) = LibraryId::resolve(library) {
                subs_data.set_library_priority_id(library_id, LibraryPriority::Permitted);
            }
        }

        subs_data
            .search_substances()
            .map_err(|e| format!("Failed to resolve substance routing: {}", e))?;
        subs_data
            .calculate_elem_composition_and_molar_mass(None)
            .map_err(|e| format!("Failed to derive molar mass and composition: {}", e))?;

        if transport_library != LibraryId::Cea {
            let molar_masses = subs_data.molar_mass_by_substance.clone();
            subs_data
                .set_M(molar_masses, Some("g/mol".to_string()))
                .map_err(|e| format!("Failed to set molar mass metadata: {}", e))?;
            subs_data
                .set_P(pressure, Some("Pa".to_string()))
                .map_err(|e| format!("Failed to set pressure: {}", e))?;
        }

        Ok(subs_data)
    }

    /// Updates available substances when library changes
    pub(crate) fn update_substances_for_library(&mut self) {
        self.selected_substance.clear();
        self.invalidate_calculation_snapshot();
        self.available_substances = self
            .transport_data
            .LibThermoData
            .get(&self.selected_library)
            .map(|lib_data| {
                let mut substances: Vec<String> = lib_data.keys().cloned().collect();
                substances.sort();
                substances
            })
            .unwrap_or_default();
    }

    /// Updates the selected substance and shows the canonical raw entry snapshot.
    pub(crate) fn search_substance_by_name(&mut self, substance_name: &str) {
        if let Some(error) = self.startup_error_message() {
            self.search_results = error.to_string();
            return;
        }
        self.invalidate_calculation_snapshot();
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

    /// Calculates transport properties at specified conditions
    ///
    /// This method integrates with the transport module to:
    /// 1. Parse temperature and pressure inputs
    /// 2. Calculate viscosity and thermal conductivity
    /// 3. Display calculated values
    pub(crate) fn calculate_properties(&mut self) {
        if let Some(error) = self.startup_error_message() {
            self.search_results = error.to_string();
            self.clear_calculated_results();
            return;
        }

        self.clear_calculated_results();

        let temp_result = parse_temperature(&self.temperature);
        let pres_result = parse_pressure(&self.pressure);

        let substance_name = &self.selected_substance;

        if substance_name.is_empty() {
            self.search_results = "Please select a substance first".to_string();
            return;
        }

        match (temp_result, pres_result) {
            (Ok(t), Ok(p)) => match self.build_calculation_subs_data(substance_name, p) {
                Ok(mut subs_data) => {
                    match self.perform_transport_calculation(&mut subs_data, substance_name, t) {
                        Ok((lambda, viscosity)) => {
                            self.calculated_lambda = Some(lambda);
                            self.calculated_viscosity = Some(viscosity);
                            self.last_calculated_temperature = Some(t);
                            self.search_results = self.format_calculation_report(
                                substance_name,
                                t,
                                lambda,
                                viscosity,
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
            (Err(err), _) | (_, Err(err)) => {
                self.search_results = err.to_string();
            }
        }
    }

    /// Performs the actual transport calculation through the canonical `SubsData` boundary.
    fn perform_transport_calculation(
        &self,
        subs_data: &mut SubsData,
        substance_name: &str,
        temperature: f64,
    ) -> Result<(f64, f64), String> {
        subs_data
            .extract_thermal_coeffs(substance_name, temperature)
            .map_err(|e| format!("Failed to extract thermochemistry coefficients: {}", e))?;
        let (cp, _, _) = subs_data
            .calculate_thermo_properties(substance_name, temperature)
            .map_err(|e| format!("Failed to calculate Cp through SubsData: {}", e))?;
        subs_data
            .extract_transport_coeffs(substance_name, temperature)
            .map_err(|e| format!("Failed to extract transport coefficients: {}", e))?;
        let (lambda, viscosity) = subs_data
            .calculate_transport_properties(substance_name, temperature, Some(cp), None)
            .map_err(|e| {
                format!(
                    "Failed to calculate transport properties through SubsData: {}",
                    e
                )
            })?;

        Ok((lambda, viscosity))
    }

    /// Main GUI rendering method for the transport properties interface
    ///
    /// ## Layout Structure:
    ///
    /// ### Top Panel:
    /// - **Substance Input**: Text field for entering substance names
    /// - **Library Selection**: Dropdown for choosing transport database
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
    /// - **Selection**: clicking a substance row updates the raw snapshot
    /// - **Calculate**: Calls `calculate_properties()` to compute transport properties
    pub fn show(&mut self, ctx: &egui::Context, open: &mut bool) {
        egui::Window::new("Transport Properties Analysis")
            .open(open)
            .default_size([1200.0, 800.0])
            .show(ctx, |ui| {
                ui.heading("Transport Properties and Calculations");
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
                                    let transport_libs =
                                        self.transport_data.transport_libs.as_ref().clone();
                                    for library in &transport_libs {
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

                                    ui.label("Pressure (Pa):");
                                    if ui.text_edit_singleline(&mut self.pressure).changed() {
                                        self.invalidate_calculation_snapshot();
                                    }

                                    if ui.button("Calculate Properties").clicked() {
                                        self.calculate_properties();
                                    }
                                });
                                ui.horizontal(|ui| {
                                    ui.label("Lambda Unit:");
                                    if ui
                                        .radio_value(
                                            &mut self.lambda_unit,
                                            LambdaUnit::WPerMK,
                                            "W/m*K",
                                        )
                                        .changed()
                                    {
                                        self.set_lambda_unit(LambdaUnit::WPerMK);
                                    }
                                    if ui
                                        .radio_value(
                                            &mut self.lambda_unit,
                                            LambdaUnit::MWPerMK,
                                            "mW/m*K",
                                        )
                                        .changed()
                                    {
                                        self.set_lambda_unit(LambdaUnit::MWPerMK);
                                    }
                                    if ui
                                        .radio_value(
                                            &mut self.lambda_unit,
                                            LambdaUnit::MKWPerMK,
                                            "uW/m*K",
                                        )
                                        .changed()
                                    {
                                        self.set_lambda_unit(LambdaUnit::MKWPerMK);
                                    }
                                });
                                ui.horizontal(|ui| {
                                    ui.label("Viscosity Unit:");
                                    if ui
                                        .radio_value(
                                            &mut self.viscosity_unit,
                                            ViscosityUnit::PaS,
                                            "Pa*s",
                                        )
                                        .changed()
                                    {
                                        self.set_viscosity_unit(ViscosityUnit::PaS);
                                    }
                                    if ui
                                        .radio_value(
                                            &mut self.viscosity_unit,
                                            ViscosityUnit::MKPaS,
                                            "uPa*s",
                                        )
                                        .changed()
                                    {
                                        self.set_viscosity_unit(ViscosityUnit::MKPaS);
                                    }
                                });
                            });

                            ui.separator();

                            ui.horizontal(|ui| {
                                if ui.button("Clear Results").clicked() {
                                    self.invalidate_calculation_snapshot();
                                }

                                if ui.button("Export Data").clicked() {
                                    self.search_results +=
                                        "\\n\\nExport functionality not yet implemented";
                                }

                                if ui.button("Load from File").clicked() {
                                    self.search_results +=
                                        "\\n\\nFile loading functionality not yet implemented";
                                }
                            });
                        });
                    });
                });
            });
    }
}
