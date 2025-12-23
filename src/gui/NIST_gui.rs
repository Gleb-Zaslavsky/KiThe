use crate::Thermodynamics::DBhandlers::NIST_parser::{Phase, SearchType};
use crate::Thermodynamics::DBhandlers::NISTdata::NISTdata;
use crate::gui::gui_plot::PlotWindow;
use eframe::egui;
use nalgebra::{DMatrix, DVector};

#[derive(Debug)]
pub struct NISTApp {
    substance_input: String,
    search_type: SearchType,
    phase: Phase,
    search_results: String,
    temperature: String,
    unit: String,
    nist_data: NISTdata,
    show_plots_window: bool,
    t0: String,
    tend: String,
    plot_window: Option<PlotWindow>,
}

impl Default for NISTApp {
    fn default() -> Self {
        Self {
            substance_input: String::new(),
            search_type: SearchType::All,
            phase: Phase::Gas,
            search_results: "No search performed yet".to_string(),
            temperature: "298.15".to_string(),
            unit: "J".to_string(),
            nist_data: NISTdata::new(),
            show_plots_window: false,
            t0: "298.15".to_string(),
            tend: "1000.0".to_string(),
            plot_window: None,
        }
    }
}

impl NISTApp {
    pub fn new() -> Self {
        Self::default()
    }

    fn search_nist_data(&mut self) {
        if self.substance_input.is_empty() {
            self.search_results = "Please enter a substance name".to_string();
            return;
        }

        match self.nist_data.get_data_from_NIST(
            self.substance_input.clone(),
            self.search_type.clone(),
            self.phase.clone(),
        ) {
            Ok(_) => {
                if let Err(e) = self.nist_data.set_unit(&self.unit) {
                    self.search_results = format!("Unit error: {}", e);
                    return;
                }

                if let Ok(temp) = self.temperature.parse::<f64>() {
                    if let Err(e) = self.nist_data.extract_coefficients(temp) {
                        self.search_results = format!("Coefficient error: {}", e);
                        return;
                    }

                    if let Err(e) = self.nist_data.calculate_cp_dh_ds(temp) {
                        self.search_results = format!("Calculation error: {}", e);
                        return;
                    }

                    let unit_suffix = if self.unit == "J" { "J" } else { "cal" };
                    self.search_results = format!(
                        "NIST data for '{}' ({:?}, {:?}):\n\nAt T = {} K:\nCp = {:.3} {}/mol·K\ndH = {:.3} k{}/mol\ndS = {:.3} {}/mol·K",
                        self.substance_input,
                        self.phase,
                        self.search_type,
                        temp,
                        self.nist_data.Cp,
                        unit_suffix,
                        self.nist_data.dh / 1000.0,
                        unit_suffix,
                        self.nist_data.ds,
                        unit_suffix
                    );
                } else {
                    self.search_results = format!("Found NIST data for '{}'", self.substance_input);
                }
            }
            Err(e) => {
                self.search_results = format!("NIST search error: {}", e);
            }
        }
    }

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
            });
    }

    fn calculate_range_data(&mut self) {
        if self.substance_input.is_empty() {
            return;
        }

        let t0_result = self.t0.parse::<f64>();
        let tend_result = self.tend.parse::<f64>();

        match (t0_result, tend_result) {
            (Ok(t0), Ok(tend)) if t0 < tend => {
                match self.calculate_temperature_range(t0, tend, 100) {
                    Ok(results) => {
                        let temperatures: Vec<f64> =
                            results.iter().map(|(t, _, _, _)| *t).collect();
                        let t_result = DVector::from_vec(temperatures);

                        let mut y_data = Vec::new();
                        for (_, cp, dh, ds) in &results {
                            y_data.push(*cp);
                            y_data.push(*dh);
                            y_data.push(*ds);
                        }

                        let y_result = DMatrix::from_row_slice(results.len(), 3, &y_data);

                        let unit_suffix = if self.unit == "J" { "J" } else { "cal" };
                        let values = vec![
                            format!("Cp ({}/mol·K)", unit_suffix),
                            format!("dH ({}/mol)", unit_suffix),
                            format!("dS ({}/mol·K)", unit_suffix),
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
            _ => {
                println!("Invalid temperature range");
            }
        }
    }

    fn calculate_temperature_range(
        &mut self,
        t0: f64,
        tend: f64,
        n_points: usize,
    ) -> Result<Vec<(f64, f64, f64, f64)>, String> {
        let mut results = Vec::new();
        let dt = (tend - t0) / (n_points - 1) as f64;

        for i in 0..n_points {
            let t = t0 + i as f64 * dt;

            if let Err(e) = self.nist_data.extract_coefficients(t) {
                continue; // Skip temperatures outside coefficient ranges
            }

            if let Err(e) = self.nist_data.calculate_cp_dh_ds(t) {
                continue; // Skip calculation errors
            }

            results.push((t, self.nist_data.Cp, self.nist_data.dh, self.nist_data.ds));
        }

        if results.is_empty() {
            return Err("No valid temperature points found in range".to_string());
        }

        Ok(results)
    }

    pub fn show(&mut self, ctx: &egui::Context, open: &mut bool) {
        egui::Window::new("NIST Database Search")
            .open(open)
            .default_size([600.0, 500.0])
            .show(ctx, |ui| {
                ui.heading("NIST Chemistry WebBook Search");
                ui.separator();

                ui.horizontal(|ui| {
                    ui.label("Substance:");
                    ui.text_edit_singleline(&mut self.substance_input);
                });

                ui.horizontal(|ui| {
                    ui.label("Phase:");
                    ui.radio_value(&mut self.phase, Phase::Gas, "Gas");
                    ui.radio_value(&mut self.phase, Phase::Liquid, "Liquid");
                    ui.radio_value(&mut self.phase, Phase::Solid, "Solid");
                });

                ui.horizontal(|ui| {
                    ui.label("Search Type:");
                    ui.radio_value(&mut self.search_type, SearchType::All, "All");
                    ui.radio_value(&mut self.search_type, SearchType::Cp, "Cp");
                    ui.radio_value(&mut self.search_type, SearchType::DeltaH, "ΔH");
                    ui.radio_value(&mut self.search_type, SearchType::DeltaS, "ΔS");
                    ui.radio_value(&mut self.search_type, SearchType::MolarMass, "Molar Mass");
                });

                ui.horizontal(|ui| {
                    ui.label("Temperature (K):");
                    ui.text_edit_singleline(&mut self.temperature);
                });

                ui.horizontal(|ui| {
                    ui.label("Units:");
                    ui.radio_value(&mut self.unit, "J".to_string(), "Joules");
                    ui.radio_value(&mut self.unit, "cal".to_string(), "Calories");
                });

                if ui.button("Search!").clicked() {
                    self.search_nist_data();
                }

                if ui.button("View Plots").clicked() {
                    self.show_plots_window = true;
                }

                if ui.button("Save to library").clicked() {
                    // TODO: Implement save to library functionality
                }

                ui.separator();

                ui.group(|ui| {
                    ui.label("Search Results:");
                    egui::ScrollArea::vertical()
                        .max_height(300.0)
                        .show(ui, |ui| {
                            ui.text_edit_multiline(&mut self.search_results);
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
    fn test_nist_app_new() {
        let app = NISTApp::new();
        assert!(app.substance_input.is_empty());
        assert_eq!(app.search_type, SearchType::All);
        assert_eq!(app.phase, Phase::Gas);
        assert_eq!(app.temperature, "298.15");
        assert_eq!(app.unit, "J");
        assert_eq!(app.t0, "298.15");
        assert_eq!(app.tend, "1000.0");
        assert!(!app.show_plots_window);
        assert!(app.plot_window.is_none());
    }

    #[test]
    fn test_search_empty_substance() {
        let mut app = NISTApp::new();
        app.search_nist_data();
        assert_eq!(app.search_results, "Please enter a substance name");
    }

    #[test]
    fn test_temperature_range_validation() {
        let mut app = NISTApp::new();
        app.substance_input = "CH4".to_string();
        app.t0 = "1000.0".to_string();
        app.tend = "500.0".to_string();

        app.calculate_range_data();
        assert!(app.plot_window.is_none());
    }

    #[test]
    fn test_unit_selection() {
        let mut app = NISTApp::new();
        assert_eq!(app.unit, "J");

        app.unit = "cal".to_string();
        assert_eq!(app.unit, "cal");
    }

    #[test]
    fn test_search_type_variants() {
        let mut app = NISTApp::new();

        app.search_type = SearchType::Cp;
        assert_eq!(app.search_type, SearchType::Cp);

        app.search_type = SearchType::DeltaH;
        assert_eq!(app.search_type, SearchType::DeltaH);
    }

    #[test]
    fn test_phase_variants() {
        let mut app = NISTApp::new();

        app.phase = Phase::Liquid;
        assert_eq!(app.phase, Phase::Liquid);

        app.phase = Phase::Solid;
        assert_eq!(app.phase, Phase::Solid);
    }

    #[test]
    fn test_temperature_parsing() {
        let app = NISTApp::new();
        assert!(app.temperature.parse::<f64>().is_ok());

        let mut app2 = NISTApp::new();
        app2.temperature = "invalid".to_string();
        assert!(app2.temperature.parse::<f64>().is_err());
    }
}
