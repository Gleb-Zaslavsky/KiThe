use eframe::egui;
use egui_plot::{Line, Plot, PlotPoints};
use log::{info, warn};
use nalgebra::{DMatrix, DVector};
use plotters::prelude::*;
use std::path::Path;
use std::path::PathBuf;

#[derive(Debug, Clone, PartialEq, Eq)]
enum PlotExportStatus {
    Saved {
        plot_name: String,
        path: PathBuf,
    },
    Failed {
        plot_name: String,
        path: PathBuf,
        message: String,
    },
}

impl PlotExportStatus {
    fn matches_plot(&self, plot_name: &str) -> bool {
        match self {
            Self::Saved {
                plot_name: stored, ..
            }
            | Self::Failed {
                plot_name: stored, ..
            } => stored == plot_name,
        }
    }

    fn status_text(&self) -> String {
        match self {
            Self::Saved { path, .. } => format!("Saved plot to {}.", path.display()),
            Self::Failed { path, message, .. } => {
                format!("Plot export failed for {}: {}", path.display(), message)
            }
        }
    }

    fn is_error(&self) -> bool {
        matches!(self, Self::Failed { .. })
    }
}
#[derive(Debug)]
pub struct PlotWindow {
    pub visible: bool,
    pub arg: String,
    pub values: Vec<String>,
    pub t_result: DVector<f64>,
    pub y_result: DMatrix<f64>,
    pub output_dir: PathBuf,
    last_export_status: Option<PlotExportStatus>,
}

impl PlotWindow {
    pub fn new(
        arg: String,
        values: Vec<String>,
        t_result: DVector<f64>,
        y_result: DMatrix<f64>,
    ) -> Self {
        Self {
            visible: true,
            arg,
            values,
            t_result,
            y_result,
            output_dir: PathBuf::from("plots"),
            last_export_status: None,
        }
    }

    /// Validate that the mesh and solution were handed off without losing rows
    /// or columns.
    fn validate_dimensions(&self) -> Result<(), String> {
        if self.t_result.len() != self.y_result.nrows() {
            return Err(format!(
                "Plot data mismatch: mesh has {} points but solution has {} rows",
                self.t_result.len(),
                self.y_result.nrows()
            ));
        }

        if self.values.len() != self.y_result.ncols() {
            return Err(format!(
                "Plot data mismatch: {} variable labels but solution has {} columns",
                self.values.len(),
                self.y_result.ncols()
            ));
        }

        Ok(())
    }

    pub fn show(&mut self, ctx: &egui::Context) {
        if !self.visible {
            return;
        }

        let validation_error = self.validate_dimensions().err();
        let mut pending_export: Option<(String, Vec<f64>, Vec<f64>)> = None;
        egui::Window::new("Simulation Results")
            .open(&mut self.visible)
            .resizable(true)
            .vscroll(true)
            .show(ctx, |ui| {
                ui.heading(format!("Profiles of {}", self.arg));
                ui.separator();

                if let Some(error) = &validation_error {
                    ui.colored_label(egui::Color32::from_rgb(180, 30, 30), error);
                    return;
                }

                let x_data: Vec<f64> = self.t_result.iter().cloned().collect();

                for col in 0..self.values.len() {
                    if col < self.y_result.ncols() {
                        let plot_name = self.values[col].clone();
                        let y_values = self.y_result.column(col).as_slice().to_vec();
                        ui.label(format!("Plot: {}", plot_name));

                        Plot::new(format!("plot_{}", col))
                            .height(200.0)
                            .show(ui, |plot_ui| {
                                let points: PlotPoints = x_data
                                    .iter()
                                    .zip(y_values.iter())
                                    .map(|(&x, &y)| [x, y])
                                    .collect();

                                let line =
                                    Line::new(plot_name.clone(), points).name(plot_name.clone());
                                plot_ui.line(line);
                            });

                        ui.horizontal(|ui| {
                            ui.add_space(20.0);
                            if ui.button("💾 Save plot").clicked() {
                                pending_export =
                                    Some((plot_name.clone(), x_data.clone(), y_values.clone()));
                            }
                        });

                        if let Some(status) = &self.last_export_status {
                            if status.matches_plot(&plot_name) {
                                let status_color = if status.is_error() {
                                    egui::Color32::from_rgb(180, 30, 30)
                                } else {
                                    egui::Color32::from_rgb(20, 120, 40)
                                };
                                ui.colored_label(status_color, status.status_text());
                            }
                        }

                        ui.add_space(15.0);
                        ui.separator();
                    }
                }
            });

        if let Some((plot_name, x_values, y_values)) = pending_export {
            if let Err(err) = self.export_plot_series(&plot_name, &x_values, &y_values) {
                warn!("Failed to save plot {}: {}", plot_name, err);
            }
        }
    }

    fn export_plot_series(
        &mut self,
        var_name: &str,
        x: &[f64],
        y: &[f64],
    ) -> Result<PathBuf, Box<dyn std::error::Error>> {
        let filename = format!("{}.png", var_name.replace('/', "_"));
        let output_path = self.output_dir.join(filename);

        match Self::save_plot_to_png(&output_path, var_name, &self.arg, x, y) {
            Ok(saved_path) => {
                self.last_export_status = Some(PlotExportStatus::Saved {
                    plot_name: var_name.to_string(),
                    path: saved_path.clone(),
                });
                Ok(saved_path)
            }
            Err(err) => {
                self.last_export_status = Some(PlotExportStatus::Failed {
                    plot_name: var_name.to_string(),
                    path: output_path,
                    message: err.to_string(),
                });
                Err(err)
            }
        }
    }

    fn save_plot_to_png(
        output_path: &Path,
        var_name: &str,
        arg: &str,
        x: &[f64],
        y: &[f64],
    ) -> Result<PathBuf, Box<dyn std::error::Error>> {
        if x.len() != y.len() {
            return Err(format!(
                "plot data mismatch: mesh has {} points but series has {} values",
                x.len(),
                y.len()
            )
            .into());
        }

        if let Some(parent) = output_path.parent() {
            std::fs::create_dir_all(parent)?;
        }

        let root = BitMapBackend::new(output_path, (800, 600)).into_drawing_area();
        root.fill(&WHITE)?;

        let x_min = *x.first().unwrap_or(&0.0);
        let x_max = *x.last().unwrap_or(&1.0);
        let (y_min, y_max) = y
            .iter()
            .fold((f64::INFINITY, f64::NEG_INFINITY), |(min, max), &val| {
                (min.min(val), max.max(val))
            });

        let mut chart = ChartBuilder::on(&root)
            .caption(format!("{} vs {}", var_name, arg), ("sans-serif", 24))
            .margin(15)
            .x_label_area_size(40)
            .y_label_area_size(60)
            .build_cartesian_2d(x_min..x_max, y_min..y_max)?;

        chart.configure_mesh().x_desc(arg).y_desc(var_name).draw()?;

        chart.draw_series(LineSeries::new(
            x.iter().zip(y.iter()).map(|(&xv, &yv)| (xv, yv)),
            &RED,
        ))?;

        chart
            .configure_series_labels()
            .border_style(&BLACK)
            .draw()?;
        root.present()?;

        info!("Saved plot to {:?}", output_path);
        Ok(output_path.to_path_buf())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::{DMatrix, DVector};

    #[test]
    fn test_plot_window_creation() {
        let arg = "time".to_string();
        let values = vec!["temperature".to_string(), "pressure".to_string()];
        let t_result = DVector::from_vec(vec![0.0, 1.0, 2.0, 3.0]);
        let y_result = DMatrix::from_row_slice(
            4,
            2,
            &[
                300.0, 101325.0, 350.0, 102000.0, 400.0, 103000.0, 450.0, 104000.0,
            ],
        );

        let plot_window = PlotWindow::new(
            arg.clone(),
            values.clone(),
            t_result.clone(),
            y_result.clone(),
        );

        assert_eq!(plot_window.arg, arg);
        assert_eq!(plot_window.values, values);
        assert_eq!(plot_window.t_result, t_result);
        assert_eq!(plot_window.y_result, y_result);
        assert!(plot_window.visible);
    }

    #[test]
    fn test_plot_window_with_chemical_data() {
        let arg = "position".to_string();
        let values = vec!["CO".to_string(), "CO2".to_string(), "H2O".to_string()];
        let t_result = DVector::from_vec(vec![0.0, 0.25, 0.5, 0.75, 1.0]);
        let y_result = DMatrix::from_row_slice(
            5,
            3,
            &[
                0.1, 0.0, 0.0, // initial concentrations
                0.08, 0.01, 0.01, 0.06, 0.02, 0.02, 0.04, 0.03, 0.03, 0.02, 0.04,
                0.04, // final concentrations
            ],
        );

        let plot_window = PlotWindow::new(arg, values.clone(), t_result, y_result);

        assert_eq!(plot_window.values.len(), 3);
        assert_eq!(plot_window.values[0], "CO");
        assert_eq!(plot_window.y_result.ncols(), 3);
        assert_eq!(plot_window.t_result.len(), 5);
    }

    #[test]
    fn test_plot_window_rejects_mesh_solution_mismatch() {
        let plot_window = PlotWindow::new(
            "time".to_string(),
            vec!["temperature".to_string()],
            DVector::from_vec(vec![0.0, 1.0, 2.0, 3.0]),
            DMatrix::from_row_slice(3, 1, &[300.0, 350.0, 400.0]),
        );

        let error = plot_window
            .validate_dimensions()
            .expect_err("mismatched mesh and solution sizes must be rejected");
        assert!(error.contains("mesh has 4 points"));
        assert!(error.contains("solution has 3 rows"));
    }

    #[test]
    fn test_plot_window_rejects_label_column_mismatch() {
        let plot_window = PlotWindow::new(
            "time".to_string(),
            vec!["temperature".to_string(), "pressure".to_string()],
            DVector::from_vec(vec![0.0, 1.0, 2.0]),
            DMatrix::from_row_slice(3, 1, &[300.0, 350.0, 400.0]),
        );

        let error = plot_window
            .validate_dimensions()
            .expect_err("mismatched label and column counts must be rejected");
        assert!(error.contains("2 variable labels"));
        assert!(error.contains("solution has 1 columns"));
    }

    #[test]
    fn test_save_plot_to_png_rejects_series_length_mismatch() {
        let error = PlotWindow::save_plot_to_png(
            Path::new("plots/temperature.png"),
            "temperature",
            "time",
            &[0.0, 1.0],
            &[300.0],
        )
        .expect_err("save_plot_to_png must reject mismatched series lengths");

        assert!(error.to_string().contains("mesh has 2 points"));
        assert!(error.to_string().contains("series has 1 values"));
    }

    #[test]
    fn test_plot_window_export_records_persistent_status() {
        let temp_dir = tempfile::tempdir().expect("temp dir should be available");
        let mut plot_window = PlotWindow::new(
            "time".to_string(),
            vec!["temperature".to_string()],
            DVector::from_vec(vec![0.0, 1.0]),
            DMatrix::from_row_slice(2, 1, &[300.0, 350.0]),
        );
        plot_window.output_dir = temp_dir.path().to_path_buf();

        let saved_path = plot_window
            .export_plot_series("temperature", &[0.0, 1.0], &[300.0, 350.0])
            .expect("plot export should succeed for aligned series");

        assert!(saved_path.exists());
        let status = plot_window
            .last_export_status
            .as_ref()
            .expect("plot export must persist its status");
        assert!(!status.is_error());
        assert!(status.status_text().contains("Saved plot to"));
        assert!(status.status_text().contains("temperature.png"));
    }
}

/// Visual test function - run with `cargo run --example plot_test`
pub fn run_plot_test() -> Result<(), eframe::Error> {
    let options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default()
            .with_inner_size([800.0, 600.0])
            .with_title("Plot Test Window"),
        ..Default::default()
    };

    eframe::run_native(
        "Plot Test",
        options,
        Box::new(|_cc| Ok(Box::new(PlotTestApp::new()))),
    )
}

struct PlotTestApp {
    plot_window: PlotWindow,
}

impl PlotTestApp {
    fn new() -> Self {
        // Create test data with sine and cosine waves
        let x_points: Vec<f64> = (0..100).map(|i| i as f64 * 0.1).collect();
        let t_result = DVector::from_vec(x_points.clone());

        let mut y_data = Vec::new();
        for &x in &x_points {
            y_data.push(x.sin()); // sine wave
            y_data.push(x.cos()); // cosine wave
            y_data.push((x * 0.5).sin() * 0.5); // slower sine wave
        }

        let y_result = DMatrix::from_row_slice(x_points.len(), 3, &y_data);
        let values = vec![
            "sin(x)".to_string(),
            "cos(x)".to_string(),
            "0.5*sin(0.5x)".to_string(),
        ];

        let plot_window = PlotWindow::new("x".to_string(), values, t_result, y_result);

        Self { plot_window }
    }
}

impl eframe::App for PlotTestApp {
    fn ui(&mut self, ui: &mut egui::Ui, _frame: &mut eframe::Frame) {
        self.plot_window.show(ui.ctx());
    }
}
