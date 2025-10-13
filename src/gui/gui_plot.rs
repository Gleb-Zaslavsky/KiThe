use eframe::egui;
use egui_plot::{Line, Plot, PlotPoints};
use nalgebra::{DMatrix, DVector};
use plotters::prelude::*;
use std::path::Path;
#[derive(Debug)]
pub struct PlotWindow {
    pub visible: bool,
    pub arg: String,
    pub values: Vec<String>,
    pub t_result: DVector<f64>,
    pub y_result: DMatrix<f64>,
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
        }
    }

    pub fn show(&mut self, ctx: &egui::Context) {
        if !self.visible {
            return;
        }

        egui::Window::new("Simulation Results")
            .open(&mut self.visible)
            .resizable(true)
            .vscroll(true)
            .show(ctx, |ui| {
                ui.heading(format!("Profiles of {}", self.arg));
                ui.separator();

                let x_data: Vec<f64> = self.t_result.iter().cloned().collect();

                for (col, name) in self.values.iter().enumerate() {
                    if col < self.y_result.ncols() {
                        ui.label(format!("Plot: {}", name));

                        Plot::new(format!("plot_{}", col))
                            .height(200.0)
                            .show(ui, |plot_ui| {
                                let y_col = self.y_result.column(col);
                                let points: PlotPoints = x_data
                                    .iter()
                                    .zip(y_col.iter())
                                    .map(|(&x, &y)| [x, y])
                                    .collect();

                                let line = Line::new(name.clone(), points).name(name.clone());
                                plot_ui.line(line);
                            });

                        ui.horizontal(|ui| {
                            ui.add_space(20.0);
                            if ui.button("💾 Save plot").clicked() {
                                if let Err(err) = Self::save_plot_to_png(
                                    name,
                                    &self.arg,
                                    &x_data,
                                    self.y_result.column(col).as_slice(),
                                ) {
                                    eprintln!("Failed to save plot {}: {}", name, err);
                                } else {
                                    ui.label(format!("Saved as {}.png", name));
                                }
                            }
                        });

                        ui.add_space(15.0);
                        ui.separator();
                    }
                }
            });
    }

    fn save_plot_to_png(
        var_name: &str,
        arg: &str,
        x: &[f64],
        y: &[f64],
    ) -> Result<(), Box<dyn std::error::Error>> {
        // Create output directory if not exists
        std::fs::create_dir_all("plots")?;

        let filename = format!("plots/{}.png", var_name.replace('/', "_"));
        let path = Path::new(&filename);

        let root = BitMapBackend::new(&path, (800, 600)).into_drawing_area();
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

        println!("Saved plot to {:?}", path);
        Ok(())
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
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        self.plot_window.show(ctx);
    }
}
