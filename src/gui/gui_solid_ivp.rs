use crate::Kinetics::solid_state_kinetics_IVP::{KineticModelIVP, KineticModelNames};
use crate::gui::gui_plot::PlotWindow;
use RustedSciThe::numerical::ODE_api2::SolverType;
use RustedSciThe::numerical::Radau::Radau_main::RadauOrder;
use egui_extras::{Column, TableBuilder};
use nalgebra::{DMatrix, DVector};
use strum::IntoEnumIterator;
use strum_macros::EnumIter;

#[derive(Debug, Clone, PartialEq, EnumIter)]
pub enum SolverTypeShotcut {
    RadauOrder3,
    RadauOrder5,
    RadauOrder7,
    NonstiffRK45,
    NonStiffDoPri,
    NonStiffAB4,
    BDF,
    BackwardEuler,
}
impl SolverTypeShotcut {
    pub fn get_solvertype(&self) -> SolverType {
        match self {
            Self::RadauOrder3 => SolverType::Radau(RadauOrder::Order3),
            Self::RadauOrder5 => SolverType::Radau(RadauOrder::Order5),
            Self::RadauOrder7 => SolverType::Radau(RadauOrder::Order7),
            Self::NonstiffRK45 => SolverType::NonStiff("RK45".to_string()),
            Self::NonStiffDoPri => SolverType::NonStiff("DOPRI".to_string()),
            Self::NonStiffAB4 => SolverType::NonStiff("AB4".to_string()),
            Self::BDF => SolverType::BDF,
            Self::BackwardEuler => SolverType::BackwardEuler,
        }
    }
}
// Local structure to hold the task data and results
pub struct SolidIVPApp {
    pub selected_solver: SolverTypeShotcut,
    pub selected_model: Option<KineticModelNames>,
    pub required_params: Vec<String>,
    t_final_string: String,
    beta_string: String,
    t0_string: String,
    e_string: String,
    a_string: String,
    pub t_final: f64,
    pub beta: f64,
    pub t0: f64,
    pub e: f64,
    pub a: f64,
    params_string: String,
    params: Vec<f64>,
    t_result: DVector<f64>,
    y_result: DMatrix<f64>,
    plot_window: Option<PlotWindow>,
}
impl Default for SolidIVPApp {
    fn default() -> Self {
        Self::new()
    }
}

impl SolidIVPApp {
    pub fn new() -> Self {
        Self {
            selected_model: None,
            selected_solver: SolverTypeShotcut::BDF,
            required_params: vec![],
            t_final_string: String::new(),
            beta_string: String::new(),
            t0_string: String::new(),
            e_string: String::new(),
            a_string: String::new(),
            t_final: 0.0,
            beta: 0.0,
            t0: 0.0,
            e: 0.0,
            a: 0.0,
            params_string: String::new(),
            params: vec![],
            t_result: DVector::zeros(0),
            y_result: DMatrix::zeros(0, 0),

            plot_window: None,
        }
    }

    pub fn run_kinetic_model(&mut self) -> Result<(), String> {
        if let Some(model) = &self.selected_model {
            let solver = self.selected_solver.get_solvertype();
            let mut kinetic_model = KineticModelIVP::new(solver);
            // For parameterized models, we would need to handle parameters
            // For now, we'll use an empty vector for parameter-free models
            kinetic_model.set_problem(self.t_final, self.beta, self.t0, self.e, self.a)?;
            kinetic_model.set_model(model.clone(), self.params.clone())?;
            kinetic_model.check_task()?;
            kinetic_model.solve()?;
            let (t_result, y_result) = kinetic_model.get_result()?;
            self.t_result = t_result;
            self.y_result = y_result;

            Ok(())
        } else {
            Err("No model selected".to_string())
        }
    }

    pub fn select_model_clicked(&mut self, model: KineticModelNames) {
        let required_params = model.required_params();
        self.required_params = required_params;
        self.selected_model = Some(model);
    }

    pub fn show(&mut self, ctx: &egui::Context, open: &mut bool) {
        // Configure text styles for better contrast and size
        let mut style = (*ctx.style()).clone();
        style.text_styles.insert(
            egui::TextStyle::Body,
            egui::FontId::new(16.0, egui::FontFamily::Proportional),
        );
        style.text_styles.insert(
            egui::TextStyle::Button,
            egui::FontId::new(16.0, egui::FontFamily::Proportional),
        );
        style.text_styles.insert(
            egui::TextStyle::Heading,
            egui::FontId::new(24.0, egui::FontFamily::Proportional),
        );
        ctx.set_style(style);

        egui::Window::new("Solid state kinetic models")
            .open(open)
            .default_size([1200.0, 800.0])
            .show(ctx, |ui| {
                // Table with 3 columns: button, model name, formula
                TableBuilder::new(ui)
                    .column(Column::auto().resizable(true)) // Button column
                    .column(Column::auto().resizable(true)) // Model name column
                    .column(Column::remainder()) // Formula column
                    .header(20.0, |mut header| {
                        header.col(|ui| {
                            ui.heading("Select");
                        });
                        header.col(|ui| {
                            ui.heading("Model name");
                        });
                        header.col(|ui| {
                            ui.heading("Formula");
                        });
                    })
                    .body(|mut body| {
                        for model in KineticModelNames::iter() {
                            body.row(20.0, |mut row| {
                                row.col(|ui| {
                                    if ui.button("Select").clicked() {
                                        self.select_model_clicked(model.clone());
                                    }
                                });
                                row.col(|ui| {
                                    ui.label(format!("{:?}", model));
                                });
                                row.col(|ui| {
                                    ui.label(model.map_of_names_and_formulas());
                                });
                            });
                        }
                    });

                ui.separator();

                egui::ComboBox::from_label("Select Solver")
                    .selected_text(format!("{:?}", self.selected_solver))
                    .show_ui(ui, |ui| {
                        for solver in SolverTypeShotcut::iter() {
                            ui.selectable_value(
                                &mut self.selected_solver,
                                solver.clone(),
                                format!("{:?}", solver),
                            );
                        }
                    });

                ui.separator();

                ui.horizontal(|ui| {
                    if self.params.len() != self.required_params.len() {
                        self.params = vec![0.0; self.required_params.len()];
                    }
                    for (idx, param) in self.required_params.iter().enumerate() {
                        ui.label(format!(" {}", param));

                        ui.text_edit_singleline(&mut self.params_string);
                        if let Ok(parsed) = self.params_string.trim().parse::<f64>() {
                            self.params[idx] = parsed;
                        }
                        ui.label(format!("params:{:?}", self.params));
                    }
                });

                // Input fields for t_final, beta, T0, E, A
                ui.separator();
                ui.horizontal(|ui| {
                    ui.label("t_final:");
                    ui.text_edit_singleline(&mut self.t_final_string);
                    if let Ok(parsed) = self.t_final_string.trim().parse::<f64>() {
                        self.t_final = parsed;
                    }
                    ui.label(self.t_final.to_string());
                });
                ui.horizontal(|ui| {
                    ui.label("beta:");
                    ui.text_edit_singleline(&mut self.beta_string);
                    if let Ok(parsed) = self.beta_string.trim().parse::<f64>() {
                        self.beta = parsed;
                    }
                    ui.label(self.beta.to_string());
                });
                ui.horizontal(|ui| {
                    ui.label("T0:");
                    ui.text_edit_singleline(&mut self.t0_string);
                    if let Ok(parsed) = self.t0_string.trim().parse::<f64>() {
                        self.t0 = parsed;
                    }
                    ui.label(self.t0.to_string());
                });
                ui.horizontal(|ui| {
                    ui.label("E:");
                    ui.text_edit_singleline(&mut self.e_string);
                    if let Ok(parsed) = self.e_string.trim().parse::<f64>() {
                        self.e = parsed;
                    }
                    ui.label(self.e.to_string());
                });
                ui.horizontal(|ui| {
                    ui.label("A:");
                    ui.text_edit_singleline(&mut self.a_string);
                    if let Ok(parsed) = self.a_string.trim().parse::<f64>() {
                        self.a = parsed;
                    }
                    ui.label(self.a.to_string());
                });
                ui.separator();
                ui.add_space(10.0);
                let text = format!(
                    "Model: {:?}, \n
                       A={}, E={}, beta={}, params = {:?}, solver {:?} ",
                    self.selected_model,
                    self.a,
                    self.e,
                    self.beta,
                    self.params,
                    self.selected_solver
                );
                ui.label(text);
                ui.add_space(10.0);
                ui.separator();
                // Update button
                if ui.button("Run Model").clicked() {
                    match self.run_kinetic_model() {
                        Ok(_) => {
                            println!("integration successful");
                            println!("y={}", self.y_result);
                            let mut plot = PlotWindow::new(
                                "t".to_string(),
                                vec!["a".to_string()],
                                self.t_result.clone(),
                                self.y_result.clone(),
                            );
                            plot.show(ctx);
                            self.plot_window = Some(plot);
                        }
                        Err(e) => {
                            println!("Error solving model: {}", e);
                        }
                    }
                }
            });
    }
}
