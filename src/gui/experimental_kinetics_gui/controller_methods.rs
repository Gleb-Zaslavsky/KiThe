//! Controller code for Kinetic Methods menu and sub-windows.
//!
//! Each menu item opens a dedicated popup window for a kinetic-analysis workflow.

use crate::Kinetics::experimental_kinetics::fitting::{
    FitColumnRequest, FitColumnResult, FittingModelName,
};
use crate::Kinetics::experimental_kinetics::kinetic_methods::KineticMethod;
use crate::Kinetics::experimental_kinetics::kinetic_methods::Kissinger::Kissinger;
use crate::Kinetics::experimental_kinetics::kinetic_methods::combined::CombinedKineticAnalysis;
use crate::Kinetics::experimental_kinetics::kinetic_methods::is_this_a_sublimation::{
    IsThisASublimationResult, SublimationMethod,
};
use crate::Kinetics::experimental_kinetics::kinetic_methods::isoconversion::{
    IsoconversionalKineticMethod, IsoconversionalMethod,
};
use crate::gui::experimental_kinetics_gui::model::{PlotModel, TGAGUIError};
use eframe::egui;

/// State for which Kinetic Methods tool windows are currently open.
pub struct KineticMethodsWindowState {
    pub isoconversional_open: bool,
    pub isoconversional_method: IsoconversionalMethod,
    pub isoconversional_status: String,
    pub isoconversional_use_selected_experiments: bool,
    pub isoconversional_selected_experiments: Vec<String>,
    pub isoconversional_output_id: String,
    pub kissinger_open: bool,
    pub kissinger_status: String,
    pub kissinger_ea_kj: Option<f64>,
    pub kissinger_r2: Option<f64>,
    pub kissinger_kinetic_expression: String,
    pub combined_kinetics_open: bool,
    pub combined_status: String,
    pub combined_n_min: f64,
    pub combined_n_max: f64,
    pub combined_n_steps: usize,
    pub combined_m_min: f64,
    pub combined_m_max: f64,
    pub combined_m_steps: usize,
    pub combined_eta_min: f64,
    pub combined_eta_max: f64,
    pub combined_refinement_steps: usize,
    pub combined_result_n: Option<f64>,
    pub combined_result_m: Option<f64>,
    pub combined_result_ea_kj: Option<f64>,
    pub combined_result_r2: Option<f64>,
    pub criado_master_curve_open: bool,
    pub fit_model_open: bool,
    pub fit_model_experiment_id: String,
    pub fit_model_x_col: String,
    pub fit_model_y_col: String,
    pub fit_model_output_col: String,
    pub fit_model_model: FittingModelName,
    pub fit_model_polynomial_degree: usize,
    pub fit_model_initial_guess: String,
    pub fit_model_tolerance: f64,
    pub fit_model_max_iter: usize,
    pub fit_model_status: String,
    pub fit_model_result: Option<FitColumnResult>,
    pub sublimation_check_open: bool,
    pub sublimation_status: String,
    pub sublimation_mean_ea_kj: Option<f64>,
    pub sublimation_mean_k: Option<f64>,
    pub sublimation_verdict: Option<String>,
    pub methods_golden_pipeline: bool,
}

impl Default for KineticMethodsWindowState {
    fn default() -> Self {
        Self {
            isoconversional_open: false,
            isoconversional_method: IsoconversionalMethod::default(),
            isoconversional_status: String::new(),
            isoconversional_use_selected_experiments: false,
            isoconversional_selected_experiments: Vec::new(),
            isoconversional_output_id: "isoconversional_result".to_string(),
            kissinger_open: false,
            kissinger_status: String::new(),
            kissinger_ea_kj: None,
            kissinger_r2: None,
            kissinger_kinetic_expression: String::new(),
            combined_kinetics_open: false,
            combined_status: String::new(),
            combined_n_min: CombinedKineticAnalysis::default().n_min,
            combined_n_max: CombinedKineticAnalysis::default().n_max,
            combined_n_steps: CombinedKineticAnalysis::default().n_steps,
            combined_m_min: CombinedKineticAnalysis::default().m_min,
            combined_m_max: CombinedKineticAnalysis::default().m_max,
            combined_m_steps: CombinedKineticAnalysis::default().m_steps,
            combined_eta_min: CombinedKineticAnalysis::default().eta_min,
            combined_eta_max: CombinedKineticAnalysis::default().eta_max,
            combined_refinement_steps: CombinedKineticAnalysis::default().refinement_steps,
            combined_result_n: None,
            combined_result_m: None,
            combined_result_ea_kj: None,
            combined_result_r2: None,
            criado_master_curve_open: false,
            fit_model_open: false,
            fit_model_experiment_id: String::new(),
            fit_model_x_col: String::new(),
            fit_model_y_col: String::new(),
            fit_model_output_col: String::new(),
            fit_model_model: FittingModelName::DecExp,
            fit_model_polynomial_degree: 2,
            fit_model_initial_guess: String::new(),
            fit_model_tolerance: 1e-6,
            fit_model_max_iter: 300,
            fit_model_status: String::new(),
            fit_model_result: None,
            sublimation_check_open: false,
            sublimation_status: String::new(),
            sublimation_mean_ea_kj: None,
            sublimation_mean_k: None,
            sublimation_verdict: None,
            methods_golden_pipeline: false,
        }
    }
}

pub struct KineticMethods;

impl KineticMethods {
    /// Render the dropdown menu entries under the "Kinetic Methods" top menu.
    pub fn show_menu(ui: &mut egui::Ui, state: &mut KineticMethodsWindowState) {
        if ui.button("Isoconversional").clicked() {
            state.isoconversional_open = true;
        }
        if ui.button("Kissinger").clicked() {
            state.kissinger_open = true;
        }
        if ui.button("Combined Kinetics Analysis").clicked() {
            state.combined_kinetics_open = true;
        }
        if ui.button("Criado Master Curve").clicked() {
            state.criado_master_curve_open = true;
        }
        if ui.button("Fit Model").clicked() {
            state.fit_model_open = true;
        }
        if ui.button("Is this a sublimation?").clicked() {
            state.sublimation_check_open = true;
        }
        if ui.button("Golden Pipeline").clicked() {
            state.methods_golden_pipeline = true;
        }
    }
    //===========================================================================================================
    /// ISOCONVERSIONAL ANALYSIS
    ///===========================================================================================================
    /// Show the isoconversional analysis window.
    pub fn show_isoconversional_window(
        ui: &mut egui::Ui,
        model: &mut PlotModel,
        state: &mut KineticMethodsWindowState,
    ) -> Result<(), TGAGUIError> {
        if !state.isoconversional_open {
            return Ok(());
        }

        egui::Window::new("Isoconversional")
            .open(&mut state.isoconversional_open)
            .default_size([540.0, 360.0])
            .show(ui.ctx(), |ui| {
                ui.vertical(|ui| {
                    ui.horizontal(|ui| {
                        ui.label("Method:");
                        egui::ComboBox::from_id_salt("isoconv_method")
                            .selected_text(state.isoconversional_method.display_name())
                            .show_ui(ui, |ui| {
                                for method in IsoconversionalMethod::all() {
                                    ui.selectable_value(
                                        &mut state.isoconversional_method,
                                        *method,
                                        method.display_name(),
                                    );
                                }
                            });
                    });

                    ui.add_space(6.0);

                    // Build a stable list of experiment ids for selection.
                    let exp_ids = model.list_of_experiments();
                    state
                        .isoconversional_selected_experiments
                        .retain(|id| exp_ids.iter().any(|exp| exp == id));

                    ui.horizontal(|ui| {
                        ui.label("Experiments:");
                        ui.checkbox(
                            &mut state.isoconversional_use_selected_experiments,
                            "Use selected only",
                        );
                    });

                    // When enabled, the user chooses a subset of experiments.
                    if state.isoconversional_use_selected_experiments {
                        if exp_ids.is_empty() {
                            ui.label("No experiments available.");
                        } else {
                            ui.vertical(|ui| {
                                for id in &exp_ids {
                                    let mut selected = state
                                        .isoconversional_selected_experiments
                                        .iter()
                                        .any(|v| v == id);
                                    if ui.checkbox(&mut selected, id).changed() {
                                        if selected {
                                            if !state
                                                .isoconversional_selected_experiments
                                                .contains(id)
                                            {
                                                state
                                                    .isoconversional_selected_experiments
                                                    .push(id.clone());
                                            }
                                        } else {
                                            state
                                                .isoconversional_selected_experiments
                                                .retain(|v| v != id);
                                        }
                                    }
                                }
                            });
                        }
                    }

                    ui.add_space(6.0);

                    ui.horizontal(|ui| {
                        ui.label("Output ID:");
                        ui.text_edit_singleline(&mut state.isoconversional_output_id);
                    });

                    ui.add_space(6.0);

                    ui.horizontal(|ui| {
                        ui.label("Status:");
                        ui.label(&state.isoconversional_status);
                    });

                    ui.add_space(8.0);

                    // Run the chosen method and store the result as a new experiment.
                    if ui.button("Calculate").clicked() {
                        let status = (|| -> Result<String, TGAGUIError> {
                            if exp_ids.is_empty() {
                                return Ok("No experiments to analyze.".to_string());
                            }

                            let selected_ids: Vec<&str> = state
                                .isoconversional_selected_experiments
                                .iter()
                                .map(|s| s.as_str())
                                .collect();

                            let what_exp_to_take = if state.isoconversional_use_selected_experiments
                            {
                                if selected_ids.is_empty() {
                                    return Ok("Select at least one experiment.".to_string());
                                }
                                Some(selected_ids.as_slice())
                            } else {
                                None
                            };

                            let output_id = if state.isoconversional_output_id.trim().is_empty() {
                                "isoconversional_result".to_string()
                            } else {
                                state.isoconversional_output_id.trim().to_string()
                            };

                            let view = model.create_kinetic_data_view_for_method(
                                what_exp_to_take,
                                &state.isoconversional_method,
                            )?;

                            let method = IsoconversionalKineticMethod {
                                method: state.isoconversional_method.clone(),
                            };

                            method.check_input(&view)?;
                            let result = method.compute(&view)?;

                            model.push_isoconversional_result(&result, &output_id)?;

                            Ok(format!("Done. Result saved as '{}'.", output_id))
                        })();

                        state.isoconversional_status = match status {
                            Ok(msg) => msg,
                            Err(err) => format!("Error: {:?}", err),
                        };
                    }
                });
            });

        Ok(())
    }
    //===========================================================================================================
    /// KISSINGER ANALYSIS
    //===========================================================================================================
    /// Show the Kissinger analysis window.
    pub fn show_kissinger_window(
        ui: &mut egui::Ui,
        model: &mut PlotModel,
        state: &mut KineticMethodsWindowState,
    ) -> Result<(), TGAGUIError> {
        if !state.kissinger_open {
            return Ok(());
        }

        egui::Window::new("Kissinger")
            .open(&mut state.kissinger_open)
            .default_size([520.0, 300.0])
            .show(ui.ctx(), |ui| {
                ui.vertical(|ui| {
                    let ea_label = state
                        .kissinger_ea_kj
                        .map(|v| format!("{:.2}", v))
                        .unwrap_or_else(|| "-".to_string());
                    let r2_label = state
                        .kissinger_r2
                        .map(|v| format!("{:.4}", v))
                        .unwrap_or_else(|| "-".to_string());

                    ui.horizontal(|ui| {
                        ui.label("Ea (kJ/mol):");
                        ui.label(ea_label);
                    });

                    ui.horizontal(|ui| {
                        ui.label("R2:");
                        ui.label(r2_label);
                    });

                    ui.add_space(6.0);

                    ui.horizontal(|ui| {
                        ui.label("System messages:");
                        ui.label(&state.kissinger_status);
                    });

                    ui.add_space(8.0);

                    // Placeholder for future work on pre-exponential factor estimation.
                    ui.label("Kinetic expression:");
                    ui.text_edit_singleline(&mut state.kissinger_kinetic_expression);
                    if ui.button("Calculate kinetic expression").clicked() {
                        todo!("Implement pre-exponential factor estimation for Kissinger method.");
                    }

                    if ui.button("Calculate pre-exponential factor").clicked() {
                        todo!("Implement pre-exponential factor estimation for Kissinger method.");
                    }

                    ui.add_space(8.0);

                    // Run the Kissinger calculation without creating new series data or curves.
                    if ui.button("Calculate").clicked() {
                        let status = (|| -> Result<String, TGAGUIError> {
                            // Kissinger needs the same columns as non-isothermal isoconversional methods.
                            let view = model.create_kinetic_data_view_for_method(
                                None,
                                &IsoconversionalMethod::OFW,
                            )?;

                            let result = Kissinger.compute(&view)?;

                            let ea_kj = result.ea / 1000.0;
                            state.kissinger_ea_kj = Some(ea_kj);
                            state.kissinger_r2 = Some(result.regression.r2);

                            Ok("Done.".to_string())
                        })();

                        state.kissinger_status = match status {
                            Ok(msg) => msg,
                            Err(err) => format!("Error: {:?}", err),
                        };
                    }
                });
            });

        Ok(())
    }
    //===========================================================================================================
    /// SUBLIMATION CHECK
    //===========================================================================================================
    /// Show the "Is this a sublimation?" window.
    pub fn show_sublimation_window(
        ui: &mut egui::Ui,
        model: &mut PlotModel,
        state: &mut KineticMethodsWindowState,
    ) -> Result<(), TGAGUIError> {
        if !state.sublimation_check_open {
            return Ok(());
        }

        egui::Window::new("Is this a sublimation?")
            .open(&mut state.sublimation_check_open)
            .default_size([560.0, 320.0])
            .show(ui.ctx(), |ui| {
                ui.vertical(|ui| {
                    let ea_label = state
                        .sublimation_mean_ea_kj
                        .map(|v| format!("{:.2}", v))
                        .unwrap_or_else(|| "-".to_string());
                    let k_label = state
                        .sublimation_mean_k
                        .map(|v| format!("{:.4e}", v))
                        .unwrap_or_else(|| "-".to_string());
                    let verdict_label = state
                        .sublimation_verdict
                        .clone()
                        .unwrap_or_else(|| "-".to_string());

                    ui.horizontal(|ui| {
                        ui.label("Mean Ea (kJ/mol):");
                        ui.label(ea_label);
                    });

                    ui.horizontal(|ui| {
                        ui.label("Mean k:");
                        ui.label(k_label);
                    });

                    ui.horizontal(|ui| {
                        ui.label("Verdict:");
                        ui.label(verdict_label);
                    });

                    ui.add_space(6.0);

                    ui.horizontal(|ui| {
                        ui.label("System messages:");
                        ui.label(&state.sublimation_status);
                    });

                    ui.add_space(8.0);

                    // Run the solver, store mean values, and push fitted curves into the series.
                    if ui.button("Calculate").clicked() {
                        let status = (|| -> Result<String, TGAGUIError> {
                            let method = SublimationMethod::default();
                            let cols = method.required_columns_by_nature();
                            let view = model.create_kinetic_data_view(None, cols)?;

                            let results = method.compute(&view)?;
                            let (mean_ea, mean_k, verdict) =
                                SublimationMethod::mean_results(&results)?;

                            state.sublimation_mean_ea_kj = Some(mean_ea / 1000.0);
                            state.sublimation_mean_k = Some(mean_k);

                            let verdict_text = match verdict {
                                IsThisASublimationResult::Yes(r2) => {
                                    format!("Yes (R2={:.4})", r2)
                                }
                                IsThisASublimationResult::MayBe(r2) => {
                                    format!("Maybe (R2={:.4})", r2)
                                }
                                IsThisASublimationResult::No(r2) => {
                                    format!("No (R2={:.4})", r2)
                                }
                            };
                            state.sublimation_verdict = Some(verdict_text);

                            model.push_sublimation_fitted_rates(&results)?;

                            Ok("Done.".to_string())
                        })();

                        state.sublimation_status = match status {
                            Ok(msg) => msg,
                            Err(err) => format!("Error: {:?}", err),
                        };
                    }
                });
            });

        Ok(())
    }
    //===========================================================================================================
    /// COMBINED KINETICS ANALYSIS
    //===========================================================================================================
    /// Show the Combined Kinetics Analysis window.
    pub fn show_combined_kinetics_window(
        ui: &mut egui::Ui,
        model: &mut PlotModel,
        state: &mut KineticMethodsWindowState,
    ) -> Result<(), TGAGUIError> {
        if !state.combined_kinetics_open {
            return Ok(());
        }

        egui::Window::new("Combined Kinetics Analysis")
            .open(&mut state.combined_kinetics_open)
            .default_size([560.0, 420.0])
            .show(ui.ctx(), |ui| {
                ui.vertical(|ui| {
                    // Manual parameter inputs for the analysis grid.
                    ui.horizontal(|ui| {
                        ui.label("n_min:");
                        ui.add(egui::DragValue::new(&mut state.combined_n_min).speed(0.1));
                        ui.label("n_max:");
                        ui.add(egui::DragValue::new(&mut state.combined_n_max).speed(0.1));
                        ui.label("n_steps:");
                        ui.add(egui::DragValue::new(&mut state.combined_n_steps).speed(1.0));
                    });

                    ui.horizontal(|ui| {
                        ui.label("m_min:");
                        ui.add(egui::DragValue::new(&mut state.combined_m_min).speed(0.1));
                        ui.label("m_max:");
                        ui.add(egui::DragValue::new(&mut state.combined_m_max).speed(0.1));
                        ui.label("m_steps:");
                        ui.add(egui::DragValue::new(&mut state.combined_m_steps).speed(1.0));
                    });

                    ui.horizontal(|ui| {
                        ui.label("alpha_min:");
                        ui.add(egui::DragValue::new(&mut state.combined_eta_min).speed(0.01));
                        ui.label("alpha_max:");
                        ui.add(egui::DragValue::new(&mut state.combined_eta_max).speed(0.01));
                        ui.label("refinement_steps:");
                        ui.add(
                            egui::DragValue::new(&mut state.combined_refinement_steps).speed(1.0),
                        );
                    });

                    ui.add_space(6.0);

                    let n_label = state
                        .combined_result_n
                        .map(|v| format!("{:.4}", v))
                        .unwrap_or_else(|| "-".to_string());
                    let m_label = state
                        .combined_result_m
                        .map(|v| format!("{:.4}", v))
                        .unwrap_or_else(|| "-".to_string());
                    let ea_label = state
                        .combined_result_ea_kj
                        .map(|v| format!("{:.2}", v))
                        .unwrap_or_else(|| "-".to_string());
                    let r2_label = state
                        .combined_result_r2
                        .map(|v| format!("{:.4}", v))
                        .unwrap_or_else(|| "-".to_string());

                    ui.horizontal(|ui| {
                        ui.label("n:");
                        ui.label(n_label);
                        ui.label("m:");
                        ui.label(m_label);
                    });

                    ui.horizontal(|ui| {
                        ui.label("Ea (kJ/mol):");
                        ui.label(ea_label);
                        ui.label("R2:");
                        ui.label(r2_label);
                    });

                    ui.add_space(6.0);

                    ui.horizontal(|ui| {
                        ui.label("System messages:");
                        ui.label(&state.combined_status);
                    });

                    ui.add_space(8.0);

                    // Run the analysis and update result labels.
                    if ui.button("Calculate").clicked() {
                        let status = (|| -> Result<String, TGAGUIError> {
                            let analysis = CombinedKineticAnalysis {
                                n_min: state.combined_n_min,
                                n_max: state.combined_n_max,
                                n_steps: state.combined_n_steps,
                                m_min: state.combined_m_min,
                                m_max: state.combined_m_max,
                                m_steps: state.combined_m_steps,
                                eta_min: state.combined_eta_min,
                                eta_max: state.combined_eta_max,
                                refinement_steps: state.combined_refinement_steps,
                            };

                            let cols = analysis.required_columns_by_nature();
                            let view = model.create_kinetic_data_view(None, cols)?;
                            let result = analysis.compute(&view)?;

                            state.combined_result_n = Some(result.n);
                            state.combined_result_m = Some(result.m);
                            state.combined_result_ea_kj = Some(result.ea / 1000.0);
                            state.combined_result_r2 = Some(result.regression.r2);

                            Ok("Done.".to_string())
                        })();

                        state.combined_status = match status {
                            Ok(msg) => msg,
                            Err(err) => format!("Error: {:?}", err),
                        };
                    }
                });
            });

        Ok(())
    }

    /// Show the column-fitting window.
    pub fn show_fit_model_window(
        ui: &mut egui::Ui,
        model: &mut PlotModel,
        state: &mut KineticMethodsWindowState,
    ) -> Result<(), TGAGUIError> {
        if !state.fit_model_open {
            return Ok(());
        }

        egui::Window::new("Fit Model")
            .open(&mut state.fit_model_open)
            .default_size([620.0, 430.0])
            .show(ui.ctx(), |ui| {
                ui.vertical(|ui| {
                    let exp_ids = model.list_of_experiments();
                    if state.fit_model_experiment_id.is_empty()
                        || !exp_ids
                            .iter()
                            .any(|id| id == &state.fit_model_experiment_id)
                    {
                        state.fit_model_experiment_id =
                            exp_ids.first().cloned().unwrap_or_default();
                    }

                    ui.horizontal(|ui| {
                        ui.label("Experiment:");
                        egui::ComboBox::from_id_salt("fit_model_experiment")
                            .selected_text(if state.fit_model_experiment_id.is_empty() {
                                "-"
                            } else {
                                &state.fit_model_experiment_id
                            })
                            .show_ui(ui, |ui| {
                                for id in &exp_ids {
                                    ui.selectable_value(
                                        &mut state.fit_model_experiment_id,
                                        id.clone(),
                                        id,
                                    );
                                }
                            });
                    });

                    let columns = if state.fit_model_experiment_id.is_empty() {
                        Ok(Vec::new())
                    } else {
                        model.list_of_columns(&state.fit_model_experiment_id)
                    };

                    match columns {
                        Ok(columns) => {
                            if state.fit_model_x_col.is_empty()
                                || !columns.iter().any(|col| col == &state.fit_model_x_col)
                            {
                                state.fit_model_x_col = columns
                                    .iter()
                                    .find(|c| c.as_str() == "time")
                                    .cloned()
                                    .or_else(|| columns.first().cloned())
                                    .unwrap_or_default();
                            }
                            if state.fit_model_y_col.is_empty()
                                || !columns.iter().any(|col| col == &state.fit_model_y_col)
                            {
                                state.fit_model_y_col = columns
                                    .iter()
                                    .find(|c| c.as_str() == "mass")
                                    .cloned()
                                    .or_else(|| {
                                        columns
                                            .iter()
                                            .find(|c| c.as_str() != state.fit_model_x_col)
                                            .cloned()
                                    })
                                    .or_else(|| columns.first().cloned())
                                    .unwrap_or_default();
                            }

                            ui.horizontal(|ui| {
                                ui.label("X:");
                                egui::ComboBox::from_id_salt("fit_model_x")
                                    .selected_text(if state.fit_model_x_col.is_empty() {
                                        "-"
                                    } else {
                                        &state.fit_model_x_col
                                    })
                                    .show_ui(ui, |ui| {
                                        for col in &columns {
                                            ui.selectable_value(
                                                &mut state.fit_model_x_col,
                                                col.clone(),
                                                col,
                                            );
                                        }
                                    });

                                ui.label("Y:");
                                egui::ComboBox::from_id_salt("fit_model_y")
                                    .selected_text(if state.fit_model_y_col.is_empty() {
                                        "-"
                                    } else {
                                        &state.fit_model_y_col
                                    })
                                    .show_ui(ui, |ui| {
                                        for col in &columns {
                                            ui.selectable_value(
                                                &mut state.fit_model_y_col,
                                                col.clone(),
                                                col,
                                            );
                                        }
                                    });
                            });

                            ui.horizontal(|ui| {
                                ui.label("Model:");
                                egui::ComboBox::from_id_salt("fit_model_kind")
                                    .selected_text(state.fit_model_model.display_name())
                                    .show_ui(ui, |ui| {
                                        for fit_model in FittingModelName::all_for_gui() {
                                            let changed = ui
                                                .selectable_value(
                                                    &mut state.fit_model_model,
                                                    fit_model,
                                                    fit_model.display_name(),
                                                )
                                                .changed();
                                            if changed {
                                                if let Some(n) = fit_model.polynomial_degree() {
                                                    state.fit_model_polynomial_degree = n;
                                                }
                                            }
                                        }
                                    });
                            });

                            if state.fit_model_model.polynomial_degree().is_some() {
                                ui.horizontal(|ui| {
                                    ui.label("Polynomial degree:");
                                    let changed = ui
                                        .add(
                                            egui::DragValue::new(
                                                &mut state.fit_model_polynomial_degree,
                                            )
                                            .speed(1.0)
                                            .range(1..=12),
                                        )
                                        .changed();
                                    if changed
                                        || state.fit_model_model.polynomial_degree()
                                            != Some(state.fit_model_polynomial_degree)
                                    {
                                        state.fit_model_model = FittingModelName::Polynom {
                                            n: state.fit_model_polynomial_degree,
                                        };
                                    }
                                });
                            }

                            ui.horizontal(|ui| {
                                ui.label("Formula:");
                                ui.monospace(state.fit_model_model.equation_str());
                            });

                            egui::CollapsingHeader::new("Model equations")
                                .default_open(true)
                                .show(ui, |ui| {
                                    egui::Grid::new("fit_model_equation_table")
                                        .striped(true)
                                        .show(ui, |ui| {
                                            ui.strong("Model");
                                            ui.strong("Formula");
                                            ui.end_row();

                                            for mut fit_model in FittingModelName::all_for_gui() {
                                                if fit_model.polynomial_degree().is_some() {
                                                    fit_model = FittingModelName::Polynom {
                                                        n: state.fit_model_polynomial_degree,
                                                    };
                                                }
                                                ui.label(fit_model.display_name());
                                                ui.monospace(fit_model.equation_str());
                                                ui.end_row();
                                            }
                                        });
                                });

                            ui.horizontal(|ui| {
                                ui.label("Output column:");
                                ui.text_edit_singleline(&mut state.fit_model_output_col);
                            });

                            ui.horizontal(|ui| {
                                ui.label("Tolerance:");
                                ui.add(
                                    egui::DragValue::new(&mut state.fit_model_tolerance)
                                        .speed(1e-7)
                                        .range(1e-12..=1e-1),
                                );
                                ui.label("Max iter:");
                                ui.add(
                                    egui::DragValue::new(&mut state.fit_model_max_iter)
                                        .speed(10.0)
                                        .range(1..=100_000),
                                );
                            });

                            ui.horizontal(|ui| {
                                ui.label("Initial guess:");
                                ui.text_edit_singleline(&mut state.fit_model_initial_guess);
                            });

                            ui.horizontal(|ui| {
                                ui.label("Status:");
                                ui.label(&state.fit_model_status);
                            });

                            ui.add_space(8.0);
                            if ui.button("Fit").clicked() {
                                let status = (|| -> Result<String, TGAGUIError> {
                                    if state.fit_model_experiment_id.is_empty()
                                        || state.fit_model_x_col.is_empty()
                                        || state.fit_model_y_col.is_empty()
                                    {
                                        return Err(TGAGUIError::BindingError(
                                            "Choose experiment, X and Y columns first.".to_string(),
                                        ));
                                    }

                                    let mut request = FitColumnRequest::new(
                                        state.fit_model_experiment_id.clone(),
                                        state.fit_model_x_col.clone(),
                                        state.fit_model_y_col.clone(),
                                        state.fit_model_model,
                                    )
                                    .with_tolerance(state.fit_model_tolerance)
                                    .with_max_iterations(state.fit_model_max_iter);

                                    let output = state.fit_model_output_col.trim();
                                    if !output.is_empty() {
                                        request = request.with_output_col(output.to_string());
                                    }

                                    if let Some(initial_guess) =
                                        parse_fit_initial_guess(&state.fit_model_initial_guess)?
                                    {
                                        request = request.with_initial_guess(initial_guess);
                                    }

                                    let result = model.fit_column_and_plot(request)?;
                                    let output_col = result.output_col.clone();
                                    state.fit_model_result = Some(result);
                                    Ok(format!("Done. Added '{}'.", output_col))
                                })();

                                state.fit_model_status = match status {
                                    Ok(msg) => msg,
                                    Err(err) => format!("Error: {:?}", err),
                                };
                            }

                            if let Some(result) = &state.fit_model_result {
                                ui.separator();
                                ui.horizontal(|ui| {
                                    ui.label("Backend:");
                                    ui.label(&result.method);
                                    ui.label("R2:");
                                    ui.label(format!("{:.6}", result.r2));
                                    ui.label("Output:");
                                    ui.label(&result.output_col);
                                });

                                egui::Grid::new("fit_model_coefficients")
                                    .striped(true)
                                    .show(ui, |ui| {
                                        ui.label("Coefficient");
                                        ui.label("Value");
                                        ui.end_row();
                                        for (name, value) in &result.coefficients {
                                            ui.label(name);
                                            ui.label(format!("{:.8e}", value));
                                            ui.end_row();
                                        }
                                    });
                            }
                        }
                        Err(err) => {
                            ui.label(format!("Cannot read columns: {:?}", err));
                        }
                    }
                });
            });

        Ok(())
    }

    //==============================================================================================
    // END OF KINETIC METHODS WINDOWS
    //==============================================================================================
    // COMBINED KINETICS ANALYSIS
    //===================================================================================================
    /// Render any open method windows.
    pub fn show_windows(
        ui: &mut egui::Ui,
        model: &mut PlotModel,
        state: &mut KineticMethodsWindowState,
    ) -> Result<(), TGAGUIError> {
        // Separate helper functions keep each window small and focused.
        Self::show_isoconversional_window(ui, model, state)?;
        Self::show_kissinger_window(ui, model, state)?;
        Self::show_combined_kinetics_window(ui, model, state)?;
        Self::show_sublimation_window(ui, model, state)?;
        Self::show_fit_model_window(ui, model, state)?;

        if state.criado_master_curve_open {
            egui::Window::new("Criado Master Curve")
                .open(&mut state.criado_master_curve_open)
                .default_size([520.0, 340.0])
                .show(ui.ctx(), |ui| {
                    ui.label("TODO: Implement Criado master curve analysis.");
                });
        }

        Ok(())
    }
}

fn parse_fit_initial_guess(raw: &str) -> Result<Option<Vec<f64>>, TGAGUIError> {
    let raw = raw.trim();
    if raw.is_empty() {
        return Ok(None);
    }

    let mut values = Vec::new();
    for token in raw
        .split(|ch: char| ch == ',' || ch == ';' || ch.is_whitespace())
        .filter(|token| !token.trim().is_empty())
    {
        let value = token.trim().parse::<f64>().map_err(|err| {
            TGAGUIError::BindingError(format!(
                "Cannot parse initial guess value '{}': {}",
                token, err
            ))
        })?;
        values.push(value);
    }

    Ok(Some(values))
}
