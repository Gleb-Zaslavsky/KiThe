//! Controller code for Kinetic Methods menu and sub-windows.
//!
//! Each menu item opens a dedicated popup window that currently contains a placeholder
//! for future implementation.

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

        if state.criado_master_curve_open {
            egui::Window::new("Criado Master Curve")
                .open(&mut state.criado_master_curve_open)
                .default_size([520.0, 340.0])
                .show(ui.ctx(), |ui| {
                    ui.label("TODO: Implement Criado master curve analysis.");
                });
        }

        if state.fit_model_open {
            egui::Window::new("Fit Model")
                .open(&mut state.fit_model_open)
                .default_size([480.0, 320.0])
                .show(ui.ctx(), |ui| {
                    ui.label("TODO: Implement model fitting.");
                });
        }

        Ok(())
    }
}
