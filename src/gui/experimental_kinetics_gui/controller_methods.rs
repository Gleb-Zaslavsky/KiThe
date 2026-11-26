//! Controller code for Kinetic Methods menu and sub-windows.
//!
//! Each menu item opens a dedicated popup window for a kinetic-analysis workflow.

use crate::Kinetics::experimental_kinetics::kinetic_methods::KineticMethod;
use crate::Kinetics::experimental_kinetics::kinetic_methods::Kissinger::Kissinger;
use crate::Kinetics::experimental_kinetics::kinetic_methods::combined::CombinedKineticAnalysis;
use crate::Kinetics::experimental_kinetics::kinetic_methods::is_this_a_sublimation::{
    IsThisASublimationResult, SublimationMethod,
};
use crate::Kinetics::experimental_kinetics::kinetic_methods::isoconversion::{
    IsoconversionalKineticMethod, IsoconversionalMethod,
};
use crate::Kinetics::experimental_kinetics::math::fitting::{
    FitColumnRequest, FitColumnResult, FitColumnSession, FitEvaluationMode, FittedEvaluationBatch,
    FittingModelName, evaluate_fit_batch,
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
    pub fit_sessions: Vec<FitColumnSession>,
    pub fit_training_range_enabled: bool,
    pub fit_training_start: f64,
    pub fit_training_end: f64,
    pub fit_gap_enabled: bool,
    pub fit_gap_start: f64,
    pub fit_gap_end: f64,
    pub fit_evaluation_start: f64,
    pub fit_evaluation_end: f64,
    pub fit_evaluation_points: usize,
    pub fit_pending_batch: Option<FittedEvaluationBatch>,
    pub fit_reconstructed_experiment_id: String,
    pub sublimation_check_open: bool,
    pub sublimation_status: String,
    pub sublimation_ea_kj: Option<f64>,
    pub sublimation_k: Option<f64>,
    pub sublimation_verdict: Option<String>,
    pub sublimation_x_min: String,
    pub sublimation_x_max: String,
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
            fit_sessions: Vec::new(),
            fit_training_range_enabled: false,
            fit_training_start: 0.0,
            fit_training_end: 1.0,
            fit_gap_enabled: false,
            fit_gap_start: 0.0,
            fit_gap_end: 1.0,
            fit_evaluation_start: 0.0,
            fit_evaluation_end: 1.0,
            fit_evaluation_points: 20,
            fit_pending_batch: None,
            fit_reconstructed_experiment_id: String::new(),
            sublimation_check_open: false,
            sublimation_status: String::new(),
            sublimation_ea_kj: None,
            sublimation_k: None,
            sublimation_verdict: None,
            sublimation_x_min: String::new(),
            sublimation_x_max: String::new(),
        }
    }
}

impl KineticMethodsWindowState {
    pub(crate) fn close_window_on_layer(&mut self, layer: Option<egui::LayerId>) -> bool {
        let Some(layer) = layer else { return false };
        let windows = [
            ("Isoconversional", &mut self.isoconversional_open),
            ("Kissinger", &mut self.kissinger_open),
            (
                "Combined Kinetics Analysis",
                &mut self.combined_kinetics_open,
            ),
            ("Is this a sublimation?", &mut self.sublimation_check_open),
            ("Fit Model", &mut self.fit_model_open),
        ];
        for (title, is_open) in windows {
            if layer.id == egui::Id::unique(Some(title)) && *is_open {
                *is_open = false;
                return true;
            }
        }
        false
    }
}

pub struct KineticMethods;

impl KineticMethods {
    fn publish_status(model: &mut PlotModel, status: &str) {
        let _ = model.push_message(status);
    }

    /// Render the dropdown menu entries under the "Kinetic Methods" top menu.
    pub fn show_menu(ui: &mut egui::Ui, state: &mut KineticMethodsWindowState) {
        if ui.button("Isoconversional").on_hover_text("Non-isothermal activation-energy analysis over conversion; can create a result experiment.").clicked() {
            state.isoconversional_open = true;
        }
        if ui.button("Kissinger").on_hover_text("Non-isothermal peak-temperature method for one activation energy; updates results only.").clicked() {
            state.kissinger_open = true;
        }
        if ui.button("Combined Kinetics Analysis").on_hover_text("Non-isothermal model-parameter search; numeric search ranges are configured in its window.").clicked() {
            state.combined_kinetics_open = true;
        }
        ui.label("Criado Master Curve: unavailable pending scientific validation");
        if ui.button("Fit Model").on_hover_text("General single-experiment curve fitting, not a kinetic method; adds a fitted output column.").clicked() {
            state.fit_model_open = true;
        }
        if ui.button("Is this a sublimation?").on_hover_text("Diagnose one selected experiment or selected conversion interval using conversion rate and temperature; adds a fitted-rate column to that experiment.").clicked() {
            state.sublimation_check_open = true;
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
                    if state.isoconversional_method == IsoconversionalMethod::FriedmanIntegral {
                        ui.label("Isothermal Friedman method: uses the setpoint temperature from experiment metadata, so a measured T column is optional. Requires conversion data from at least three isothermal experiments.");
                    } else {
                        ui.label("Non-isothermal method: estimates activation energy versus conversion from experiments run at different heating rates. Requires a measured temperature column and heating-rate metadata; isothermal experiments are rejected.");
                    }
                    ui.add_space(4.0);
                    ui.horizontal(|ui| {
                        let label = ui.label("Method:");
                        let combo = egui::ComboBox::from_id_salt("isoconv_method")
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
                        let _ = combo.response.labelled_by(label.id);
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
                        let label = ui.label("Output ID:");
                        let input = ui.text_edit_singleline(&mut state.isoconversional_output_id);
                        let _ = input.labelled_by(label.id);
                    });

                    ui.add_space(6.0);

                    ui.horizontal(|ui| {
                        ui.label("Status:");
                        ui.label(&state.isoconversional_status);
                    });

                    ui.add_space(8.0);

                    // Run the chosen method and store the result as a new experiment.
                    let calculate_help = if state.isoconversional_method == IsoconversionalMethod::FriedmanIntegral {
                        "Requires at least three isothermal experiments with a finite setpoint temperature and conversion data. Adds a result experiment; source data are unchanged."
                    } else {
                        "Requires non-isothermal experiments with a measured temperature column and heating-rate metadata. Isothermal experiments are rejected. Adds a result experiment; source data are unchanged."
                    };
                    if ui.button("Calculate").on_hover_text(calculate_help).clicked() {
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
                            let warnings = result.warnings.clone();

                            model.push_isoconversional_result(&result, &output_id)?;

                            if warnings.is_empty() {
                                Ok(format!("Done. Result saved as '{}'.", output_id))
                            } else {
                                Ok(format!(
                                    "Result saved as '{}'; {} diagnostic warning(s): {}",
                                    output_id,
                                    warnings.len(),
                                    warnings.join(" | ")
                                ))
                            }
                        })();

                        state.isoconversional_status = match status {
                            Ok(msg) => msg,
                            Err(err) => format!("Error: {err}"),
                        };
                        Self::publish_status(model, &state.isoconversional_status);
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
                    ui.label("Non-isothermal Kissinger peak method: estimates one activation energy from peak temperatures measured at different heating rates. No numeric or text input is required; results appear below and source data is unchanged.");
                    ui.add_space(4.0);
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

                    // Run the Kissinger calculation without creating new series data or curves.
                    if ui.button("Calculate").on_hover_text("Requires multiple non-isothermal experiments with different heating rates. Updates Ea/R2 labels only; no column or experiment is created.").clicked() {
                        state.kissinger_ea_kj = None;
                        state.kissinger_r2 = None;
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
                            Err(err) => format!("Error: {err}"),
                        };
                        Self::publish_status(model, &state.kissinger_status);
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
                    ui.label("Diagnoses one selected experiment. Leave Conversion X bounds empty to use the whole experiment, enter both bounds to restrict the fit, or drag-select an X interval on a conversion-vs-rate curve. Entered bounds take precedence over a plot selection. A used selection touching multiple curves is rejected. The fit uses only the chosen X interval; fitted rates are evaluated across the full experiment and stored in a new column.");
                    ui.horizontal(|ui| {
                        ui.label("Conversion X min (optional):");
                        ui.add(egui::TextEdit::singleline(&mut state.sublimation_x_min)
                            .hint_text("min conversion; empty = full range"));
                        ui.label("max:");
                        ui.add(egui::TextEdit::singleline(&mut state.sublimation_x_max)
                            .hint_text("max conversion; empty = full range"));
                    });
                    ui.add_space(4.0);
                    let ea_label = state
                        .sublimation_ea_kj
                        .map(|v| format!("{:.2}", v))
                        .unwrap_or_else(|| "-".to_string());
                    let k_label = state
                        .sublimation_k
                        .map(|v| format!("{:.4e}", v))
                        .unwrap_or_else(|| "-".to_string());
                    let verdict_label = state
                        .sublimation_verdict
                        .clone()
                        .unwrap_or_else(|| "-".to_string());

                    ui.horizontal(|ui| {
                        ui.label("Ea (kJ/mol):");
                        ui.label(ea_label);
                    });

                    ui.horizontal(|ui| {
                        ui.label("Pre-exponential factor k:");
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

                    // Fit the selected experiment, then add its full-length predicted rate column.
                    if ui.button("Calculate").on_hover_text("Requires one selected experiment with conversion, conversion rate and temperature. Optional bounds or the selected X interval restrict the fit; fitted rates are added over the full experiment.").clicked() {
                        let status = (|| -> Result<String, TGAGUIError> {
                            let selected_index = model.get_selected_curve_index().ok_or_else(|| TGAGUIError::BindingError("Select one curve for the experiment to diagnose.".into()))?;
                            let curve = model.plots.get(selected_index).ok_or_else(|| TGAGUIError::BindingError("Selected curve is no longer available.".into()))?;
                            let exp_id = curve.experiment_id.clone();
                            let x_name = curve.x_name.clone();
                            let parse_bound = |text: &str| -> Result<Option<f64>, TGAGUIError> {
                                if text.trim().is_empty() { Ok(None) } else {
                                    let value = text.trim().parse::<f64>().map_err(|_| TGAGUIError::BindingError("Conversion X bounds must be numbers.".into()))?;
                                    if value.is_finite() { Ok(Some(value)) } else { Err(TGAGUIError::BindingError("Conversion X bounds must be finite.".into())) }
                                }
                            };
                            let lower = parse_bound(&state.sublimation_x_min)?;
                            let upper = parse_bound(&state.sublimation_x_max)?;
                            if lower.is_some() != upper.is_some() {
                                return Err(TGAGUIError::BindingError("Enter both conversion X bounds, or leave both empty.".into()));
                            }
                                    let selected_range = if lower.is_none() {
                                if let Some(rect) = model.interaction.selection_rect {
                                    let (x_min, x_max, _, _) = rect.bounds();
                                    let intersecting = model.plots.iter().filter(|plot| plot.highlighted).count();
                                    if intersecting != 1 || !curve.highlighted {
                                        return Err(TGAGUIError::BindingError(format!("The selected region intersects {intersecting} plotted curves; narrow it so exactly the selected curve is included.")));
                                    }
                                    let experiment = model.series.get_experiment_by_id(&exp_id)?;
                                    if experiment.dataset.schema.eta.as_deref() != Some(x_name.as_str()) {
                                        return Err(TGAGUIError::BindingError("For a plot selection, the curve X column must be the experiment's bound conversion column. Otherwise enter conversion X bounds above.".into()));
                                    }
                                    Some((x_min, x_max))
                                } else {
                                    None
                                }
                            } else {
                                None
                            };
                            let cols = SublimationMethod::default().required_columns_by_nature();
                            let view = model.series.create_kinetic_data_view(Some(&[exp_id.as_str()]), cols)?;
                            let exp_data = view.experiments.first().ok_or_else(|| TGAGUIError::BindingError("Selected experiment has no usable kinetic data.".into()))?;
                            let (x_min, x_max) = if let (Some(a), Some(b)) = (lower, upper) {
                                (a, b)
                            } else if let Some((a, b)) = selected_range {
                                (a, b)
                            } else {
                                (
                                    exp_data.conversion.iter().copied().fold(f64::INFINITY, f64::min),
                                    exp_data.conversion.iter().copied().fold(f64::NEG_INFINITY, f64::max),
                                )
                            };
                            let method = SublimationMethod::new(x_min, x_max);
                            let results = method.compute(&view)?;
                            let result = results.first().ok_or_else(|| TGAGUIError::BindingError("The selected experiment produced no sublimation result.".into()))?;

                            state.sublimation_ea_kj = Some(result.ea / 1000.0);
                            state.sublimation_k = Some(result.k);

                            let verdict_text = match result.resume {
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
                            Err(TGAGUIError::PostCommitWarning(warning)) => {
                                format!("Warning: {warning}")
                            }
                            Err(err) => format!("Error: {err}"),
                        };
                        Self::publish_status(model, &state.sublimation_status);
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
                    ui.label("Non-isothermal combined analysis: searches model exponents and activation energy using conversion, conversion rate and temperature. Set the numeric search ranges below; calculation updates result labels without adding columns.");
                    ui.add_space(4.0);
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
                    if ui.button("Calculate").on_hover_text("Requires non-isothermal temperature, conversion and conversion-rate data. Uses the numeric search settings; updates result labels only, with no new column.").clicked() {
                        state.combined_result_n = None;
                        state.combined_result_m = None;
                        state.combined_result_ea_kj = None;
                        state.combined_result_r2 = None;
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
                            Err(err) => format!("Error: {err}"),
                        };
                        Self::publish_status(model, &state.combined_status);
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
            .max_height(520.0)
            .vscroll(true)
            .show(ui.ctx(), |ui| {
                ui.vertical(|ui| {
                    ui.label("General curve fitting, not a kinetic method. Fit adds a new fitted-Y column to this experiment and plots it; measured data remain unchanged. Evaluate separately reconstructs a gap or extrapolates onto a new X grid.");
                    ui.add_space(4.0);
                    let exp_ids = model.list_of_experiments();
                    if state.fit_model_experiment_id.is_empty()
                        || !exp_ids
                            .iter()
                            .any(|id| id == &state.fit_model_experiment_id)
                    {
                        state.fit_model_experiment_id = model
                            .get_experiment_by_selected_curve()
                            .ok()
                            .filter(|id| exp_ids.is_empty() || exp_ids.contains(id))
                            .or_else(|| exp_ids.first().cloned())
                            .unwrap_or_default();
                    }

                    ui.horizontal(|ui| {
                        let label = ui.label("Experiment:");
                        let combo = egui::ComboBox::from_id_salt("fit_model_experiment")
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
                        let _ = combo.response.labelled_by(label.id);
                    });

                    let columns = if state.fit_model_experiment_id.is_empty() {
                        Ok(Vec::new())
                    } else {
                        model.list_of_columns(&state.fit_model_experiment_id)
                    };
                    let selected_curve_axes = model
                        .get_experiment_by_selected_curve()
                        .ok()
                        .filter(|id| id == &state.fit_model_experiment_id)
                        .and_then(|_| Some((model.this_is_x().ok()?, model.this_is_y().ok()?)));

                    match columns {
                        Ok(columns) => {
                            if state.fit_model_x_col.is_empty()
                                || !columns.iter().any(|col| col == &state.fit_model_x_col)
                            {
                                state.fit_model_x_col = selected_curve_axes
                                    .as_ref()
                                    .map(|(x, _)| x)
                                    .filter(|x| columns.contains(x))
                                    .cloned()
                                    .or_else(|| {
                                        columns.iter().find(|c| c.as_str() == "time").cloned()
                                    })
                                    .or_else(|| columns.iter().min().cloned())
                                    .unwrap_or_default();
                            }
                            if state.fit_model_y_col.is_empty()
                                || !columns.iter().any(|col| col == &state.fit_model_y_col)
                            {
                                state.fit_model_y_col = selected_curve_axes
                                    .as_ref()
                                    .map(|(_, y)| y)
                                    .filter(|y| **y != state.fit_model_x_col && columns.contains(y))
                                    .cloned()
                                    .or_else(|| {
                                        columns
                                            .iter()
                                            .find(|c| {
                                                c.as_str() == "mass" && **c != state.fit_model_x_col
                                            })
                                            .cloned()
                                    })
                                    .or_else(|| {
                                        columns
                                            .iter()
                                            .filter(|c| c.as_str() != state.fit_model_x_col)
                                            .min()
                                            .cloned()
                                    })
                                    .or_else(|| columns.iter().min().cloned())
                                    .unwrap_or_default();
                            }

                            ui.horizontal(|ui| {
                                let label = ui.label("X:");
                                let combo = egui::ComboBox::from_id_salt("fit_model_x")
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
                                let _ = combo.response.labelled_by(label.id);

                                let label = ui.label("Y:");
                                let combo = egui::ComboBox::from_id_salt("fit_model_y")
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
                                let _ = combo.response.labelled_by(label.id);
                            });

                            ui.horizontal(|ui| {
                                let label = ui.label("Model:");
                                let combo = egui::ComboBox::from_id_salt("fit_model_kind")
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
                                let _ = combo.response.labelled_by(label.id);
                            });

                            if state.fit_model_model.polynomial_degree().is_some() {
                                ui.horizontal(|ui| {
                                    let label = ui.label("Polynomial degree:");
                                    let response = ui.add(
                                        egui::DragValue::new(
                                            &mut state.fit_model_polynomial_degree,
                                        )
                                        .speed(1.0)
                                        .range(1..=12),
                                    );
                                    let changed = response.changed();
                                    let _ = response.labelled_by(label.id);
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
                                .default_open(false)
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
                                let label = ui.label("Output column:");
                                let input =
                                    ui.text_edit_singleline(&mut state.fit_model_output_col);
                                let _ = input.labelled_by(label.id);
                            });

                            ui.horizontal(|ui| {
                                let label = ui.label("Tolerance:");
                                let input = ui.add(
                                    egui::DragValue::new(&mut state.fit_model_tolerance)
                                        .speed(1e-7)
                                        .range(1e-12..=1e-1),
                                );
                                let _ = input.labelled_by(label.id);
                                let label = ui.label("Max iter:");
                                let input = ui.add(
                                    egui::DragValue::new(&mut state.fit_model_max_iter)
                                        .speed(10.0)
                                        .range(1..=100_000),
                                );
                                let _ = input.labelled_by(label.id);
                            });

                            ui.horizontal(|ui| {
                                let label = ui.label("Initial guess:");
                                let input =
                                    ui.text_edit_singleline(&mut state.fit_model_initial_guess);
                                let _ = input.labelled_by(label.id);
                            });

                            ui.separator();
                            ui.strong("Fit data");
                            ui.horizontal(|ui| {
                                ui.checkbox(&mut state.fit_training_range_enabled, "Limit fit range");
                                if state.fit_training_range_enabled {
                                    let label = ui.label("X start:");
                                    let response = ui.add(egui::DragValue::new(&mut state.fit_training_start).speed(0.1));
                                    let _ = response.labelled_by(label.id);
                                    let label = ui.label("end:");
                                    let response = ui.add(egui::DragValue::new(&mut state.fit_training_end).speed(0.1));
                                    let _ = response.labelled_by(label.id);
                                }
                            });
                            ui.horizontal(|ui| {
                                ui.checkbox(&mut state.fit_gap_enabled, "Exclude missing interval from fit");
                                if state.fit_gap_enabled {
                                    let label = ui.label("Gap start:");
                                    let response = ui.add(egui::DragValue::new(&mut state.fit_gap_start).speed(0.1));
                                    let _ = response.labelled_by(label.id);
                                    let label = ui.label("end:");
                                    let response = ui.add(egui::DragValue::new(&mut state.fit_gap_end).speed(0.1));
                                    let _ = response.labelled_by(label.id);
                                }
                            });

                            ui.horizontal(|ui| {
                                ui.label("Status:");
                                ui.label(&state.fit_model_status);
                            });

                            ui.add_space(8.0);
                            if ui.button("Fit").on_hover_text("Requires an experiment and X/Y columns. Polynomial degree is required for polynomial models; initial guess is optional. Adds a new fitted column to the same experiment, plots it, and keeps the model available for Evaluate. Existing columns are not overwritten.").clicked() {
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

                                    if state.fit_training_range_enabled {
                                        request = request.with_training_range(
                                            state.fit_training_start,
                                            state.fit_training_end,
                                        );
                                    }

                                    if state.fit_gap_enabled {
                                        request = request.excluding_range(
                                            state.fit_gap_start,
                                            state.fit_gap_end,
                                        );
                                    }

                                    let output = state.fit_model_output_col.trim();
                                    if !output.is_empty() {
                                        request = request.with_output_col(output.to_string());
                                    }

                                    if let Some(initial_guess) =
                                        parse_fit_initial_guess(&state.fit_model_initial_guess)?
                                    {
                                        request = request.with_initial_guess(initial_guess);
                                    }

                                    let session = model.fit_session_and_plot(request)?;
                                    let result = session.result.clone();
                                    state.fit_model_result = Some(result.clone());
                                    state.fit_sessions.retain(|existing| {
                                        existing.result.experiment_id != result.experiment_id
                                            || existing.result.y_col != result.y_col
                                            || existing.result.x_col != result.x_col
                                    });
                                    state.fit_sessions.push(session);
                                    state.fit_pending_batch = None;
                                    Ok(format!("Fit added '{}' and plotted. R2={:.6}; measured data unchanged.", result.output_col, result.r2))
                                })();

                                state.fit_model_status = match status {
                                    Ok(msg) => msg,
                                    Err(err) => format!("Error: {err}"),
                                };
                                Self::publish_status(model, &state.fit_model_status);
                            }

                            if !state.fit_sessions.is_empty() {
                                ui.separator();
                                ui.strong("Stored fits");
                                egui::Grid::new("fit_model_session_results")
                                    .striped(true)
                                    .show(ui, |ui| {
                                        ui.strong("Y");
                                        ui.strong("Experiment");
                                        ui.strong("Model");
                                        ui.strong("R2");
                                        ui.strong("Coefficients");
                                        ui.strong("Side RMSE");
                                        ui.end_row();
                                        for session in &state.fit_sessions {
                                            let result = &session.result;
                                            let coefficients = result.coefficients.iter()
                                                .map(|(name, value)| format!("{name}={value:.5e}"))
                                                .collect::<Vec<_>>().join(", ");
                                            let residuals = match (session.left_rmse, session.right_rmse) {
                                                (Some(left), Some(right)) => format!("{left:.3e} / {right:.3e}"),
                                                _ => "-".to_string(),
                                            };
                                            ui.label(&result.y_col);
                                            ui.label(&result.experiment_id);
                                            ui.label(result.kinetic_model.display_name());
                                            ui.label(format!("{:.6}", result.r2));
                                            ui.monospace(coefficients);
                                            ui.label(residuals);
                                            ui.end_row();
                                        }
                                    });
                            }

                            ui.separator();
                            ui.strong("Evaluate on a new X grid");
                            ui.label("This is optional reconstruction, not required to inspect the fitted curve.");
                            ui.horizontal(|ui| {
                                if state.fit_gap_enabled {
                                    ui.label(format!("Fill excluded interval [{:.4}, {:.4}]", state.fit_gap_start, state.fit_gap_end));
                                } else {
                                    let label = ui.label("Extrapolation X start:");
                                    let response = ui.add(egui::DragValue::new(&mut state.fit_evaluation_start).speed(0.1));
                                    let _ = response.labelled_by(label.id);
                                    let label = ui.label("end:");
                                    let response = ui.add(egui::DragValue::new(&mut state.fit_evaluation_end).speed(0.1));
                                    let _ = response.labelled_by(label.id);
                                }
                                let label = ui.label("Points:");
                                let response = ui.add(egui::DragValue::new(&mut state.fit_evaluation_points).range(1..=100_000));
                                let _ = response.labelled_by(label.id);
                            });
                            ui.horizontal(|ui| {
                                if ui.button("Evaluate fitted columns")
                                    .on_hover_text("Evaluates the retained fit models on a shared new X grid for extrapolation or gap reconstruction. This does not replace the fitted column; review the preview and explicitly create a reconstructed experiment.")
                                    .clicked()
                                {
                                    // A failed evaluation must not leave an older preview available
                                    // for accidental commit from Column Manager.
                                    state.fit_pending_batch = None;
                                    let matching_sessions: Vec<_> = state
                                        .fit_sessions
                                        .iter()
                                        .filter(|session| {
                                            session.result.experiment_id
                                                == state.fit_model_experiment_id
                                                && session.result.x_col == state.fit_model_x_col
                                        })
                                        .cloned()
                                        .collect();
                                    let mode = if state.fit_gap_enabled {
                                        FitEvaluationMode::FillGap {
                                            start: state.fit_gap_start,
                                            end: state.fit_gap_end,
                                            missing_points: state.fit_evaluation_points,
                                        }
                                    } else {
                                        FitEvaluationMode::Extrapolate {
                                            start: state.fit_evaluation_start,
                                            end: state.fit_evaluation_end,
                                            points: state.fit_evaluation_points,
                                        }
                                    };
                                    if matching_sessions.is_empty() {
                                        state.fit_model_status = format!(
                                            "Evaluation failed: no fitted model for experiment '{}' and X column '{}'.",
                                            state.fit_model_experiment_id,
                                            state.fit_model_x_col
                                        );
                                    } else {
                                        match evaluate_fit_batch(&matching_sessions, mode) {
                                            Ok(batch) => {
                                                state.fit_reconstructed_experiment_id = format!(
                                                    "{}_reconstructed",
                                                    batch.source_experiment_id
                                                );
                                                state.fit_pending_batch = Some(batch);
                                                state.fit_model_status = "Evaluation ready. Review and commit it in Column Manager.".to_string();
                                            }
                                            Err(error) => {
                                                state.fit_model_status =
                                                    format!("Evaluation failed: {error}")
                                            }
                                        }
                                    }
                                }
                                if ui.button("Clear fit sessions").clicked() {
                                    state.fit_sessions.clear();
                                    state.fit_model_result = None;
                                    state.fit_pending_batch = None;
                                }
                                ui.label(format!("{} fitted column(s)", state.fit_sessions.len()));
                            });
                        }
                        Err(err) => {
                            ui.label(format!("Cannot read columns: {err}"));
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
