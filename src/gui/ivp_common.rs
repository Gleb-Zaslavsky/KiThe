use crate::Kinetics::experimental_kinetics::one_experiment_dataset::{
    ColumnMeta, ColumnNature, History, TGADataset, TGADomainError, TGASchema, Unit,
};
use crate::Kinetics::solid_state_kinetics_IVP::{KineticModelIVP, KineticModelNames};
use RustedSciThe::numerical::StiffAPI::StiffSolverKind;
use egui_extras::{Column as TableColumn, TableBuilder};
use nalgebra::{DMatrix, DVector};
use polars::prelude::{Column as PolarsColumn, DataFrame, IntoLazy};
use std::collections::HashMap;
use strum::IntoEnumIterator;
use strum_macros::EnumIter;

#[derive(Debug, Clone, PartialEq, EnumIter)]
pub enum SolverTypeShotcut {
    Radau,
    LSODE2,
    BDF,
    BackwardEuler,
}

impl SolverTypeShotcut {
    pub fn stiff_solver_kind(&self) -> StiffSolverKind {
        match self {
            Self::Radau => StiffSolverKind::Radau,
            Self::LSODE2 => StiffSolverKind::Lsode2,
            Self::BDF => StiffSolverKind::Bdf,
            Self::BackwardEuler => StiffSolverKind::BackwardEuler,
        }
    }

    fn label(&self) -> &'static str {
        match self {
            Self::Radau => "Radau",
            Self::LSODE2 => "LSODE2",
            Self::BDF => "BDF",
            Self::BackwardEuler => "Backward Euler",
        }
    }
}

#[derive(Debug, Clone)]
pub struct IvpResult {
    pub time: DVector<f64>,
    pub values: DMatrix<f64>,
}

impl IvpResult {
    pub fn primary_series_name(&self) -> &'static str {
        "conversion"
    }
}

#[derive(Debug, Clone)]
pub struct IvpTaskState {
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
    param_strings: Vec<String>,
    params: Vec<f64>,
}

impl Default for IvpTaskState {
    fn default() -> Self {
        Self::new()
    }
}

impl IvpTaskState {
    pub fn new() -> Self {
        Self {
            selected_solver: SolverTypeShotcut::BDF,
            selected_model: None,
            required_params: Vec::new(),
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
            param_strings: Vec::new(),
            params: Vec::new(),
        }
    }

    pub fn select_model_clicked(&mut self, model: KineticModelNames) {
        self.required_params = model.required_params();
        self.param_strings = vec![String::new(); self.required_params.len()];
        self.params = vec![0.0; self.required_params.len()];
        self.selected_model = Some(model);
    }

    fn sync_param_buffers(&mut self) {
        let required = self.required_params.len();
        if self.param_strings.len() != required {
            self.param_strings.resize(required, String::new());
        }
        if self.params.len() != required {
            self.params.resize(required, 0.0);
        }
    }

    pub fn show_model_table(&mut self, ui: &mut egui::Ui) {
        let catalog_height = (ui.available_height() - 300.0).clamp(220.0, 560.0);
        TableBuilder::new(ui)
            .vscroll(true)
            .min_scrolled_height(180.0)
            .max_scroll_height(catalog_height)
            .column(TableColumn::auto().resizable(true))
            .column(TableColumn::auto().resizable(true))
            .column(TableColumn::remainder())
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
                    body.row(24.0, |mut row| {
                        row.col(|ui| {
                            if ui.button("Use model").clicked() {
                                self.select_model_clicked(model.clone());
                            }
                        });
                        row.col(|ui| {
                            ui.label(format!("{:?}", model));
                        });
                        row.col(|ui| {
                            let formula = model.map_of_names_and_formulas();
                            ui.label(&formula).on_hover_text(formula);
                        });
                    });
                }
            });
    }

    pub fn show_solver_combo(&mut self, ui: &mut egui::Ui) {
        egui::ComboBox::from_label("Select Solver")
            .selected_text(self.selected_solver.label())
            .show_ui(ui, |ui| {
                for solver in SolverTypeShotcut::iter() {
                    ui.selectable_value(&mut self.selected_solver, solver.clone(), solver.label());
                }
            });
    }

    pub fn show_parameter_inputs(&mut self, ui: &mut egui::Ui) {
        self.sync_param_buffers();
        ui.horizontal_wrapped(|ui| {
            for (idx, param) in self.required_params.iter().enumerate() {
                let label = ui.label(param);
                ui.label_inputs_by(label.id, |ui| {
                    if ui
                        .text_edit_singleline(&mut self.param_strings[idx])
                        .changed()
                    {
                        if let Ok(parsed) = self.param_strings[idx].trim().parse::<f64>() {
                            self.params[idx] = parsed;
                        }
                    }
                });
                ui.label(format!("{}", self.params[idx]));
            }
        });
    }

    pub fn show_problem_inputs(&mut self, ui: &mut egui::Ui) {
        ui.horizontal_wrapped(|ui| {
            let input = |ui: &mut egui::Ui,
                         label_text: &str,
                         value: &mut String,
                         parsed_value: &mut f64| {
                ui.horizontal(|ui| {
                    let label = ui.label(format!("{label_text}:"));
                    ui.label_inputs_by(label.id, |ui| {
                        ui.add_sized(
                            [100.0, ui.spacing().interact_size.y],
                            egui::TextEdit::singleline(value),
                        );
                    });
                });
                if let Ok(parsed) = value.trim().parse::<f64>() {
                    *parsed_value = parsed;
                }
            };

            input(
                ui,
                "t_final [s]",
                &mut self.t_final_string,
                &mut self.t_final,
            );
            input(ui, "beta [K/s]", &mut self.beta_string, &mut self.beta);
            input(ui, "T0 [K]", &mut self.t0_string, &mut self.t0);
            input(ui, "E [J/mol]", &mut self.e_string, &mut self.e);
            input(ui, "A [1/s]", &mut self.a_string, &mut self.a);
        });
    }

    pub fn show_summary(&self, ui: &mut egui::Ui) {
        let text = format!(
            "Model: {:?}, A={}, E={}, beta={}, params={:?}, solver {:?}",
            self.selected_model, self.a, self.e, self.beta, self.params, self.selected_solver
        );
        ui.label(text);
    }

    pub fn run_kinetic_model(&mut self) -> Result<IvpResult, String> {
        let model = self
            .selected_model
            .clone()
            .ok_or_else(|| "No model selected".to_string())?;
        self.sync_param_buffers();

        let params = self
            .param_strings
            .iter()
            .enumerate()
            .map(|(idx, value)| {
                value.trim().parse::<f64>().map_err(|_| {
                    format!(
                        "Failed to parse parameter {} for {:?}: '{}'",
                        idx + 1,
                        model,
                        value
                    )
                })
            })
            .collect::<Result<Vec<_>, _>>()?;
        self.params = params.clone();

        let solver = self.selected_solver.stiff_solver_kind();
        let mut kinetic_model = KineticModelIVP::new(solver);
        kinetic_model.set_problem(self.t_final, self.beta, self.t0, self.e, self.a)?;
        kinetic_model.set_model(model, params)?;
        kinetic_model.check_task()?;
        kinetic_model.solve()?;
        let (time, values) = kinetic_model.get_result()?;

        Ok(IvpResult { time, values })
    }
}

#[cfg(test)]
mod solver_selector_tests {
    use super::*;

    #[test]
    fn selector_exposes_each_stiff_api_route_once() {
        let choices: Vec<_> = SolverTypeShotcut::iter().collect();

        assert_eq!(
            choices,
            [
                SolverTypeShotcut::Radau,
                SolverTypeShotcut::LSODE2,
                SolverTypeShotcut::BDF,
                SolverTypeShotcut::BackwardEuler,
            ]
        );
        assert_eq!(choices[0].label(), "Radau");
        assert_eq!(choices[1].label(), "LSODE2");
        assert_eq!(choices[2].label(), "BDF");
        assert_eq!(choices[3].label(), "Backward Euler");

        assert!(matches!(
            choices[0].stiff_solver_kind(),
            StiffSolverKind::Radau
        ));
        assert!(matches!(
            choices[1].stiff_solver_kind(),
            StiffSolverKind::Lsode2
        ));
        assert!(matches!(
            choices[2].stiff_solver_kind(),
            StiffSolverKind::Bdf
        ));
        assert!(matches!(
            choices[3].stiff_solver_kind(),
            StiffSolverKind::BackwardEuler
        ));
    }
}

pub fn ivp_result_to_experiment(
    result: &IvpResult,
    experiment_id: String,
) -> Result<
    crate::Kinetics::experimental_kinetics::experiment_series_main::TGAExperiment,
    TGADomainError,
> {
    let time_values: Vec<f64> = result.time.iter().copied().collect();
    let mut columns: Vec<PolarsColumn> =
        vec![PolarsColumn::new("time".into(), time_values.as_slice())];
    let mut schema_columns = HashMap::new();
    schema_columns.insert(
        "time".to_string(),
        ColumnMeta::raw("time".to_string(), Unit::Second, ColumnNature::Time),
    );

    for col_idx in 0..result.values.ncols() {
        let values: Vec<f64> = result.values.column(col_idx).iter().copied().collect();
        let (name, nature) = if col_idx == 0 {
            (
                result.primary_series_name().to_string(),
                ColumnNature::Conversion,
            )
        } else {
            (format!("state_{}", col_idx + 1), ColumnNature::Unknown)
        };
        columns.push(PolarsColumn::new(name.clone().into(), values.as_slice()));
        schema_columns.insert(
            name.clone(),
            ColumnMeta::raw(name, Unit::Dimensionless, nature),
        );
    }

    let frame = DataFrame::new(time_values.len(), columns)?.lazy();
    let mut dataset = TGADataset {
        frame,
        schema: TGASchema {
            columns: schema_columns,
            time: Some("time".to_string()),
            temperature: None,
            mass: None,
            alpha: None,
            dm_dt: None,
            eta: None,
            deta_dt: None,
            dalpha_dt: None,
            dT_dt: None,
            E: None,
            R2: None,
        },
        oneframeplot: None,
        history_of_operations: History::new(),
        undo_stack: Vec::new(),
        undo_snapshot_latch: false,
    };
    dataset.initialize_column_provenance();

    Ok(
        crate::Kinetics::experimental_kinetics::experiment_series_main::TGAExperiment::new(dataset)
            .with_id(experiment_id),
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::{DMatrix, DVector};

    #[test]
    fn ivp_result_to_experiment_creates_conversion_column() {
        let result = IvpResult {
            time: DVector::from_vec(vec![0.0, 1.0, 2.0]),
            values: DMatrix::from_vec(3, 1, vec![0.0, 0.4, 1.0]),
        };

        let experiment = ivp_result_to_experiment(&result, "demo".to_string()).unwrap();
        assert_eq!(experiment.meta.id, "demo");
        assert!(experiment.dataset.get_column("time").is_ok());
        assert!(experiment.dataset.get_column("conversion").is_ok());
    }
}
