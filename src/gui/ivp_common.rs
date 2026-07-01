use crate::Kinetics::experimental_kinetics::one_experiment_dataset::{
    ColumnMeta, ColumnNature, ColumnOrigin, History, TGADataset, TGADomainError, TGASchema, Unit,
};
use crate::Kinetics::solid_state_kinetics_IVP::{KineticModelIVP, KineticModelNames};
use RustedSciThe::numerical::ODE_api2::SolverType;
use RustedSciThe::numerical::Radau::Radau_main::RadauOrder;
use egui_extras::{Column as TableColumn, TableBuilder};
use nalgebra::{DMatrix, DVector};
use polars::prelude::{Column as PolarsColumn, DataFrame, IntoLazy};
use std::collections::HashMap;
use strum::IntoEnumIterator;
use strum_macros::EnumIter;

#[derive(Debug, Clone, PartialEq, EnumIter)]
pub enum SolverTypeShotcut {
    RadauOrder3,
    RadauOrder5,
    RadauOrder7,
    LSODE,
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
            Self::LSODE => SolverType::LSODE2,
            Self::NonstiffRK45 => SolverType::NonStiff("RK45".to_string()),
            Self::NonStiffDoPri => SolverType::NonStiff("DOPRI".to_string()),
            Self::NonStiffAB4 => SolverType::NonStiff("AB4".to_string()),
            Self::BDF => SolverType::BDF,
            Self::BackwardEuler => SolverType::BackwardEuler,
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
        TableBuilder::new(ui)
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
    }

    pub fn show_solver_combo(&mut self, ui: &mut egui::Ui) {
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
    }

    pub fn show_parameter_inputs(&mut self, ui: &mut egui::Ui) {
        self.sync_param_buffers();
        ui.horizontal_wrapped(|ui| {
            for (idx, param) in self.required_params.iter().enumerate() {
                ui.label(param);
                if ui
                    .text_edit_singleline(&mut self.param_strings[idx])
                    .changed()
                {
                    if let Ok(parsed) = self.param_strings[idx].trim().parse::<f64>() {
                        self.params[idx] = parsed;
                    }
                }
                ui.label(format!("{}", self.params[idx]));
            }
        });
    }

    pub fn show_problem_inputs(&mut self, ui: &mut egui::Ui) {
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

        let solver = self.selected_solver.get_solvertype();
        let mut kinetic_model = KineticModelIVP::new(solver);
        kinetic_model.set_problem(self.t_final, self.beta, self.t0, self.e, self.a)?;
        kinetic_model.set_model(model, params)?;
        kinetic_model.check_task()?;
        kinetic_model.solve()?;
        let (time, values) = kinetic_model.get_result()?;

        Ok(IvpResult { time, values })
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
