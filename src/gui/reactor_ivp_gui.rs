//! Condensed reactor IVP GUI.
//!
//! This screen is the reactor-IVP counterpart to the BVP combustion editor.
//! It keeps the user-facing solver surface typed and explicit:
//! - physics and initial state are edited separately from the solver backend;
//! - Lambdify and AOT are mutually exclusive execution routes;
//! - Sparse and Banded are the two matrix-storage choices shown to the user;
//! - postprocessing controls choose how solved results are exported;
//! - preview prints a tabulated snapshot plus the pretty-printed symbolic RHS.

use crate::Kinetics::mechfinder_api::ReactionData;

use crate::ReactorsBVP::reactor_BVP_utils::ScalingConfig;
use crate::ReactorsIVP::SimpleReactorIVP::{
    IvpError, IvpSolveSnapshot, ReactorIvpInitialStateConfig, ReactorIvpPhysicalConfig,
    SimpleReactorTask,
};
use crate::ReactorsIVP::solver_backend::{
    ReactorIvpExecutionBackend, ReactorIvpMatrixBackend, ReactorIvpMethod, ReactorIvpSolverConfig,
    ReactorIvpStopCondition, ReactorIvpSymbolicBackend,
};
use crate::gui::document_lifecycle::DocumentLifecycleState;
use crate::gui::gui_plot::PlotWindow;
use RustedSciThe::Utils::postprocessing::{PostprocessDataset, PostprocessPlan, PostprocessStatus};
use egui::{ComboBox, DragValue, Grid, TextEdit};
use log::{info, warn};

use RustedSciThe::numerical::LSODE2::{
    Lsode2AotProfile, Lsode2AotToolchain, Lsode2NativeExecutionConfig,
};
use std::collections::HashMap;
use std::collections::hash_map::DefaultHasher;
use std::hash::{Hash, Hasher};
use std::path::PathBuf;
use std::sync::mpsc;
use std::thread;
#[cfg(test)]
use std::time::Duration;
use tabled::{Table, Tabled};

#[path = "reactor_ivp_lifecycle.rs"]
mod reactor_ivp_lifecycle;

/// UI-only execution selector that keeps Lambdify and AOT mutually exclusive.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum ReactorIvpGuiExecutionMode {
    Lambdify,
    Aot,
}

/// Typed export controls for the reactor-IVP postprocessing section.
#[derive(Clone, Debug, PartialEq)]
pub struct IvpPostprocessingConfig {
    pub save_txt: bool,
    pub txt_path: String,
    pub save_csv: bool,
    pub csv_path: String,
    pub write_report: bool,
    pub report_path: String,
    pub plotters_png: bool,
    pub plotters_dir: String,
    pub gnuplot_png: bool,
    pub gnuplot_dir: String,
    pub terminal_plot: bool,
}

impl Default for IvpPostprocessingConfig {
    fn default() -> Self {
        Self {
            save_txt: false,
            txt_path: String::new(),
            save_csv: false,
            csv_path: String::new(),
            write_report: false,
            report_path: String::new(),
            plotters_png: false,
            plotters_dir: String::new(),
            gnuplot_png: false,
            gnuplot_dir: String::new(),
            terminal_plot: false,
        }
    }
}

impl IvpPostprocessingConfig {
    /// Return true when at least one export action is configured.
    fn has_actions(&self) -> bool {
        self.save_txt
            || self.save_csv
            || self.write_report
            || self.plotters_png
            || self.gnuplot_png
            || self.terminal_plot
    }

    /// Build a declarative LSODE2 postprocessing plan from the typed editor state.
    fn to_plan(&self, default_stem: &str) -> Result<Option<PostprocessPlan>, IvpError> {
        let mut plan = PostprocessPlan::new();

        if self.save_txt {
            plan = plan.save_txt(self.resolve_file_path(default_stem, "txt"));
        }
        if self.save_csv {
            plan = plan.save_csv(self.resolve_file_path(default_stem, "csv"));
        }
        if self.write_report {
            plan = plan.write_report(self.resolve_file_path(default_stem, "md"));
        }
        if self.plotters_png {
            plan = plan.plotters_png(self.resolve_dir_path(default_stem, "plotters"));
        }
        if self.gnuplot_png {
            plan = plan.gnuplot_png(self.resolve_dir_path(default_stem, "gnuplot"));
        }
        if self.terminal_plot {
            plan = plan.terminal_plot();
        }

        if plan.actions.is_empty() {
            Ok(None)
        } else {
            Ok(Some(plan))
        }
    }

    /// Resolve a file path, falling back to a stable task-specific default.
    fn resolve_file_path(&self, default_stem: &str, extension: &str) -> PathBuf {
        let explicit = match extension {
            "txt" => &self.txt_path,
            "csv" => &self.csv_path,
            "md" => &self.report_path,
            _ => "",
        };
        if !explicit.trim().is_empty() {
            return PathBuf::from(explicit.trim());
        }
        PathBuf::from(format!("{default_stem}.{extension}"))
    }

    /// Resolve an output directory, falling back to a stable task-specific default.
    fn resolve_dir_path(&self, default_stem: &str, suffix: &str) -> PathBuf {
        let explicit = match suffix {
            "plotters" => &self.plotters_dir,
            "gnuplot" => &self.gnuplot_dir,
            _ => "",
        };
        if !explicit.trim().is_empty() {
            return PathBuf::from(explicit.trim());
        }
        PathBuf::from(format!("{default_stem}_{suffix}"))
    }

    /// Format the configured export surface as compact preview rows.
    fn preview_rows(&self) -> Vec<ReactorIvpPreviewRow> {
        vec![
            preview_row("postprocessing", "save_txt", self.save_txt.to_string()),
            preview_row(
                "postprocessing",
                "txt_path",
                self.display_path(&self.txt_path),
            ),
            preview_row("postprocessing", "save_csv", self.save_csv.to_string()),
            preview_row(
                "postprocessing",
                "csv_path",
                self.display_path(&self.csv_path),
            ),
            preview_row(
                "postprocessing",
                "write_report",
                self.write_report.to_string(),
            ),
            preview_row(
                "postprocessing",
                "report_path",
                self.display_path(&self.report_path),
            ),
            preview_row(
                "postprocessing",
                "plotters_png",
                self.plotters_png.to_string(),
            ),
            preview_row(
                "postprocessing",
                "plotters_dir",
                self.display_path(&self.plotters_dir),
            ),
            preview_row(
                "postprocessing",
                "gnuplot_png",
                self.gnuplot_png.to_string(),
            ),
            preview_row(
                "postprocessing",
                "gnuplot_dir",
                self.display_path(&self.gnuplot_dir),
            ),
            preview_row(
                "postprocessing",
                "terminal_plot",
                self.terminal_plot.to_string(),
            ),
        ]
    }

    fn display_path(&self, value: &str) -> String {
        if value.trim().is_empty() {
            "auto".to_string()
        } else {
            value.to_string()
        }
    }
}

/// Explicit worker lifecycle for the reactor-IVP calculation button.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ReactorIvpCalculationState {
    Idle,
    Running { run_id: u64 },
    Completed,
    Failed,
}

impl ReactorIvpCalculationState {
    fn active_run_id(self) -> Option<u64> {
        match self {
            Self::Running { run_id } => Some(run_id),
            Self::Idle | Self::Completed | Self::Failed => None,
        }
    }

    fn is_active(self) -> bool {
        self.active_run_id().is_some()
    }
}

/// Owned worker result that can cross the background-thread boundary.
#[derive(Debug)]
struct ReactorIvpWorkerMessage {
    run_id: u64,
    config_fingerprint: u64,
    outcome: Result<
        (
            SimpleReactorTask,
            Vec<String>,
            nalgebra::DVector<f64>,
            nalgebra::DMatrix<f64>,
            ReactorIvpPreviewSnapshot,
            String,
        ),
        IvpError,
    >,
}

/// Visible consent request for an AOT reactor-IVP run.
#[derive(Clone, Debug, PartialEq, Eq)]
struct ReactorIvpAotRequest {
    config_fingerprint: u64,
    summary: String,
}

impl ReactorIvpAotRequest {
    fn from_config(config: &IvpGuiConfig) -> Self {
        let execution_backend = match config.execution_mode {
            ReactorIvpGuiExecutionMode::Lambdify => "Lambdify".to_string(),
            ReactorIvpGuiExecutionMode::Aot => {
                format!("AOT ({:?}, {:?})", config.aot_toolchain, config.aot_profile)
            }
        };
        Self {
            config_fingerprint: ReactorIvpApp::fingerprint_config(config),
            summary: format!(
                "Execution backend: {execution_backend}; matrix backend: {:?}; symbolic backend: {:?}",
                config.matrix_backend, config.symbolic_backend
            ),
        }
    }
}

impl Default for ReactorIvpGuiExecutionMode {
    fn default() -> Self {
        Self::Lambdify
    }
}

/// Typed editor model for the condensed reactor IVP screen.
#[derive(Clone, Debug)]
pub struct IvpGuiConfig {
    pub problem_name: String,
    pub problem_description: String,
    pub physical: ReactorIvpPhysicalConfig,
    pub scaling: ScalingConfig,
    pub initial_state: ReactorIvpInitialStateConfig,
    pub thermal_effects: Vec<f64>,
    pub stop_condition_enabled: bool,
    pub stop_condition_species_index: usize,
    pub stop_condition_threshold: f64,
    pub method: ReactorIvpMethod,
    pub symbolic_backend: ReactorIvpSymbolicBackend,
    pub matrix_backend: ReactorIvpMatrixBackend,
    pub execution_mode: ReactorIvpGuiExecutionMode,
    pub aot_toolchain: Lsode2AotToolchain,
    pub aot_profile: Lsode2AotProfile,
    pub x0: f64,
    pub x_bound: f64,
    pub first_step: Option<f64>,
    pub max_step: f64,
    pub rtol: f64,
    pub atol: f64,
    pub postprocessing: IvpPostprocessingConfig,
    pub native_execution: Lsode2NativeExecutionConfig,
}

impl Default for IvpGuiConfig {
    fn default() -> Self {
        Self {
            problem_name: String::new(),
            problem_description: String::new(),
            physical: ReactorIvpPhysicalConfig {
                ro: 1200.0,
                Cp: 1000.0,
                Lambda: 0.25,
                m: 0.015,
                L: 0.05,
            },
            scaling: ScalingConfig::new(100.0, 0.05, 100.0),
            initial_state: ReactorIvpInitialStateConfig {
                temperature: 450.0,
                heat_flux: 0.15,
                species: HashMap::new(),
            },
            thermal_effects: vec![-2.0e5],
            stop_condition_enabled: false,
            stop_condition_species_index: 0,
            stop_condition_threshold: 1e-4,
            method: ReactorIvpMethod::default(),
            symbolic_backend: ReactorIvpSymbolicBackend::default(),
            matrix_backend: ReactorIvpMatrixBackend::default(),
            execution_mode: ReactorIvpGuiExecutionMode::default(),
            aot_toolchain: Lsode2AotToolchain::default(),
            aot_profile: Lsode2AotProfile::default(),
            x0: 0.0,
            x_bound: 0.01,
            first_step: None,
            max_step: 1.0e-3,
            rtol: 1e-6,
            atol: 1e-8,
            postprocessing: IvpPostprocessingConfig::default(),
            // Interactive GUI defaults should favor the stable bridge path.
            native_execution: Lsode2NativeExecutionConfig::bridge_solve(),
        }
    }
}

impl IvpGuiConfig {
    /// Fold the typed GUI model into the solver facade config.
    pub fn to_solver_config(&self) -> ReactorIvpSolverConfig {
        let execution_backend = match self.execution_mode {
            ReactorIvpGuiExecutionMode::Lambdify => ReactorIvpExecutionBackend::Lambdify,
            ReactorIvpGuiExecutionMode::Aot => ReactorIvpExecutionBackend::Aot {
                toolchain: self.aot_toolchain,
                profile: self.aot_profile,
            },
        };

        ReactorIvpSolverConfig {
            method: self.method,
            symbolic_backend: self.symbolic_backend,
            matrix_backend: self.matrix_backend,
            execution_backend,
            stop_condition: self.stop_condition_enabled.then(|| {
                ReactorIvpStopCondition::new(
                    self.stop_condition_species_index,
                    self.stop_condition_threshold,
                )
            }),
            x0: self.x0,
            x_bound: self.x_bound,
            first_step: self.first_step,
            max_step: self.max_step,
            rtol: self.rtol,
            atol: self.atol,
            native_execution: self.native_execution,
        }
    }

    /// Update the config from a task snapshot.
    pub fn from_task(task: &SimpleReactorTask) -> Self {
        let mut config = Self::default();
        config.problem_name = task.problem_name.clone().unwrap_or_default();
        config.problem_description = task.problem_description.clone().unwrap_or_default();
        config.physical.ro = task.ro;
        config.physical.Cp = task.Cp;
        config.physical.Lambda = task.Lambda;
        config.physical.m = task.m;
        config.physical.L = task.L;
        config.scaling = task.scaling.clone();
        config.initial_state.temperature = *task.initial_conditions.get("T").unwrap_or(&450.0);
        config.initial_state.heat_flux = *task.initial_conditions.get("q").unwrap_or(&0.15);
        config.initial_state.species = task
            .kindata
            .substances
            .iter()
            .filter_map(|name| {
                task.initial_conditions
                    .get(name)
                    .map(|value| (name.clone(), *value))
            })
            .collect();
        config.thermal_effects = task.thermal_effects.clone();
        if let Some(stop_condition) = task.solver_backend_config.stop_condition {
            config.stop_condition_enabled = true;
            config.stop_condition_species_index = stop_condition.species_index;
            config.stop_condition_threshold = stop_condition.threshold;
        }
        config.method = task.solver_backend_config.method;
        config.symbolic_backend = task.solver_backend_config.symbolic_backend;
        config.matrix_backend = task.solver_backend_config.matrix_backend;
        config.execution_mode = match task.solver_backend_config.execution_backend {
            ReactorIvpExecutionBackend::Lambdify => ReactorIvpGuiExecutionMode::Lambdify,
            ReactorIvpExecutionBackend::Aot { toolchain, profile } => {
                config.aot_toolchain = toolchain;
                config.aot_profile = profile;
                ReactorIvpGuiExecutionMode::Aot
            }
        };
        config.x0 = task.solver_backend_config.x0;
        config.x_bound = task.solver_backend_config.x_bound;
        config.first_step = task.solver_backend_config.first_step;
        config.max_step = task.solver_backend_config.max_step;
        config.rtol = task.solver_backend_config.rtol;
        config.atol = task.solver_backend_config.atol;
        config.native_execution = task.solver_backend_config.native_execution;
        config
    }

    /// Enable the default combustion stop condition using the first available species.
    ///
    /// This is used by the built-in demo task only. Loaded documents keep their
    /// authored state until the user changes the toggle explicitly.
    pub fn enable_default_stop_condition_if_possible(&mut self, species_count: usize) {
        if species_count == 0 {
            return;
        }
        self.stop_condition_enabled = true;
        self.stop_condition_species_index = 0;
        self.stop_condition_threshold = 1e-4;
    }

    /// Apply the typed GUI config to a reactor task.
    pub fn apply_to_task(&self, task: &mut SimpleReactorTask) -> Result<(), IvpError> {
        if self.problem_name.trim().is_empty() {
            task.problem_name = None;
        } else {
            task.problem_name = Some(self.problem_name.clone());
        }

        if self.problem_description.trim().is_empty() {
            task.problem_description = None;
        } else {
            task.problem_description = Some(self.problem_description.clone());
        }

        task.set_density(self.physical.ro)?;
        task.set_transport_properties(self.physical.Lambda, self.physical.Cp);
        task.m = self.physical.m;
        task.L = self.physical.L;
        task.set_scaling(self.scaling.clone())
            .map_err(IvpError::from)?;

        let mut initial_conditions = HashMap::with_capacity(2 + task.kindata.substances.len());
        initial_conditions.insert("T".to_string(), self.initial_state.temperature);
        initial_conditions.insert("q".to_string(), self.initial_state.heat_flux);
        for substance in &task.kindata.substances {
            let value = self.initial_state.species.get(substance).ok_or_else(|| {
                IvpError::MissingData(format!(
                    "Missing initial condition for substance `{}`",
                    substance
                ))
            })?;
            initial_conditions.insert(substance.clone(), *value);
        }
        task.set_initial_conditions(initial_conditions);
        task.set_thermal_effects(self.thermal_effects.clone());
        task.set_solver_backend_config(self.to_solver_config());
        Ok(())
    }

    /// Build an LSODE2 postprocessing plan from the current editor settings.
    pub fn postprocessing_plan(
        &self,
        task: &SimpleReactorTask,
    ) -> Result<Option<PostprocessPlan>, IvpError> {
        let mut stem = self
            .problem_name
            .trim()
            .chars()
            .map(|ch| {
                if ch.is_ascii_alphanumeric() || ch == '-' || ch == '_' {
                    ch
                } else {
                    '_'
                }
            })
            .collect::<String>();
        stem = stem.trim_matches('_').to_string();
        if stem.is_empty() {
            stem = task
                .problem_name
                .as_deref()
                .unwrap_or("reactor_ivp")
                .chars()
                .map(|ch| {
                    if ch.is_ascii_alphanumeric() || ch == '-' || ch == '_' {
                        ch
                    } else {
                        '_'
                    }
                })
                .collect::<String>();
            stem = stem.trim_matches('_').to_string();
        }
        if stem.is_empty() {
            stem = "reactor_ivp".to_string();
        }
        self.postprocessing.to_plan(&stem)
    }
}

/// One row in the preview summary table.
#[derive(Clone, Debug, Tabled)]
pub struct ReactorIvpPreviewRow {
    pub section: String,
    pub field: String,
    pub value: String,
}

/// One row in the pretty-printed equation table.
#[derive(Clone, Debug, Tabled)]
pub struct ReactorIvpEquationRow {
    pub index: usize,
    pub variable: String,
    pub equation: String,
}

/// Structured preview snapshot shown before a real solve.
#[derive(Clone, Debug)]
pub struct ReactorIvpPreviewSnapshot {
    pub summary_rows: Vec<ReactorIvpPreviewRow>,
    pub equation_rows: Vec<ReactorIvpEquationRow>,
}

impl ReactorIvpPreviewSnapshot {
    /// Print the preview in a human-readable form for console logs.
    pub fn print_to_console(&self) {
        info!(
            "Reactor IVP preview:\n{}",
            Table::new(self.summary_rows.clone())
        );
        info!(
            "Reactor IVP equations:\n{}",
            Table::new(self.equation_rows.clone())
        );
    }
}

/// Build a read-only preview snapshot without mutating the source task.
pub fn build_reactor_ivp_task_preview_snapshot(
    task: &SimpleReactorTask,
    config: &IvpGuiConfig,
) -> Result<ReactorIvpPreviewSnapshot, IvpError> {
    let mut preview_task = task.clone();
    config.apply_to_task(&mut preview_task)?;
    preview_task.setup_condensed_ivp()?;
    preview_task.check_condensed_before_solution()?;

    let mut summary_rows = vec![
        preview_row(
            "problem",
            "name",
            preview_task.problem_name.clone().unwrap_or_default(),
        ),
        preview_row(
            "problem",
            "description",
            preview_task.problem_description.clone().unwrap_or_default(),
        ),
        preview_row("physics", "ro", preview_task.ro.to_string()),
        preview_row("physics", "Cp", preview_task.Cp.to_string()),
        preview_row("physics", "Lambda", preview_task.Lambda.to_string()),
        preview_row("physics", "m", preview_task.m.to_string()),
        preview_row("physics", "L", preview_task.L.to_string()),
        preview_row("scaling", "dT", preview_task.scaling.dT.to_string()),
        preview_row(
            "scaling",
            "T_scale",
            preview_task.scaling.T_scale.to_string(),
        ),
        preview_row(
            "initial_state",
            "temperature",
            preview_task
                .initial_conditions
                .get("T")
                .map(|value| value.to_string())
                .unwrap_or_else(|| "missing".to_string()),
        ),
        preview_row(
            "initial_state",
            "heat_flux",
            preview_task
                .initial_conditions
                .get("q")
                .map(|value| value.to_string())
                .unwrap_or_else(|| "missing".to_string()),
        ),
        preview_row(
            "initial_state",
            "species_count",
            preview_task.kindata.substances.len().to_string(),
        ),
        preview_row(
            "integration",
            "x0",
            preview_task.solver_backend_config.x0.to_string(),
        ),
        preview_row(
            "integration",
            "x_bound",
            preview_task.solver_backend_config.x_bound.to_string(),
        ),
        preview_row(
            "integration",
            "first_step",
            preview_task
                .solver_backend_config
                .first_step
                .map(|value| value.to_string())
                .unwrap_or_else(|| "None".to_string()),
        ),
        preview_row(
            "integration",
            "max_step",
            preview_task.solver_backend_config.max_step.to_string(),
        ),
        preview_row(
            "integration",
            "rtol",
            preview_task.solver_backend_config.rtol.to_string(),
        ),
        preview_row(
            "integration",
            "atol",
            preview_task.solver_backend_config.atol.to_string(),
        ),
        preview_row(
            "solver",
            "method",
            format!("{:?}", preview_task.solver_backend_config.method),
        ),
        preview_row(
            "solver",
            "symbolic_backend",
            format!("{:?}", preview_task.solver_backend_config.symbolic_backend),
        ),
        preview_row(
            "solver",
            "matrix_backend",
            format!("{:?}", preview_task.solver_backend_config.matrix_backend),
        ),
        preview_row(
            "solver",
            "execution_backend",
            match preview_task.solver_backend_config.execution_backend {
                ReactorIvpExecutionBackend::Lambdify => "Lambdify".to_string(),
                ReactorIvpExecutionBackend::Aot { toolchain, profile } => {
                    format!("AOT ({:?}, {:?})", toolchain, profile)
                }
            },
        ),
        preview_row(
            "solver",
            "stop_condition",
            match preview_task.solver_backend_config.stop_condition {
                Some(stop_condition) => preview_task
                    .kindata
                    .substances
                    .get(stop_condition.species_index)
                    .map(|species| format!("{species} <= {}", stop_condition.threshold))
                    .unwrap_or_else(|| {
                        format!(
                            "index {} <= {}",
                            stop_condition.species_index, stop_condition.threshold
                        )
                    }),
                None => "disabled".to_string(),
            },
        ),
        preview_row(
            "chemistry",
            "reaction_count",
            preview_task.kindata.vec_of_equations.len().to_string(),
        ),
        preview_row(
            "chemistry",
            "species_count",
            preview_task.kindata.substances.len().to_string(),
        ),
        preview_row(
            "chemistry",
            "thermal_effects",
            preview_task.thermal_effects.len().to_string(),
        ),
    ];
    summary_rows.extend(config.postprocessing.preview_rows());

    let equation_rows = preview_task
        .solver
        .unknowns
        .iter()
        .enumerate()
        .filter_map(|(index, variable)| {
            preview_task
                .solver
                .eq_system
                .get(index)
                .map(|equation| ReactorIvpEquationRow {
                    index,
                    variable: variable.clone(),
                    equation: equation.pretty_print(),
                })
        })
        .collect::<Vec<_>>();

    Ok(ReactorIvpPreviewSnapshot {
        summary_rows,
        equation_rows,
    })
}

fn preview_row(section: &str, field: &str, value: impl Into<String>) -> ReactorIvpPreviewRow {
    ReactorIvpPreviewRow {
        section: section.to_string(),
        field: field.to_string(),
        value: value.into(),
    }
}

fn run_reactor_ivp_postprocessing(
    task: &SimpleReactorTask,
    config: &IvpGuiConfig,
    variable_names: &[String],
    axis: &nalgebra::DVector<f64>,
    values: &nalgebra::DMatrix<f64>,
) -> Result<Option<String>, IvpError> {
    let Some(plan) = config.postprocessing_plan(task)? else {
        return Ok(None);
    };

    let dataset = PostprocessDataset::new(
        task.solver.arg_name.clone(),
        variable_names.to_vec(),
        axis.clone(),
        values.clone(),
    )
    .map_err(|err| IvpError::CalculationError(err.to_string()))?;
    let dataset = dataset
        .with_metadata(
            "problem_name",
            task.problem_name.clone().unwrap_or_default(),
        )
        .with_metadata(
            "solver_backend",
            format!("{:?}", task.solver_backend_config.execution_backend),
        );
    let report = plan
        .execute(&dataset)
        .map_err(|err| IvpError::CalculationError(err.to_string()))?;

    if report.entries.is_empty() {
        return Ok(None);
    }

    let done = report
        .entries
        .iter()
        .filter(|entry| matches!(entry.status, PostprocessStatus::Done))
        .count();
    let skipped = report
        .entries
        .iter()
        .filter(|entry| matches!(entry.status, PostprocessStatus::Skipped))
        .count();
    let summary = if skipped > 0 {
        format!("Postprocessing completed: {done} done, {skipped} skipped")
    } else {
        format!("Postprocessing completed: {done} actions")
    };
    Ok(Some(summary))
}

/// Normalize the user-facing run status after LSODE2 finishes.
///
/// The native faithful route reports `finished_native_faithful` together with a
/// native termination kind. For the GUI, a stop-condition hit should still read
/// like a stop condition to the user rather than a transport-layer detail.
fn summarize_reactor_ivp_status(
    result: &IvpSolveSnapshot,
    postprocessing_summary: Option<String>,
) -> String {
    let base_status = if result
        .summary
        .native_termination_kind
        .as_deref()
        .is_some_and(|kind| kind == "reached_stop_condition")
    {
        "stopped_by_condition".to_string()
    } else {
        result.summary.status.clone()
    };

    match postprocessing_summary {
        Some(summary) => format!("{base_status}; {summary}"),
        None => base_status,
    }
}

/// Reactor-IVP GUI application shell.
///
/// The app keeps a typed editor model and commits it into a task only on
/// preview or run, so failed validation does not poison the stored task.
pub struct ReactorIvpApp {
    pub task: SimpleReactorTask,
    pub config: IvpGuiConfig,
    pub current_file_path: Option<PathBuf>,
    pub(crate) document_lifecycle: DocumentLifecycleState,
    plot_window: Option<PlotWindow>,
    calculation_state: ReactorIvpCalculationState,
    calculation_receiver: Option<mpsc::Receiver<ReactorIvpWorkerMessage>>,
    pub(crate) pending_document_load: Option<PathBuf>,
    pub(crate) pending_window_close: bool,
    pending_aot_confirmation: Option<ReactorIvpAotRequest>,
    last_error: Option<String>,
    last_status: Option<String>,
    last_preview: Option<ReactorIvpPreviewSnapshot>,
    next_run_id: u64,
}

impl Default for ReactorIvpApp {
    fn default() -> Self {
        Self::new()
    }
}

impl ReactorIvpApp {
    /// Create a reactor-IVP GUI with a compact demo task.
    pub fn new() -> Self {
        match build_demo_task() {
            Ok(task) => {
                let mut config = IvpGuiConfig::from_task(&task);
                config.enable_default_stop_condition_if_possible(task.kindata.substances.len());
                let mut app = Self {
                    task,
                    config,
                    current_file_path: None,
                    document_lifecycle: DocumentLifecycleState::default(),
                    plot_window: None,
                    calculation_state: ReactorIvpCalculationState::Idle,
                    calculation_receiver: None,
                    pending_document_load: None,
                    pending_window_close: false,
                    pending_aot_confirmation: None,
                    last_error: None,
                    last_status: Some("Ready".to_string()),
                    last_preview: None,
                    next_run_id: 1,
                };
                app.document_lifecycle
                    .mark_clean_snapshot(app.current_document_fingerprint().unwrap_or_default());
                app
            }
            Err(err) => {
                warn!("failed to seed demo reactor IVP task: {err}");
                Self {
                    task: SimpleReactorTask::new(),
                    config: IvpGuiConfig::default(),
                    current_file_path: None,
                    document_lifecycle: DocumentLifecycleState::default(),
                    plot_window: None,
                    calculation_state: ReactorIvpCalculationState::Idle,
                    calculation_receiver: None,
                    pending_document_load: None,
                    pending_window_close: false,
                    pending_aot_confirmation: None,
                    last_error: Some(err.to_string()),
                    last_status: None,
                    last_preview: None,
                    next_run_id: 1,
                }
            }
        }
    }

    /// Render the reactor-IVP window.
    pub fn show(&mut self, ctx: &egui::Context, open: &mut bool) {
        self.poll_calculation_worker(ctx);
        let was_open = *open;
        let mut window_open = *open;
        egui::Window::new("Condensed reactor IVP")
            .open(&mut window_open)
            .default_size([1280.0, 920.0])
            .vscroll(true)
            .show(ctx, |ui| {
                ui.heading("Condensed reactor IVP");
                ui.separator();

                self.render_document_toolbar(ui);
                ui.separator();

                self.render_problem_header(ui);
                ui.separator();
                self.render_physics_section(ui);
                ui.separator();
                self.render_initial_state_section(ui);
                ui.separator();
                self.render_solver_backend_section(ui);
                ui.separator();
                self.render_integration_section(ui);
                ui.separator();
                self.render_postprocessing_section(ui);
                ui.separator();
                self.render_conservation_section(ui);
                ui.separator();
                self.render_solve_report_section(ui);
                ui.separator();
                self.render_actions(ui);
                ui.separator();

                if let Some(status) = &self.last_status {
                    ui.label(status);
                }
                if let Some(error) = &self.last_error {
                    ui.colored_label(egui::Color32::from_rgb(170, 30, 30), error);
                }
            });

        if was_open && !window_open && !self.request_window_close() {
            *open = true;
        } else {
            *open = window_open;
        }
        self.show_window_close_confirmation(ctx, open);
        self.show_document_load_confirmation(ctx);
        self.show_aot_confirmation(ctx);

        if let Some(plot) = &mut self.plot_window {
            plot.show(ctx);
        }
    }

    /// Return the last visible status message for tests and higher-level GUI assertions.
    #[allow(dead_code)]
    pub(crate) fn last_status(&self) -> Option<&str> {
        self.last_status.as_deref()
    }

    /// Return the last visible error message for tests and higher-level GUI assertions.
    #[allow(dead_code)]
    pub(crate) fn last_error(&self) -> Option<&str> {
        self.last_error.as_deref()
    }

    /// Return the most recent preview snapshot, if any.
    #[allow(dead_code)]
    pub(crate) fn last_preview(&self) -> Option<&ReactorIvpPreviewSnapshot> {
        self.last_preview.as_ref()
    }

    /// Return the underlying task snapshot for read-only regression tests.
    #[allow(dead_code)]
    pub(crate) fn task(&self) -> &SimpleReactorTask {
        &self.task
    }

    /// Return whether the calculation worker is currently active.
    #[allow(dead_code)]
    pub(crate) fn calculation_state(&self) -> ReactorIvpCalculationState {
        self.calculation_state
    }

    /// Replace the active reactor task and resynchronize the typed editor state.
    ///
    /// This is the safe problem-switching primitive for future load paths and
    /// programmatic task replacement. It clears stale result surfaces so the new
    /// task starts from a clean editor/result boundary.
    #[allow(dead_code)]
    pub(crate) fn replace_task(&mut self, task: SimpleReactorTask) -> Result<(), IvpError> {
        if self.calculation_state.is_active() {
            return Err(IvpError::InvalidConfiguration(
                "cannot replace the active reactor IVP task while a calculation is running"
                    .to_string(),
            ));
        }

        self.config = IvpGuiConfig::from_task(&task);
        self.task = task;
        self.set_document_path(None);
        self.document_lifecycle
            .mark_clean_snapshot(self.current_document_fingerprint().unwrap_or_default());
        self.plot_window = None;
        self.last_preview = None;
        self.last_error = None;
        self.last_status = Some("Task replaced.".to_string());
        self.pending_document_load = None;
        self.pending_window_close = false;
        Ok(())
    }

    fn render_problem_header(&mut self, ui: &mut egui::Ui) {
        ui.horizontal_wrapped(|ui| {
            ui.label("Problem name");
            ui.add(TextEdit::singleline(&mut self.config.problem_name).desired_width(220.0));
            ui.label("Description");
            ui.add(TextEdit::singleline(&mut self.config.problem_description).desired_width(520.0));
        });
    }

    fn render_physics_section(&mut self, ui: &mut egui::Ui) {
        ui.heading("Physics");
        Grid::new("ivp_physics_grid")
            .num_columns(4)
            .spacing([16.0, 8.0])
            .show(ui, |ui| {
                ui.label("Density ro");
                ui.add(DragValue::new(&mut self.config.physical.ro).speed(1.0));
                ui.label("Heat capacity Cp");
                ui.add(DragValue::new(&mut self.config.physical.Cp).speed(1.0));
                ui.end_row();

                ui.label("Conductivity Lambda");
                ui.add(DragValue::new(&mut self.config.physical.Lambda).speed(0.01));
                ui.label("Mass flux m");
                ui.add(DragValue::new(&mut self.config.physical.m).speed(0.001));
                ui.end_row();

                ui.label("Length L");
                ui.add(DragValue::new(&mut self.config.physical.L).speed(0.001));
                ui.label("dT");
                ui.add(DragValue::new(&mut self.config.scaling.dT).speed(1.0));
                ui.end_row();

                ui.label("T_scale");
                ui.add(DragValue::new(&mut self.config.scaling.T_scale).speed(1.0));
                ui.label("Scaling L");
                ui.add(DragValue::new(&mut self.config.scaling.L).speed(0.001));
                ui.end_row();

                ui.label("Thermal effects");
                if self.config.thermal_effects.is_empty() {
                    ui.label("No thermal effects configured");
                } else {
                    ui.add(DragValue::new(&mut self.config.thermal_effects[0]).speed(1000.0));
                }
                ui.label("Thermal effects count");
                ui.label(self.config.thermal_effects.len().to_string());
                ui.end_row();
            });

        ui.horizontal_wrapped(|ui| {
            ui.label("Stop condition");
            ui.checkbox(
                &mut self.config.stop_condition_enabled,
                "Enable fuel depletion stop",
            );

            let species_count = self.task.kindata.substances.len();
            if species_count == 0 {
                ui.label("No species available");
            } else {
                if self.config.stop_condition_species_index >= species_count {
                    self.config.stop_condition_species_index = 0;
                }

                ui.add_enabled_ui(self.config.stop_condition_enabled, |ui| {
                    ui.label("Species");
                    ComboBox::from_id_salt("ivp_stop_condition_species")
                        .selected_text(
                            self.task
                                .kindata
                                .substances
                                .get(self.config.stop_condition_species_index)
                                .cloned()
                                .unwrap_or_else(|| "unknown".to_string()),
                        )
                        .show_ui(ui, |ui| {
                            for (index, species) in self.task.kindata.substances.iter().enumerate()
                            {
                                ui.selectable_value(
                                    &mut self.config.stop_condition_species_index,
                                    index,
                                    species,
                                );
                            }
                        });
                    ui.label("Threshold");
                    ui.add(
                        DragValue::new(&mut self.config.stop_condition_threshold)
                            .speed(0.000001)
                            .range(0.0..=f64::INFINITY),
                    );
                });
            }
        });
    }

    fn render_initial_state_section(&mut self, ui: &mut egui::Ui) {
        ui.heading("Initial state");
        Grid::new("ivp_initial_state_grid")
            .num_columns(4)
            .spacing([16.0, 8.0])
            .show(ui, |ui| {
                ui.label("Temperature");
                ui.add(DragValue::new(&mut self.config.initial_state.temperature).speed(1.0));
                ui.label("Heat flux");
                ui.add(DragValue::new(&mut self.config.initial_state.heat_flux).speed(0.01));
                ui.end_row();

                ui.label("Species");
                ui.label(format!("{}", self.task.kindata.substances.len()));
                ui.label("State order");
                ui.label(self.task.kindata.substances.join(", "));
                ui.end_row();
            });

        if self.task.kindata.substances.is_empty() {
            ui.label("No species are loaded yet.");
            return;
        }

        Grid::new("ivp_species_grid")
            .num_columns(3)
            .striped(true)
            .show(ui, |ui| {
                ui.label("Species");
                ui.label("Initial value");
                ui.label("Canonical");
                ui.end_row();

                for species in &self.task.kindata.substances {
                    let value = self
                        .config
                        .initial_state
                        .species
                        .entry(species.clone())
                        .or_insert(0.0);
                    ui.label(species);
                    ui.add(DragValue::new(value).speed(0.01));
                    ui.label("yes");
                    ui.end_row();
                }
            });
    }

    fn render_solver_backend_section(&mut self, ui: &mut egui::Ui) {
        ui.heading("Solver backend");
        Grid::new("ivp_solver_grid")
            .num_columns(4)
            .spacing([16.0, 8.0])
            .show(ui, |ui| {
                ui.label("Method");
                ComboBox::from_id_salt("ivp_method")
                    .selected_text(format!("{:?}", self.config.method))
                    .show_ui(ui, |ui| {
                        ui.selectable_value(
                            &mut self.config.method,
                            ReactorIvpMethod::Auto,
                            "Auto",
                        );
                        ui.selectable_value(&mut self.config.method, ReactorIvpMethod::Bdf, "BDF");
                        ui.selectable_value(
                            &mut self.config.method,
                            ReactorIvpMethod::Adams,
                            "Adams",
                        );
                    });
                ui.label("Residual/Jacobian");
                ComboBox::from_id_salt("ivp_symbolic_backend")
                    .selected_text(format!("{:?}", self.config.symbolic_backend))
                    .show_ui(ui, |ui| {
                        ui.selectable_value(
                            &mut self.config.symbolic_backend,
                            ReactorIvpSymbolicBackend::AtomView,
                            "AtomView",
                        );
                        ui.selectable_value(
                            &mut self.config.symbolic_backend,
                            ReactorIvpSymbolicBackend::ExprLegacy,
                            "ExprLegacy",
                        );
                    });
                ui.end_row();

                ui.label("Matrix backend");
                ComboBox::from_id_salt("ivp_matrix_backend")
                    .selected_text(format!("{:?}", self.config.matrix_backend))
                    .show_ui(ui, |ui| {
                        ui.selectable_value(
                            &mut self.config.matrix_backend,
                            ReactorIvpMatrixBackend::Sparse,
                            "Sparse",
                        );
                        ui.selectable_value(
                            &mut self.config.matrix_backend,
                            ReactorIvpMatrixBackend::Banded,
                            "Banded",
                        );
                    });
                ui.label("Execution backend");
                ComboBox::from_id_salt("ivp_execution_mode")
                    .selected_text(match self.config.execution_mode {
                        ReactorIvpGuiExecutionMode::Lambdify => "Lambdify",
                        ReactorIvpGuiExecutionMode::Aot => "AOT",
                    })
                    .show_ui(ui, |ui| {
                        ui.selectable_value(
                            &mut self.config.execution_mode,
                            ReactorIvpGuiExecutionMode::Lambdify,
                            "Lambdify",
                        );
                        ui.selectable_value(
                            &mut self.config.execution_mode,
                            ReactorIvpGuiExecutionMode::Aot,
                            "AOT",
                        );
                    });
                ui.end_row();

                if matches!(self.config.execution_mode, ReactorIvpGuiExecutionMode::Aot) {
                    ui.end_row();
                    ui.label("AOT settings");
                    ui.label("Compiler and build profile are only relevant for AOT.");
                    ui.end_row();

                    ui.label("toolchain");
                    ComboBox::from_id_salt("ivp_aot_toolchain")
                        .selected_text(format!("{:?}", self.config.aot_toolchain))
                        .show_ui(ui, |ui| {
                            ui.selectable_value(
                                &mut self.config.aot_toolchain,
                                Lsode2AotToolchain::CTcc,
                                "CTcc",
                            );
                            ui.selectable_value(
                                &mut self.config.aot_toolchain,
                                Lsode2AotToolchain::CGcc,
                                "CGcc",
                            );
                            ui.selectable_value(
                                &mut self.config.aot_toolchain,
                                Lsode2AotToolchain::Zig,
                                "Zig",
                            );
                            ui.selectable_value(
                                &mut self.config.aot_toolchain,
                                Lsode2AotToolchain::Rust,
                                "Rust",
                            );
                        });
                    ui.label("profile");
                    ComboBox::from_id_salt("ivp_aot_profile")
                        .selected_text(format!("{:?}", self.config.aot_profile))
                        .show_ui(ui, |ui| {
                            ui.selectable_value(
                                &mut self.config.aot_profile,
                                Lsode2AotProfile::Release,
                                "Release",
                            );
                            ui.selectable_value(
                                &mut self.config.aot_profile,
                                Lsode2AotProfile::Debug,
                                "Debug",
                            );
                        });
                    ui.end_row();
                    ui.label("AOT confirmation");
                    ui.label("AOT runs must be explicitly confirmed before execution.");
                    ui.end_row();
                }
            });
    }

    fn render_integration_section(&mut self, ui: &mut egui::Ui) {
        ui.heading("Integration");
        Grid::new("ivp_integration_grid")
            .num_columns(4)
            .spacing([16.0, 8.0])
            .show(ui, |ui| {
                ui.label("x0");
                ui.add(DragValue::new(&mut self.config.x0).speed(0.01));
                ui.label("x_bound");
                ui.add(DragValue::new(&mut self.config.x_bound).speed(0.01));
                ui.end_row();

                ui.label("first_step");
                let mut first_step_enabled = self.config.first_step.is_some();
                ui.checkbox(&mut first_step_enabled, "");
                if first_step_enabled {
                    let value = self.config.first_step.get_or_insert(1e-6);
                    ui.add(DragValue::new(value).speed(0.000_001));
                } else {
                    self.config.first_step = None;
                    ui.label("None");
                }
                ui.end_row();

                ui.label("max_step");
                ui.add(DragValue::new(&mut self.config.max_step).speed(0.01));
                ui.label("rtol");
                ui.add(DragValue::new(&mut self.config.rtol).speed(0.000001));
                ui.end_row();

                ui.label("atol");
                ui.add(DragValue::new(&mut self.config.atol).speed(0.0000001));
                ui.end_row();
            });
    }

    fn render_postprocessing_section(&mut self, ui: &mut egui::Ui) {
        ui.heading("Postprocessing");
        Grid::new("ivp_postprocessing_grid")
            .num_columns(4)
            .spacing([16.0, 8.0])
            .show(ui, |ui| {
                ui.label("Save TXT");
                ui.checkbox(&mut self.config.postprocessing.save_txt, "");
                ui.label("TXT path");
                ui.add(
                    TextEdit::singleline(&mut self.config.postprocessing.txt_path)
                        .desired_width(360.0)
                        .hint_text("auto"),
                );
                ui.end_row();

                ui.label("Save CSV");
                ui.checkbox(&mut self.config.postprocessing.save_csv, "");
                ui.label("CSV path");
                ui.add(
                    TextEdit::singleline(&mut self.config.postprocessing.csv_path)
                        .desired_width(360.0)
                        .hint_text("auto"),
                );
                ui.end_row();

                ui.label("Write report");
                ui.checkbox(&mut self.config.postprocessing.write_report, "");
                ui.label("Report path");
                ui.add(
                    TextEdit::singleline(&mut self.config.postprocessing.report_path)
                        .desired_width(360.0)
                        .hint_text("auto"),
                );
                ui.end_row();

                ui.label("Plotters PNG");
                ui.checkbox(&mut self.config.postprocessing.plotters_png, "");
                ui.label("Plotters dir");
                ui.add(
                    TextEdit::singleline(&mut self.config.postprocessing.plotters_dir)
                        .desired_width(360.0)
                        .hint_text("auto"),
                );
                ui.end_row();

                ui.label("Gnuplot PNG");
                ui.checkbox(&mut self.config.postprocessing.gnuplot_png, "");
                ui.label("Gnuplot dir");
                ui.add(
                    TextEdit::singleline(&mut self.config.postprocessing.gnuplot_dir)
                        .desired_width(360.0)
                        .hint_text("auto"),
                );
                ui.end_row();

                ui.label("Terminal plot");
                ui.checkbox(&mut self.config.postprocessing.terminal_plot, "");
                ui.label("Selected exports");
                ui.label(if self.config.postprocessing.has_actions() {
                    "configured"
                } else {
                    "none"
                });
                ui.end_row();
            });

        ui.separator();
        Grid::new("ivp_postprocessing_status_grid")
            .num_columns(4)
            .spacing([16.0, 8.0])
            .show(ui, |ui| {
                ui.label("Last status");
                ui.label(self.last_status.as_deref().unwrap_or("Not run yet"));
                ui.label("Last error");
                ui.label(self.last_error.as_deref().unwrap_or("None"));
                ui.end_row();

                ui.label("Preview rows");
                ui.label(
                    self.last_preview
                        .as_ref()
                        .map(|snapshot| snapshot.summary_rows.len().to_string())
                        .unwrap_or_else(|| "0".to_string()),
                );
                ui.label("Equation rows");
                ui.label(
                    self.last_preview
                        .as_ref()
                        .map(|snapshot| snapshot.equation_rows.len().to_string())
                        .unwrap_or_else(|| "0".to_string()),
                );
                ui.end_row();
            });
    }

    fn render_conservation_section(&mut self, ui: &mut egui::Ui) {
        ui.heading("Conservation");
        match self.task.latest_conservation_report() {
            Ok(report) => {
                ui.label("Latest solve conservation summary");
                let table = Table::new(report.rows()).to_string();
                ui.monospace(table);
            }
            Err(err) => {
                ui.label("Latest solve conservation summary");
                ui.colored_label(
                    egui::Color32::from_rgb(120, 120, 120),
                    format!("Not available yet: {err}"),
                );
            }
        }
    }

    fn render_solve_report_section(&mut self, ui: &mut egui::Ui) {
        ui.heading("Solve report");
        match self.task.latest_solve_report_rows() {
            Ok(rows) => {
                ui.label("Latest solve diagnostics");
                let table = Table::new(rows).to_string();
                ui.monospace(table);
            }
            Err(err) => {
                ui.label("Latest solve diagnostics");
                ui.colored_label(
                    egui::Color32::from_rgb(120, 120, 120),
                    format!("Not available yet: {err}"),
                );
            }
        }
    }

    fn render_actions(&mut self, ui: &mut egui::Ui) {
        ui.horizontal(|ui| {
            if ui.button("Preview Task").clicked() {
                match self.preview_task() {
                    Ok(snapshot) => {
                        snapshot.print_to_console();
                        self.last_preview = Some(snapshot);
                        self.last_status = Some("Task preview printed to console.".to_string());
                        self.last_error = None;
                    }
                    Err(err) => {
                        self.last_error = Some(err.to_string());
                    }
                }
            }

            let run_enabled = !self.calculation_state.is_active();
            let run_label = if matches!(self.config.execution_mode, ReactorIvpGuiExecutionMode::Aot)
            {
                "Run AOT calculation"
            } else {
                "Run Calculation"
            };
            if ui
                .add_enabled(run_enabled, egui::Button::new(run_label))
                .clicked()
            {
                self.request_calculation();
            }
        });

        if let Some(snapshot) = &self.last_preview {
            ui.label(format!(
                "Preview: {} summary rows, {} equation rows",
                snapshot.summary_rows.len(),
                snapshot.equation_rows.len()
            ));
        }
    }

    fn request_calculation(&mut self) {
        if self.calculation_state.is_active() {
            self.last_error = Some("A calculation is already running.".to_string());
            return;
        }

        if matches!(self.config.execution_mode, ReactorIvpGuiExecutionMode::Aot) {
            self.request_aot_calculation();
        } else {
            self.start_calculation_worker();
        }
    }

    fn request_aot_calculation(&mut self) {
        let mut preview_task = self.task.clone();
        if let Err(err) = self.config.apply_to_task(&mut preview_task) {
            self.last_error = Some(err.to_string());
            self.last_status = Some("AOT calculation rejected before confirmation.".to_string());
            return;
        }

        let request = ReactorIvpAotRequest::from_config(&self.config);
        self.last_status = Some(format!("AOT run pending confirmation: {}", request.summary));
        self.last_error = None;
        self.pending_aot_confirmation = Some(request);
    }

    fn confirm_pending_aot_calculation(&mut self) {
        let Some(request) = self.pending_aot_confirmation.take() else {
            return;
        };

        let current_fingerprint = Self::fingerprint_config(&self.config);
        if current_fingerprint != request.config_fingerprint {
            self.last_status =
                Some("AOT confirmation was cancelled because the editor changed.".to_string());
            self.last_error = Some("AOT request became stale before confirmation.".to_string());
            return;
        }

        self.last_status = Some(format!("AOT confirmed: {}", request.summary));
        self.start_calculation_worker();
    }

    fn cancel_pending_aot_calculation(&mut self) {
        if self.pending_aot_confirmation.take().is_some() {
            self.last_status = Some("AOT calculation cancelled before execution.".to_string());
            self.last_error = None;
        }
    }

    fn show_aot_confirmation(&mut self, ctx: &egui::Context) {
        let Some(request) = self.pending_aot_confirmation.clone() else {
            return;
        };

        let mut confirm = false;
        let mut cancel = false;
        egui::Window::new("Confirm AOT run")
            .collapsible(false)
            .resizable(false)
            .show(ctx, |ui| {
                ui.label("AOT execution will lower the symbolic system into generated code.");
                ui.label("Confirm only if you want the external compiler route to run now.");
                ui.monospace(request.summary.clone());
                ui.horizontal(|ui| {
                    if ui.button("Cancel AOT run").clicked() {
                        cancel = true;
                    }
                    if ui.button("Confirm AOT run").clicked() {
                        confirm = true;
                    }
                });
            });

        if cancel {
            self.cancel_pending_aot_calculation();
        } else if confirm {
            self.confirm_pending_aot_calculation();
        }
    }

    fn preview_task(&self) -> Result<ReactorIvpPreviewSnapshot, IvpError> {
        build_reactor_ivp_task_preview_snapshot(&self.task, &self.config)
    }

    fn start_calculation_worker(&mut self) {
        if self.calculation_state.is_active() {
            self.last_error = Some("A calculation is already running.".to_string());
            return;
        }

        let mut task = self.task.clone();
        let config = self.config.clone();
        let config_fingerprint = Self::fingerprint_config(&config);
        let run_id = self.next_run_id;
        self.next_run_id = self.next_run_id.wrapping_add(1);
        let (sender, receiver) = mpsc::channel();
        self.last_status = Some("Running calculation...".to_string());
        self.last_error = None;
        self.calculation_state = ReactorIvpCalculationState::Running { run_id };
        self.calculation_receiver = Some(receiver);

        thread::Builder::new()
            .name(format!("kithe-reactor-ivp-{run_id}"))
            .spawn(move || {
                let outcome = (|| -> Result<_, IvpError> {
                    config.apply_to_task(&mut task)?;
                    let result = task.solve()?;
                    let render_data = task.solution_render_data()?;
                    let preview_rows = render_data.solution_preview_rows()?;
                    let postprocessing_summary = run_reactor_ivp_postprocessing(
                        &task,
                        &config,
                        &render_data.unknowns,
                        &render_data.x_mesh,
                        &render_data.solution,
                    )?;
                    let status = summarize_reactor_ivp_status(&result, postprocessing_summary);
                    let preview = ReactorIvpPreviewSnapshot {
                        summary_rows: vec![preview_row("solve", "status", status.clone())],
                        equation_rows: preview_rows
                            .into_iter()
                            .enumerate()
                            .map(|(index, row)| ReactorIvpEquationRow {
                                index,
                                variable: row.variable,
                                equation: format!(
                                    "first={} middle={} last={}",
                                    row.first, row.middle, row.last
                                ),
                            })
                            .collect(),
                    };
                    Ok((
                        task,
                        render_data.unknowns,
                        render_data.x_mesh,
                        render_data.solution,
                        preview,
                        status,
                    ))
                })();
                let _ = sender.send(ReactorIvpWorkerMessage {
                    run_id,
                    config_fingerprint,
                    outcome,
                });
            })
            .expect("reactor IVP worker thread should start");
    }

    /// Start a deterministic worker used only by tests.
    #[cfg(test)]
    pub(crate) fn start_test_calculation_worker(&mut self, delay: Duration, succeed: bool) {
        if self.calculation_state.is_active() {
            self.last_error = Some("A calculation is already running.".to_string());
            return;
        }

        let mut task = self.task.clone();
        let config = self.config.clone();
        let config_fingerprint = Self::fingerprint_config(&config);
        let run_id = self.next_run_id;
        self.next_run_id = self.next_run_id.wrapping_add(1);
        let (sender, receiver) = mpsc::channel();
        self.last_status = Some("Running calculation...".to_string());
        self.last_error = None;
        self.calculation_state = ReactorIvpCalculationState::Running { run_id };
        self.calculation_receiver = Some(receiver);

        thread::Builder::new()
            .name(format!("kithe-reactor-ivp-test-{run_id}"))
            .spawn(move || {
                thread::sleep(delay);
                let outcome = if succeed {
                    (|| -> Result<_, IvpError> {
                        config.apply_to_task(&mut task)?;
                        let result = task.solve()?;
                        let render_data = task.solution_render_data()?;
                        let preview_rows = render_data.solution_preview_rows()?;
                        let postprocessing_summary = run_reactor_ivp_postprocessing(
                            &task,
                            &config,
                            &render_data.unknowns,
                            &render_data.x_mesh,
                            &render_data.solution,
                        )?;
                        let status = summarize_reactor_ivp_status(&result, postprocessing_summary);
                        let preview = ReactorIvpPreviewSnapshot {
                            summary_rows: vec![preview_row("solve", "status", status.clone())],
                            equation_rows: preview_rows
                                .into_iter()
                                .enumerate()
                                .map(|(index, row)| ReactorIvpEquationRow {
                                    index,
                                    variable: row.variable,
                                    equation: format!(
                                        "first={} middle={} last={}",
                                        row.first, row.middle, row.last
                                    ),
                                })
                                .collect(),
                        };
                        Ok((
                            task,
                            render_data.unknowns,
                            render_data.x_mesh,
                            render_data.solution,
                            preview,
                            status,
                        ))
                    })()
                } else {
                    Err(IvpError::CalculationError(
                        "synthetic reactor IVP worker failure".to_string(),
                    ))
                };
                let _ = sender.send(ReactorIvpWorkerMessage {
                    run_id,
                    config_fingerprint,
                    outcome,
                });
            })
            .expect("reactor IVP test worker thread should start");
    }

    /// Compute a stable fingerprint for the currently edited solver inputs.
    ///
    /// This is used to reject stale worker results if the user modifies the
    /// editor while a background solve is still in flight.
    fn fingerprint_config(config: &IvpGuiConfig) -> u64 {
        let mut hasher = DefaultHasher::new();

        config.problem_name.hash(&mut hasher);
        config.problem_description.hash(&mut hasher);

        hasher.write_u64(config.physical.ro.to_bits());
        hasher.write_u64(config.physical.Cp.to_bits());
        hasher.write_u64(config.physical.Lambda.to_bits());
        hasher.write_u64(config.physical.m.to_bits());
        hasher.write_u64(config.physical.L.to_bits());

        hasher.write_u64(config.scaling.dT.to_bits());
        hasher.write_u64(config.scaling.L.to_bits());
        hasher.write_u64(config.scaling.T_scale.to_bits());

        hasher.write_u64(config.initial_state.temperature.to_bits());
        hasher.write_u64(config.initial_state.heat_flux.to_bits());

        let mut species_entries = config
            .initial_state
            .species
            .iter()
            .map(|(name, value)| (name.clone(), value.to_bits()))
            .collect::<Vec<_>>();
        species_entries.sort_by(|left, right| left.0.cmp(&right.0));
        for (name, value_bits) in species_entries {
            name.hash(&mut hasher);
            hasher.write_u64(value_bits);
        }

        for effect in &config.thermal_effects {
            hasher.write_u64(effect.to_bits());
        }

        config.method.hash(&mut hasher);
        config.symbolic_backend.hash(&mut hasher);
        config.matrix_backend.hash(&mut hasher);
        config.execution_mode.hash(&mut hasher);
        format!("{:?}", config.aot_toolchain).hash(&mut hasher);
        format!("{:?}", config.aot_profile).hash(&mut hasher);

        hasher.write_u64(config.x0.to_bits());
        hasher.write_u64(config.x_bound.to_bits());
        match config.first_step {
            Some(value) => {
                1u8.hash(&mut hasher);
                hasher.write_u64(value.to_bits());
            }
            None => 0u8.hash(&mut hasher),
        }
        hasher.write_u64(config.max_step.to_bits());
        hasher.write_u64(config.rtol.to_bits());
        hasher.write_u64(config.atol.to_bits());
        config.postprocessing.save_txt.hash(&mut hasher);
        config.postprocessing.txt_path.hash(&mut hasher);
        config.postprocessing.save_csv.hash(&mut hasher);
        config.postprocessing.csv_path.hash(&mut hasher);
        config.postprocessing.write_report.hash(&mut hasher);
        config.postprocessing.report_path.hash(&mut hasher);
        config.postprocessing.plotters_png.hash(&mut hasher);
        config.postprocessing.plotters_dir.hash(&mut hasher);
        config.postprocessing.gnuplot_png.hash(&mut hasher);
        config.postprocessing.gnuplot_dir.hash(&mut hasher);
        config.postprocessing.terminal_plot.hash(&mut hasher);
        format!("{:?}", config.native_execution).hash(&mut hasher);

        hasher.finish()
    }

    fn apply_worker_outcome(
        &mut self,
        current_fingerprint: u64,
        expected_fingerprint: u64,
        outcome: Result<
            (
                SimpleReactorTask,
                Vec<String>,
                nalgebra::DVector<f64>,
                nalgebra::DMatrix<f64>,
                ReactorIvpPreviewSnapshot,
                String,
            ),
            IvpError,
        >,
    ) {
        if current_fingerprint != expected_fingerprint {
            self.last_status = Some(
                "Calculation result ignored because the editor changed during the run.".to_string(),
            );
            self.last_error =
                Some("Stale calculation result discarded after editor changes.".to_string());
            self.calculation_state = ReactorIvpCalculationState::Failed;
            return;
        }

        match outcome {
            Ok((task, variable_names, axis, values, preview, status)) => {
                self.task = task;
                self.plot_window = Some(PlotWindow::new(
                    self.task.solver.arg_name.clone(),
                    variable_names,
                    axis,
                    values,
                ));
                self.last_preview = Some(preview);
                self.last_status = Some(format!("Calculation completed with status: {status}"));
                self.last_error = None;
                self.calculation_state = ReactorIvpCalculationState::Completed;
            }
            Err(err) => {
                self.last_status = Some(format!("Calculation failed: {err}"));
                self.last_error = Some(err.to_string());
                self.calculation_state = ReactorIvpCalculationState::Failed;
            }
        }
    }

    fn poll_calculation_worker(&mut self, ctx: &egui::Context) {
        let Some(receiver) = &self.calculation_receiver else {
            return;
        };

        match receiver.try_recv() {
            Ok(message) => {
                let active_run_id = self.calculation_state.active_run_id();
                self.calculation_receiver = None;
                if active_run_id == Some(message.run_id) {
                    let current_fingerprint = Self::fingerprint_config(&self.config);
                    self.apply_worker_outcome(
                        current_fingerprint,
                        message.config_fingerprint,
                        message.outcome,
                    );
                }
                ctx.request_repaint();
            }
            Err(mpsc::TryRecvError::Empty) => {
                ctx.request_repaint_after(std::time::Duration::from_millis(100));
            }
            Err(mpsc::TryRecvError::Disconnected) => {
                self.calculation_receiver = None;
                self.calculation_state = ReactorIvpCalculationState::Failed;
                self.last_status =
                    Some("Calculation failed: worker disconnected without a result.".to_string());
                self.last_error =
                    Some("Calculation worker disconnected without a result.".to_string());
                ctx.request_repaint();
            }
        }
    }
}

fn build_demo_task() -> Result<SimpleReactorTask, IvpError> {
    let mut task = SimpleReactorTask::new();
    task.kindata
        .set_reaction_data_directly(
            vec![ReactionData::new_elementary(
                "A=>B".to_string(),
                vec![1.0e3, 0.0, 1.0e4],
                None,
            )],
            None,
        )
        .map_err(|err| IvpError::InvalidConfiguration(err.to_string()))?;
    task.kindata.stecheodata.vec_of_molmasses = Some(vec![10.0, 20.0]);
    task.kindata.substances = vec!["A".to_string(), "B".to_string()];
    task.problem_name = Some("CondensedBurn".to_string());
    task.problem_description = Some("Canonical condensed-phase IVP demo".to_string());
    task.set_density(1200.0)?;
    task.set_transport_properties(0.25, 1000.0);
    task.m = 0.015;
    task.L = 0.05;
    task.set_scaling(ScalingConfig::new(100.0, 0.05, 100.0))
        .map_err(IvpError::from)?;
    task.set_initial_conditions(HashMap::from([
        ("T".to_string(), 450.0),
        ("q".to_string(), 0.15),
        ("A".to_string(), 0.7),
        ("B".to_string(), 0.3),
    ]));
    task.set_thermal_effects(vec![-2.0e5]);
    task.set_solver_backend_config(
        ReactorIvpSolverConfig::default()
            .with_integration_domain(0.0, 0.01)
            .with_max_step(1.0e-3)
            .with_native_execution(Lsode2NativeExecutionConfig::bridge_solve()),
    );
    Ok(task)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn demo_preview_snapshot_prints_equations_and_solver_rows() {
        let task = build_demo_task().expect("demo task should seed");
        let config = IvpGuiConfig::from_task(&task);
        let snapshot =
            build_reactor_ivp_task_preview_snapshot(&task, &config).expect("preview should build");

        assert!(
            snapshot
                .summary_rows
                .iter()
                .any(|row| row.section == "solver" && row.field == "matrix_backend")
        );
        assert!(
            snapshot
                .equation_rows
                .iter()
                .any(|row| row.variable == "Teta")
        );
        assert!(
            snapshot
                .summary_rows
                .iter()
                .any(|row| row.section == "initial_state" && row.field == "temperature")
        );
        assert!(
            snapshot
                .summary_rows
                .iter()
                .any(|row| row.section == "integration" && row.field == "x_bound")
        );
    }

    #[test]
    fn gui_config_roundtrips_solver_backend_choices() {
        let mut config = IvpGuiConfig::default();
        config.method = ReactorIvpMethod::Adams;
        config.matrix_backend = ReactorIvpMatrixBackend::Banded;
        config.symbolic_backend = ReactorIvpSymbolicBackend::ExprLegacy;
        config.execution_mode = ReactorIvpGuiExecutionMode::Aot;
        config.aot_toolchain = Lsode2AotToolchain::CTcc;
        config.aot_profile = Lsode2AotProfile::Debug;

        let solver = config.to_solver_config();
        assert_eq!(solver.method, ReactorIvpMethod::Adams);
        assert_eq!(solver.matrix_backend, ReactorIvpMatrixBackend::Banded);
        assert_eq!(
            solver.symbolic_backend,
            ReactorIvpSymbolicBackend::ExprLegacy
        );
        assert!(matches!(
            solver.execution_backend,
            ReactorIvpExecutionBackend::Aot {
                toolchain: Lsode2AotToolchain::CTcc,
                profile: Lsode2AotProfile::Debug
            }
        ));
    }
}
