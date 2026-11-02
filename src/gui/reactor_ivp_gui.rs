//! Condensed reactor IVP GUI.
//!
//! This screen is the reactor-IVP counterpart to the BVP combustion editor.
//! It keeps the user-facing solver surface typed and explicit:
//! - physics and initial state are edited separately from the solver backend;
//! - Lambdify and AOT are mutually exclusive execution routes;
//! - Sparse and Banded are the two matrix-storage choices shown to the user;
//! - postprocessing controls choose how solved results are exported;
//! - preview prints a tabulated snapshot plus the pretty-printed symbolic RHS.

use crate::Kinetics::mechfinder_api::{ReactionData, ReactionKinetics};

use crate::ReactorsBVP::reactor_BVP_utils::ScalingConfig;
use crate::ReactorsIVP::SimpleReactorIVP::{
    IvpError, IvpSolveSnapshot, ReactorIvpInitialStateConfig, ReactorIvpPhysicalConfig,
    SimpleReactorTask,
};
use crate::ReactorsIVP::solver_backend::{
    ReactorIvpExecutionBackend, ReactorIvpFamilyOptions, ReactorIvpMatrixBackend, ReactorIvpMethod,
    ReactorIvpSolverConfig, ReactorIvpSolverFamily, ReactorIvpStopCondition,
    ReactorIvpSymbolicBackend,
};
use crate::gui::document_lifecycle::DocumentLifecycleState;
use crate::gui::gui_plot::PlotWindow;
use RustedSciThe::Utils::postprocessing::{PostprocessDataset, PostprocessPlan, PostprocessStatus};
use egui::{ComboBox, DragValue, Grid, Popup, TextEdit};
use log::{info, warn};

use RustedSciThe::numerical::LSODE2::{
    Lsode2AotProfile, Lsode2AotToolchain, Lsode2NativeExecutionConfig,
};
use RustedSciThe::numerical::SolverTelemetryMode;
use RustedSciThe::numerical::logging::SolverLogLevel;
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

fn handle_combo_keyboard<T: Copy + PartialEq>(
    ctx: &egui::Context,
    response: &egui::Response,
    selected: &mut T,
    options: &[T],
) {
    if !response.has_focus() || options.is_empty() {
        return;
    }
    response.clone().highlight();
    let (activate, previous, next) = ctx.input(|input| {
        (
            input.key_pressed(egui::Key::Enter) || input.key_pressed(egui::Key::Space),
            input.key_pressed(egui::Key::ArrowUp),
            input.key_pressed(egui::Key::ArrowDown),
        )
    });
    let popup_id = response.id.with("popup");
    if activate {
        if response.clicked() {
            // egui already toggled the popup from its native activation.
        } else if Popup::is_id_open(ctx, popup_id) {
            Popup::close_id(ctx, popup_id);
        } else {
            Popup::open_id(ctx, popup_id);
        }
    } else if previous || next {
        let index = options
            .iter()
            .position(|option| option == selected)
            .unwrap_or(0);
        let next_index = if next {
            (index + 1) % options.len()
        } else {
            (index + options.len() - 1) % options.len()
        };
        *selected = options[next_index];
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub enum ReactorIvpGuiCallbackPolicy {
    #[default]
    Sequential,
    Parallel,
    Auto,
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
    pub substances: Vec<String>,
    pub reactions: Vec<ReactionData>,
    pub thermal_effects: Vec<f64>,
    pub stop_condition_enabled: bool,
    pub stop_condition_species_index: usize,
    pub stop_condition_threshold: f64,
    pub solver_family: ReactorIvpSolverFamily,
    pub method: ReactorIvpMethod,
    pub symbolic_backend: ReactorIvpSymbolicBackend,
    pub matrix_backend: ReactorIvpMatrixBackend,
    pub execution_mode: ReactorIvpGuiExecutionMode,
    pub callback_policy: ReactorIvpGuiCallbackPolicy,
    pub callback_min_work: usize,
    pub aot_toolchain: Lsode2AotToolchain,
    pub aot_profile: Lsode2AotProfile,
    pub bdf_max_steps: usize,
    pub radau_max_steps: usize,
    pub radau_max_newton_iterations: usize,
    pub x0: f64,
    pub x_bound: f64,
    pub first_step: Option<f64>,
    pub max_step: f64,
    pub rtol: f64,
    pub atol: f64,
    pub telemetry: SolverTelemetryMode,
    pub logging: SolverLogLevel,
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
            substances: Vec::new(),
            reactions: Vec::new(),
            thermal_effects: vec![-2.0e5],
            stop_condition_enabled: false,
            stop_condition_species_index: 0,
            stop_condition_threshold: 1e-4,
            solver_family: ReactorIvpSolverFamily::Lsode2,
            method: ReactorIvpMethod::default(),
            symbolic_backend: ReactorIvpSymbolicBackend::default(),
            matrix_backend: ReactorIvpMatrixBackend::default(),
            execution_mode: ReactorIvpGuiExecutionMode::default(),
            callback_policy: ReactorIvpGuiCallbackPolicy::default(),
            callback_min_work: 128,
            aot_toolchain: Lsode2AotToolchain::default(),
            aot_profile: Lsode2AotProfile::default(),
            bdf_max_steps: 200_000,
            radau_max_steps: 200_000,
            radau_max_newton_iterations: 6,
            x0: 0.0,
            x_bound: 0.01,
            first_step: None,
            max_step: 1.0e-3,
            rtol: 1e-6,
            atol: 1e-8,
            telemetry: SolverTelemetryMode::Off,
            logging: SolverLogLevel::Off,
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

        let family_options = match self.solver_family {
            ReactorIvpSolverFamily::Lsode2 => ReactorIvpFamilyOptions::Lsode2 {
                controller: self.method,
                native_execution: self.native_execution,
            },
            ReactorIvpSolverFamily::Bdf => ReactorIvpFamilyOptions::Bdf {
                max_steps: self.bdf_max_steps,
            },
            ReactorIvpSolverFamily::Radau => ReactorIvpFamilyOptions::Radau {
                max_steps: self.radau_max_steps,
                max_newton_iterations: self.radau_max_newton_iterations,
            },
        };

        let callback_policy = match self.callback_policy {
            ReactorIvpGuiCallbackPolicy::Sequential => {
                crate::ReactorsIVP::solver_backend::ReactorIvpCallbackPolicy::Sequential
            }
            ReactorIvpGuiCallbackPolicy::Parallel => {
                crate::ReactorsIVP::solver_backend::ReactorIvpCallbackPolicy::Parallel {
                    min_work: self.callback_min_work,
                }
            }
            ReactorIvpGuiCallbackPolicy::Auto => {
                crate::ReactorsIVP::solver_backend::ReactorIvpCallbackPolicy::Auto {
                    min_work: self.callback_min_work,
                }
            }
        };

        ReactorIvpSolverConfig {
            family_options,
            symbolic_backend: self.symbolic_backend,
            matrix_backend: self.matrix_backend,
            execution_backend,
            callback_policy: (self.execution_mode == ReactorIvpGuiExecutionMode::Lambdify)
                .then_some(callback_policy),
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
            telemetry: self.telemetry,
            logging: self.logging,
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
        config.substances = task.kindata.substances().to_vec();
        config.reactions = task.kindata.reaction_data().unwrap_or_default().to_vec();
        config.thermal_effects = task.thermal_effects.clone();
        config.solver_family = task.solver_backend_config.family();
        match task.solver_backend_config.family_options {
            ReactorIvpFamilyOptions::Bdf { max_steps } => config.bdf_max_steps = max_steps,
            ReactorIvpFamilyOptions::Radau {
                max_steps,
                max_newton_iterations,
            } => {
                config.radau_max_steps = max_steps;
                config.radau_max_newton_iterations = max_newton_iterations;
            }
            ReactorIvpFamilyOptions::Lsode2 { .. } => {}
        }
        if let Some(policy) = task.solver_backend_config.callback_policy {
            match policy {
                crate::ReactorsIVP::solver_backend::ReactorIvpCallbackPolicy::Sequential => {
                    config.callback_policy = ReactorIvpGuiCallbackPolicy::Sequential;
                }
                crate::ReactorsIVP::solver_backend::ReactorIvpCallbackPolicy::Parallel {
                    min_work,
                } => {
                    config.callback_policy = ReactorIvpGuiCallbackPolicy::Parallel;
                    config.callback_min_work = min_work;
                }
                crate::ReactorsIVP::solver_backend::ReactorIvpCallbackPolicy::Auto { min_work } => {
                    config.callback_policy = ReactorIvpGuiCallbackPolicy::Auto;
                    config.callback_min_work = min_work;
                }
            }
        }
        if let Some(stop_condition) = task.solver_backend_config.stop_condition {
            config.stop_condition_enabled = true;
            config.stop_condition_species_index = stop_condition.species_index;
            config.stop_condition_threshold = stop_condition.threshold;
        }
        if let Some(controller) = task.solver_backend_config.lsode2_controller() {
            config.method = controller;
        }
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
        config.telemetry = task.solver_backend_config.telemetry;
        config.logging = task.solver_backend_config.logging;
        if let ReactorIvpFamilyOptions::Lsode2 {
            native_execution, ..
        } = task.solver_backend_config.family_options
        {
            config.native_execution = native_execution;
        }
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
        if self.reactions.len() != self.thermal_effects.len() {
            return Err(IvpError::InvalidConfiguration(format!(
                "reaction count ({}) must match thermal-effect count ({})",
                self.reactions.len(),
                self.thermal_effects.len()
            )));
        }
        let groups = task.kindata.groups().cloned();
        task.kindata
            .set_reaction_data_directly(self.reactions.clone(), groups)
            .map_err(|error| {
                IvpError::InvalidConfiguration(format!("invalid reaction set: {error}"))
            })?;
        task.kindata.substances = self.substances.clone();
        task.kindata
            .stecheodata
            .set_substances(self.substances.clone());
        if !self.reactions.is_empty() {
            task.kindata.analyze_reactions().map_err(|error| {
                IvpError::InvalidConfiguration(format!("reaction analysis failed: {error}"))
            })?;
        }

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

        let mut initial_conditions = HashMap::with_capacity(2 + self.substances.len());
        initial_conditions.insert("T".to_string(), self.initial_state.temperature);
        initial_conditions.insert("q".to_string(), self.initial_state.heat_flux);
        for substance in &self.substances {
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
        task.set_hybrid_solver_settings(None);
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
    /// Format the preview as deterministic console-ready tables.
    pub fn format_console_report(&self) -> String {
        format!(
            "Reactor IVP preview:\n{}\n\nReactor IVP equations:\n{}",
            Table::new(self.summary_rows.clone()),
            Table::new(self.equation_rows.clone())
        )
    }

    /// Print the preview directly to stdout regardless of logger configuration.
    pub fn print_to_console(&self) {
        let report = self.format_console_report();
        println!("{report}");
        info!("{report}");
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
            "family",
            format!("{:?}", preview_task.solver_backend_config.family()),
        ),
        preview_row(
            "solver",
            "telemetry",
            format!("{:?}", preview_task.solver_backend_config.telemetry),
        ),
        preview_row(
            "solver",
            "logging",
            format!("{:?}", preview_task.solver_backend_config.logging),
        ),
        preview_row(
            "solver",
            "method",
            preview_task
                .solver_backend_config
                .lsode2_controller()
                .map(|method| format!("{method:?}"))
                .unwrap_or_else(|| "family-native".to_string()),
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
            "substances",
            preview_task.kindata.substances.join(", "),
        ),
        preview_row(
            "chemistry",
            "reactions",
            preview_task.kindata.vec_of_equations.join("; "),
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

fn equation_references_species(equation: &str, species: &str) -> bool {
    if species.is_empty() {
        return false;
    }
    equation.match_indices(species).any(|(start, matched)| {
        let end = start + matched.len();
        let before = equation[..start].chars().next_back();
        let after = equation[end..].chars().next();
        let valid_before = before.is_none_or(|ch| {
            ch.is_whitespace() || ch.is_ascii_digit() || matches!(ch, '.' | '+' | '=' | '<' | '>')
        });
        let valid_after =
            after.is_none_or(|ch| ch.is_whitespace() || matches!(ch, '+' | '=' | '<' | '>'));
        valid_before && valid_after
    })
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
        .as_ref()
        .and_then(|summary| summary.native_termination_kind.as_deref())
        .is_some_and(|kind| kind == "reached_stop_condition")
    {
        "stopped_by_condition".to_string()
    } else {
        result.status.clone()
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
    new_substance_name: String,
    new_reaction_equation: String,
    new_reaction_arrhenius: [f64; 3],
    new_reaction_heat: f64,
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
                    new_substance_name: String::new(),
                    new_reaction_equation: String::new(),
                    new_reaction_arrhenius: [1.0e3, 0.0, 1.0e4],
                    new_reaction_heat: 0.0,
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
                    new_substance_name: String::new(),
                    new_reaction_equation: String::new(),
                    new_reaction_arrhenius: [1.0e3, 0.0, 1.0e4],
                    new_reaction_heat: 0.0,
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

    #[allow(dead_code)]
    pub(crate) fn has_plot_window(&self) -> bool {
        self.plot_window.is_some()
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
            let name_label = ui.label("Problem name");
            let name_input =
                ui.add(TextEdit::singleline(&mut self.config.problem_name).desired_width(220.0));
            let _ = name_input.labelled_by(name_label.id);
            let description_label = ui.label("Description");
            let description_input = ui.add(
                TextEdit::singleline(&mut self.config.problem_description).desired_width(520.0),
            );
            let _ = description_input.labelled_by(description_label.id);
        });
    }

    fn render_physics_section(&mut self, ui: &mut egui::Ui) {
        ui.heading("Physics");
        Grid::new("ivp_physics_grid")
            .num_columns(4)
            .spacing([16.0, 8.0])
            .show(ui, |ui| {
                let label = ui.label("Density ro");
                let input = ui.add(DragValue::new(&mut self.config.physical.ro).speed(1.0));
                let _ = input.labelled_by(label.id);
                let label = ui.label("Heat capacity Cp");
                let input = ui.add(DragValue::new(&mut self.config.physical.Cp).speed(1.0));
                let _ = input.labelled_by(label.id);
                ui.end_row();

                let label = ui.label("Conductivity Lambda");
                let input = ui.add(DragValue::new(&mut self.config.physical.Lambda).speed(0.01));
                let _ = input.labelled_by(label.id);
                let label = ui.label("Mass flux m");
                let input = ui.add(DragValue::new(&mut self.config.physical.m).speed(0.001));
                let _ = input.labelled_by(label.id);
                ui.end_row();

                let label = ui.label("Length L");
                let input = ui.add(DragValue::new(&mut self.config.physical.L).speed(0.001));
                let _ = input.labelled_by(label.id);
                let label = ui.label("dT");
                let input = ui.add(DragValue::new(&mut self.config.scaling.dT).speed(1.0));
                let _ = input.labelled_by(label.id);
                ui.end_row();

                let label = ui.label("T_scale");
                let input = ui.add(DragValue::new(&mut self.config.scaling.T_scale).speed(1.0));
                let _ = input.labelled_by(label.id);
                let label = ui.label("Scaling L");
                let input = ui.add(DragValue::new(&mut self.config.scaling.L).speed(0.001));
                let _ = input.labelled_by(label.id);
                ui.end_row();

                let thermal_label = ui.label("Thermal effects");
                if self.config.thermal_effects.is_empty() {
                    ui.label("No thermal effects configured");
                } else {
                    let input =
                        ui.add(DragValue::new(&mut self.config.thermal_effects[0]).speed(1000.0));
                    let _ = input.labelled_by(thermal_label.id);
                }
                ui.label("Thermal effects count");
                ui.label(self.config.thermal_effects.len().to_string());
                ui.end_row();
            });

        self.render_mechanism_section(ui);

        ui.horizontal_wrapped(|ui| {
            ui.label("Stop condition")
                .on_hover_text("Stop when the selected species reaches its threshold. This control is supported by LSODE2 only.");
            let stop_supported = self.config.solver_family == ReactorIvpSolverFamily::Lsode2;
            ui.add_enabled_ui(stop_supported, |ui| {
                ui.checkbox(
                    &mut self.config.stop_condition_enabled,
                    "Enable fuel depletion stop",
                );
            });
            if !stop_supported {
                ui.label("Available for LSODE2 only");
            }

            let species_count = self.config.substances.len();
            if species_count == 0 {
                ui.label("No species available");
            } else {
                if self.config.stop_condition_species_index >= species_count {
                    self.config.stop_condition_species_index = 0;
                }

                ui.add_enabled_ui(self.config.stop_condition_enabled, |ui| {
                    let label = ui.label("Species");
                    let response = ComboBox::from_id_salt("ivp_stop_condition_species")
                        .selected_text(
                            self.config
                                .substances
                                .get(self.config.stop_condition_species_index)
                                .cloned()
                                .unwrap_or_else(|| "unknown".to_string()),
                        )
                        .show_ui(ui, |ui| {
                            for (index, species) in self.config.substances.iter().enumerate()
                            {
                                ui.selectable_value(
                                    &mut self.config.stop_condition_species_index,
                                    index,
                                    species,
                                );
                            }
                        });
                    let _ = response.response.clone().labelled_by(label.id);
                    handle_combo_keyboard(
                        ui.ctx(),
                        &response.response,
                        &mut self.config.stop_condition_species_index,
                        &(0..species_count).collect::<Vec<_>>(),
                    );
                    let threshold_label = ui.label("Threshold");
                    let threshold = ui.add(
                        DragValue::new(&mut self.config.stop_condition_threshold)
                            .speed(0.000001)
                            .range(0.0..=f64::INFINITY),
                    );
                    let _ = threshold.labelled_by(threshold_label.id);
                });
            }
        });
    }

    fn render_mechanism_section(&mut self, ui: &mut egui::Ui) {
        egui::Frame::group(ui.style()).show(ui, |ui| {
            ui.heading("Substances");
            if self.config.substances.is_empty() {
                ui.label("No substances configured.");
            }
            let mut remove_substance = None;
            for (index, species) in self.config.substances.iter().enumerate() {
                ui.horizontal(|ui| {
                    ui.label(species);
                    let used = self
                        .config
                        .reactions
                        .iter()
                        .any(|reaction| equation_references_species(&reaction.eq, species));
                    let remove = ui
                        .add_enabled(!used, egui::Button::new("Remove substance"))
                        .on_hover_text(if used {
                            "Remove reactions that use this substance before deleting it."
                        } else {
                            "Delete this substance and its initial-state value."
                        });
                    if remove.clicked() {
                        remove_substance = Some(index);
                    }
                });
            }
            if let Some(index) = remove_substance {
                let species = self.config.substances.remove(index);
                self.config.initial_state.species.remove(&species);
                if self.config.stop_condition_species_index >= self.config.substances.len() {
                    self.config.stop_condition_species_index =
                        self.config.substances.len().saturating_sub(1);
                }
                self.last_status = Some(format!("Removed substance `{species}`."));
                self.last_error = None;
            }

            ui.horizontal(|ui| {
                let label = ui.label("New substance");
                let input = ui.text_edit_singleline(&mut self.new_substance_name);
                let _ = input.clone().labelled_by(label.id);
                if ui.button("Add substance").clicked() {
                    let name = self.new_substance_name.trim().to_string();
                    if name.is_empty() {
                        self.last_error = Some("Substance name cannot be empty.".to_string());
                    } else if self
                        .config
                        .substances
                        .iter()
                        .any(|existing| existing == &name)
                    {
                        self.last_error = Some(format!("Substance `{name}` already exists."));
                    } else {
                        self.config.substances.push(name.clone());
                        self.config.initial_state.species.insert(name.clone(), 0.0);
                        self.new_substance_name.clear();
                        self.last_status = Some(format!("Added substance `{name}`."));
                        self.last_error = None;
                    }
                }
            });
        });

        ui.add_space(6.0);
        egui::Frame::group(ui.style()).show(ui, |ui| {
            ui.heading("Reactions");
            if self.config.reactions.is_empty() {
                ui.label("No reactions configured.");
            }
            let mut remove_reaction = None;
            for index in 0..self.config.reactions.len() {
                ui.push_id(("ivp_reaction", index), |ui| {
                    egui::Frame::group(ui.style()).show(ui, |ui| {
                        let reaction = &mut self.config.reactions[index];
                        ui.horizontal(|ui| {
                            let label = ui.label(format!("Reaction {}", index + 1));
                            let input = ui.text_edit_singleline(&mut reaction.eq);
                            let _ = input.clone().labelled_by(label.id);
                            if ui.button("Remove reaction").clicked() {
                                remove_reaction = Some(index);
                            }
                        });
                        if let ReactionKinetics::Elementary(data) = &mut reaction.data {
                            ui.horizontal(|ui| {
                                for (parameter_index, name) in
                                    ["A", "n", "E (J/mol)"].into_iter().enumerate()
                                {
                                    if let Some(value) = data.Arrenius.get_mut(parameter_index) {
                                        let label = ui.label(name);
                                        let input = ui.add(DragValue::new(value).speed(0.1));
                                        let _ = input.labelled_by(label.id);
                                    }
                                }
                            });
                        } else {
                            ui.label(format!(
                                "Rate model: {:?}; model-specific parameters are preserved.",
                                reaction.reaction_type
                            ));
                        }
                        if let Some(heat) = self.config.thermal_effects.get_mut(index) {
                            let label = ui.label("Heat effect (J/mol)");
                            let input = ui.add(DragValue::new(heat).speed(100.0));
                            let _ = input.labelled_by(label.id);
                        }
                    });
                });
            }
            if let Some(index) = remove_reaction {
                let equation = self.config.reactions.remove(index).eq;
                if index < self.config.thermal_effects.len() {
                    self.config.thermal_effects.remove(index);
                }
                self.last_status = Some(format!("Removed reaction `{equation}`."));
                self.last_error = None;
            }

            ui.separator();
            ui.horizontal_wrapped(|ui| {
                let label = ui.label("New reaction equation");
                let input = ui.text_edit_singleline(&mut self.new_reaction_equation);
                let _ = input.clone().labelled_by(label.id);
            });
            ui.horizontal_wrapped(|ui| {
                for (parameter, name) in
                    self.new_reaction_arrhenius
                        .iter_mut()
                        .zip(["A", "n", "E (J/mol)"])
                {
                    let label = ui.label(name);
                    let input = ui.add(DragValue::new(parameter).speed(0.1));
                    let _ = input.labelled_by(label.id);
                }
                let label = ui.label("Heat effect (J/mol)");
                let input = ui.add(DragValue::new(&mut self.new_reaction_heat).speed(100.0));
                let _ = input.labelled_by(label.id);
                if ui.button("Add reaction").clicked() {
                    let equation = self.new_reaction_equation.trim().to_string();
                    if equation.is_empty() {
                        self.last_error = Some("Reaction equation cannot be empty.".to_string());
                    } else if self
                        .config
                        .reactions
                        .iter()
                        .any(|reaction| reaction.eq == equation)
                    {
                        self.last_error = Some(format!("Reaction `{equation}` already exists."));
                    } else {
                        match ReactionData::try_new_elementary(
                            equation.clone(),
                            self.new_reaction_arrhenius.to_vec(),
                            None,
                        ) {
                            Ok(reaction) => {
                                self.config.reactions.push(reaction);
                                self.config.thermal_effects.push(self.new_reaction_heat);
                                self.new_reaction_equation.clear();
                                self.last_status = Some(format!("Added reaction `{equation}`."));
                                self.last_error = None;
                            }
                            Err(error) => self.last_error = Some(error.to_string()),
                        }
                    }
                }
            });
        });
    }

    fn render_initial_state_section(&mut self, ui: &mut egui::Ui) {
        ui.heading("Initial state");
        Grid::new("ivp_initial_state_grid")
            .num_columns(4)
            .spacing([16.0, 8.0])
            .show(ui, |ui| {
                let label = ui.label("Temperature");
                let input =
                    ui.add(DragValue::new(&mut self.config.initial_state.temperature).speed(1.0));
                let _ = input.labelled_by(label.id);
                let label = ui.label("Heat flux");
                let input =
                    ui.add(DragValue::new(&mut self.config.initial_state.heat_flux).speed(0.01));
                let _ = input.labelled_by(label.id);
                ui.end_row();

                ui.label("Species");
                ui.label(format!("{}", self.config.substances.len()));
                ui.label("State order");
                ui.label(self.config.substances.join(", "));
                ui.end_row();
            });

        if self.config.substances.is_empty() {
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

                for species in &self.config.substances {
                    let value = self
                        .config
                        .initial_state
                        .species
                        .entry(species.clone())
                        .or_insert(0.0);
                    let label = ui.label(species);
                    let input = ui.add(DragValue::new(value).speed(0.01));
                    let _ = input.labelled_by(label.id);
                    ui.label("yes");
                    ui.end_row();
                }
            });
    }

    fn render_solver_backend_section(&mut self, ui: &mut egui::Ui) {
        ui.heading("Solver backend");

        ui.strong("Solver route");
        Grid::new("ivp_solver_route_grid")
            .num_columns(4)
            .spacing([16.0, 8.0])
            .show(ui, |ui| {
                let family_label = ui.label("Solver family").on_hover_text(
                    "LSODE2 is the legacy-compatible controller route; standalone BDF and Radau are separate solver families.",
                );
                let family_combo = ComboBox::from_id_salt("ivp_solver_family")
                    .selected_text(match self.config.solver_family {
                        ReactorIvpSolverFamily::Lsode2 => "LSODE2",
                        ReactorIvpSolverFamily::Bdf => "BDF (standalone)",
                        ReactorIvpSolverFamily::Radau => "Radau",
                    })
                    .show_ui(ui, |ui| {
                        ui.selectable_value(
                            &mut self.config.solver_family,
                            ReactorIvpSolverFamily::Lsode2,
                            "LSODE2",
                        );
                        ui.selectable_value(
                            &mut self.config.solver_family,
                            ReactorIvpSolverFamily::Bdf,
                            "BDF (standalone)",
                        );
                        ui.selectable_value(
                            &mut self.config.solver_family,
                            ReactorIvpSolverFamily::Radau,
                            "Radau",
                        );
                    });
                let _ = family_combo.response.clone().labelled_by(family_label.id);
                handle_combo_keyboard(
                    ui.ctx(),
                    &family_combo.response,
                    &mut self.config.solver_family,
                    &[
                        ReactorIvpSolverFamily::Lsode2,
                        ReactorIvpSolverFamily::Bdf,
                        ReactorIvpSolverFamily::Radau,
                    ],
                );
                if self.config.solver_family != ReactorIvpSolverFamily::Lsode2 {
                    self.config.matrix_backend = ReactorIvpMatrixBackend::Dense;
                    self.config.execution_mode = ReactorIvpGuiExecutionMode::Lambdify;
                    self.config.stop_condition_enabled = false;
                }
                let symbolic_label = ui.label("Symbolic frontend").on_hover_text(
                    "Select the symbolic expression frontend used to prepare residual and Jacobian callbacks.",
                );
                let symbolic_combo = ComboBox::from_id_salt("ivp_symbolic_backend")
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
                let _ = symbolic_combo.response.clone().labelled_by(symbolic_label.id);
                handle_combo_keyboard(
                    ui.ctx(),
                    &symbolic_combo.response,
                    &mut self.config.symbolic_backend,
                    &[
                        ReactorIvpSymbolicBackend::AtomView,
                        ReactorIvpSymbolicBackend::ExprLegacy,
                    ],
                );
                ui.end_row();
            });

        ui.strong("Execution");
        Grid::new("ivp_solver_execution_grid")
            .num_columns(4)
            .spacing([16.0, 8.0])
            .show(ui, |ui| {
                if self.config.solver_family == ReactorIvpSolverFamily::Lsode2 {
                    let method_label = ui.label("LSODE2 controller");
                    let method_combo = ComboBox::from_id_salt("ivp_method")
                        .selected_text(format!("{:?}", self.config.method))
                        .show_ui(ui, |ui| {
                            ui.selectable_value(
                                &mut self.config.method,
                                ReactorIvpMethod::Auto,
                                "Auto (Adams/BDF)",
                            );
                            ui.selectable_value(
                                &mut self.config.method,
                                ReactorIvpMethod::Bdf,
                                "BDF only",
                            );
                            ui.selectable_value(
                                &mut self.config.method,
                                ReactorIvpMethod::Adams,
                                "Adams only",
                            );
                        });
                    let _ = method_combo.response.clone().labelled_by(method_label.id);
                    handle_combo_keyboard(
                        ui.ctx(),
                        &method_combo.response,
                        &mut self.config.method,
                        &[ReactorIvpMethod::Auto, ReactorIvpMethod::Bdf, ReactorIvpMethod::Adams],
                    );
                    let matrix_label = ui.label("Matrix backend").on_hover_text(
                        "LSODE2 supports sparse and banded layouts. Standalone BDF currently requires Dense.",
                    );
                    let matrix_combo = ComboBox::from_id_salt("ivp_matrix_backend")
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
                    let _ = matrix_combo.response.clone().labelled_by(matrix_label.id);
                    handle_combo_keyboard(
                        ui.ctx(),
                        &matrix_combo.response,
                        &mut self.config.matrix_backend,
                        &[
                            ReactorIvpMatrixBackend::Sparse,
                            ReactorIvpMatrixBackend::Banded,
                        ],
                    );
                    ui.end_row();

                    let execution_label = ui.label("Execution backend").on_hover_text(
                        "AOT is currently available only for LSODE2 and requires explicit confirmation before invoking a toolchain.",
                    );
                    let execution_combo = ComboBox::from_id_salt("ivp_execution_mode")
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
                    let _ = execution_combo.response.clone().labelled_by(execution_label.id);
                    handle_combo_keyboard(
                        ui.ctx(),
                        &execution_combo.response,
                        &mut self.config.execution_mode,
                        &[
                            ReactorIvpGuiExecutionMode::Lambdify,
                            ReactorIvpGuiExecutionMode::Aot,
                        ],
                    );
                    ui.label("AOT / stop controls are LSODE2-only");
                    ui.end_row();
                } else {
                    ui.label("Matrix backend").on_hover_text(
                        "Standalone BDF and Radau currently use the Dense matrix layout through this GUI.",
                    );
                    ui.label("Dense (required by selected family)");
                    ui.label("Execution backend");
                    ui.label("Lambdify (supported route)");
                    ui.end_row();
                }

                match self.config.solver_family {
                    ReactorIvpSolverFamily::Bdf => {
                        let label = ui.label("BDF maximum steps").on_hover_text(
                            "Maximum accepted integration steps for standalone BDF.",
                        );
                        let input = ui.add(
                            DragValue::new(&mut self.config.bdf_max_steps)
                                .range(1..=usize::MAX),
                        );
                        let _ = input.labelled_by(label.id);
                        ui.end_row();
                    }
                    ReactorIvpSolverFamily::Radau => {
                        let label = ui.label("Radau maximum steps").on_hover_text(
                            "Maximum accepted integration steps for standalone Radau.",
                        );
                        let input = ui.add(
                            DragValue::new(&mut self.config.radau_max_steps)
                                .range(1..=usize::MAX),
                        );
                        let _ = input.labelled_by(label.id);
                        let label = ui.label("Radau maximum Newton iterations").on_hover_text(
                            "Maximum Newton iterations per standalone Radau step.",
                        );
                        let input = ui.add(
                            DragValue::new(&mut self.config.radau_max_newton_iterations)
                                .range(1..=usize::MAX),
                        );
                        let _ = input.labelled_by(label.id);
                        ui.end_row();
                    }
                    ReactorIvpSolverFamily::Lsode2 => {}
                }
            });

        if self.config.execution_mode == ReactorIvpGuiExecutionMode::Lambdify {
            ui.label("Callback execution");
            Grid::new("ivp_callback_execution_grid")
                .num_columns(4)
                .spacing([16.0, 8.0])
                .show(ui, |ui| {
                    let label = ui.label("Callback policy").on_hover_text(
                        "Controls sequential or parallel evaluation of symbolic residual/Jacobian callbacks. Sequential is the conservative default.",
                    );
                    let combo = ComboBox::from_id_salt("ivp_callback_policy")
                        .selected_text(match self.config.callback_policy {
                            ReactorIvpGuiCallbackPolicy::Sequential => "Sequential",
                            ReactorIvpGuiCallbackPolicy::Parallel => "Parallel",
                            ReactorIvpGuiCallbackPolicy::Auto => "Auto",
                        })
                        .show_ui(ui, |ui| {
                            ui.selectable_value(
                                &mut self.config.callback_policy,
                                ReactorIvpGuiCallbackPolicy::Sequential,
                                "Sequential",
                            );
                            ui.selectable_value(
                                &mut self.config.callback_policy,
                                ReactorIvpGuiCallbackPolicy::Parallel,
                                "Parallel",
                            );
                            ui.selectable_value(
                                &mut self.config.callback_policy,
                                ReactorIvpGuiCallbackPolicy::Auto,
                                "Auto",
                            );
                        });
                    let _ = combo.response.clone().labelled_by(label.id);
                    handle_combo_keyboard(
                        ui.ctx(),
                        &combo.response,
                        &mut self.config.callback_policy,
                        &[
                            ReactorIvpGuiCallbackPolicy::Sequential,
                            ReactorIvpGuiCallbackPolicy::Parallel,
                            ReactorIvpGuiCallbackPolicy::Auto,
                        ],
                    );
                    if self.config.callback_policy != ReactorIvpGuiCallbackPolicy::Sequential {
                        let label = ui.label("Minimum callback work").on_hover_text(
                            "Work threshold used by Parallel or Auto callback scheduling.",
                        );
                        let input = ui.add(
                            DragValue::new(&mut self.config.callback_min_work)
                                .range(0..=usize::MAX),
                        );
                        let _ = input.labelled_by(label.id);
                    } else {
                        ui.label("Threshold unused");
                    }
                    ui.end_row();
                });
        }

        ui.separator();
        ui.strong("Diagnostics");
        Grid::new("ivp_solver_diagnostics_grid")
            .num_columns(4)
            .spacing([16.0, 8.0])
            .show(ui, |ui| {
                let telemetry_label = ui.label("Telemetry").on_hover_text(
                    "Telemetry collection is independent from solver logging; Off avoids collecting counters and timings.",
                );
                let telemetry_combo = ComboBox::from_id_salt("ivp_telemetry_mode")
                    .selected_text(format!("{:?}", self.config.telemetry))
                    .show_ui(ui, |ui| {
                        ui.selectable_value(
                            &mut self.config.telemetry,
                            SolverTelemetryMode::Off,
                            "Off",
                        );
                        ui.selectable_value(
                            &mut self.config.telemetry,
                            SolverTelemetryMode::Counters,
                            "Counters",
                        );
                        ui.selectable_value(
                            &mut self.config.telemetry,
                            SolverTelemetryMode::Timings,
                            "Timings",
                        );
                        ui.selectable_value(
                            &mut self.config.telemetry,
                            SolverTelemetryMode::Detailed,
                            "Detailed",
                        );
                    });
                let _ = telemetry_combo.response.clone().labelled_by(telemetry_label.id);
                handle_combo_keyboard(
                    ui.ctx(),
                    &telemetry_combo.response,
                    &mut self.config.telemetry,
                    &[
                        SolverTelemetryMode::Off,
                        SolverTelemetryMode::Counters,
                        SolverTelemetryMode::Timings,
                        SolverTelemetryMode::Detailed,
                    ],
                );
                let logging_label = ui.label("Solver logging").on_hover_text(
                    "Logging controls solver messages only and does not enable telemetry collection.",
                );
                let logging_combo = ComboBox::from_id_salt("ivp_logging_level")
                    .selected_text(format!("{:?}", self.config.logging))
                    .show_ui(ui, |ui| {
                        ui.selectable_value(&mut self.config.logging, SolverLogLevel::Off, "Off");
                        ui.selectable_value(
                            &mut self.config.logging,
                            SolverLogLevel::Error,
                            "Error",
                        );
                        ui.selectable_value(&mut self.config.logging, SolverLogLevel::Warn, "Warn");
                        ui.selectable_value(&mut self.config.logging, SolverLogLevel::Info, "Info");
                        ui.selectable_value(
                            &mut self.config.logging,
                            SolverLogLevel::Debug,
                            "Debug",
                        );
                    });
                let _ = logging_combo.response.clone().labelled_by(logging_label.id);
                handle_combo_keyboard(
                    ui.ctx(),
                    &logging_combo.response,
                    &mut self.config.logging,
                    &[
                        SolverLogLevel::Off,
                        SolverLogLevel::Error,
                        SolverLogLevel::Warn,
                        SolverLogLevel::Info,
                        SolverLogLevel::Debug,
                    ],
                );
                ui.end_row();
            });

        if self.config.solver_family == ReactorIvpSolverFamily::Lsode2
            && matches!(self.config.execution_mode, ReactorIvpGuiExecutionMode::Aot)
        {
            ui.separator();
            ui.strong("AOT build");
            ui.label("Compiler and build profile are only relevant for AOT.");
            Grid::new("ivp_aot_settings_grid")
                .num_columns(4)
                .spacing([16.0, 8.0])
                .show(ui, |ui| {
                    let toolchain_label = ui.label("AOT toolchain");
                    let toolchain_combo = ComboBox::from_id_salt("ivp_aot_toolchain")
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
                    let _ = toolchain_combo
                        .response
                        .clone()
                        .labelled_by(toolchain_label.id);
                    handle_combo_keyboard(
                        ui.ctx(),
                        &toolchain_combo.response,
                        &mut self.config.aot_toolchain,
                        &[
                            Lsode2AotToolchain::CTcc,
                            Lsode2AotToolchain::CGcc,
                            Lsode2AotToolchain::Zig,
                            Lsode2AotToolchain::Rust,
                        ],
                    );
                    let profile_label = ui.label("AOT profile");
                    let profile_combo = ComboBox::from_id_salt("ivp_aot_profile")
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
                    let _ = profile_combo.response.clone().labelled_by(profile_label.id);
                    handle_combo_keyboard(
                        ui.ctx(),
                        &profile_combo.response,
                        &mut self.config.aot_profile,
                        &[Lsode2AotProfile::Release, Lsode2AotProfile::Debug],
                    );
                    ui.end_row();
                    ui.label("AOT confirmation");
                    ui.label("AOT runs must be explicitly confirmed before execution.");
                    ui.end_row();
                });
        }
    }

    fn render_integration_section(&mut self, ui: &mut egui::Ui) {
        ui.heading("Integration");
        Grid::new("ivp_integration_grid")
            .num_columns(4)
            .spacing([16.0, 8.0])
            .show(ui, |ui| {
                let label = ui.label("x0");
                let input = ui.add(DragValue::new(&mut self.config.x0).speed(0.01));
                let _ = input.labelled_by(label.id);
                let label = ui.label("x_bound");
                let input = ui.add(DragValue::new(&mut self.config.x_bound).speed(0.01));
                let _ = input.labelled_by(label.id);
                ui.end_row();

                let first_step_label = ui.label("first_step");
                let mut first_step_enabled = self.config.first_step.is_some();
                let enabled_response = ui.checkbox(&mut first_step_enabled, "enabled");
                let _ = enabled_response.labelled_by(first_step_label.id);
                if first_step_enabled {
                    let value = self.config.first_step.get_or_insert(1e-6);
                    let input = ui.add(DragValue::new(value).speed(0.000_001));
                    let _ = input.labelled_by(first_step_label.id);
                } else {
                    self.config.first_step = None;
                    ui.label("None");
                }
                ui.end_row();

                let label = ui.label("max_step");
                let input = ui.add(DragValue::new(&mut self.config.max_step).speed(0.01));
                let _ = input.labelled_by(label.id);
                let label = ui.label("rtol");
                let input = ui.add(DragValue::new(&mut self.config.rtol).speed(0.000001));
                let _ = input.labelled_by(label.id);
                ui.end_row();

                let label = ui.label("atol");
                let input = ui.add(DragValue::new(&mut self.config.atol).speed(0.0000001));
                let _ = input.labelled_by(label.id);
                ui.end_row();
            });
    }

    fn render_postprocessing_section(&mut self, ui: &mut egui::Ui) {
        ui.heading("Postprocessing");
        Grid::new("ivp_postprocessing_grid")
            .num_columns(4)
            .spacing([16.0, 8.0])
            .show(ui, |ui| {
                ui.checkbox(&mut self.config.postprocessing.save_txt, "Save TXT");
                let label = ui.label("TXT path");
                let input = ui.add(
                    TextEdit::singleline(&mut self.config.postprocessing.txt_path)
                        .desired_width(360.0)
                        .hint_text("auto"),
                );
                let _ = input.labelled_by(label.id);
                ui.end_row();

                ui.checkbox(&mut self.config.postprocessing.save_csv, "Save CSV");
                let label = ui.label("CSV path");
                let input = ui.add(
                    TextEdit::singleline(&mut self.config.postprocessing.csv_path)
                        .desired_width(360.0)
                        .hint_text("auto"),
                );
                let _ = input.labelled_by(label.id);
                ui.end_row();

                ui.checkbox(&mut self.config.postprocessing.write_report, "Write report");
                let label = ui.label("Report path");
                let input = ui.add(
                    TextEdit::singleline(&mut self.config.postprocessing.report_path)
                        .desired_width(360.0)
                        .hint_text("auto"),
                );
                let _ = input.labelled_by(label.id);
                ui.end_row();

                ui.checkbox(&mut self.config.postprocessing.plotters_png, "Plotters PNG");
                let label = ui.label("Plotters dir");
                let input = ui.add(
                    TextEdit::singleline(&mut self.config.postprocessing.plotters_dir)
                        .desired_width(360.0)
                        .hint_text("auto"),
                );
                let _ = input.labelled_by(label.id);
                ui.end_row();

                ui.checkbox(&mut self.config.postprocessing.gnuplot_png, "Gnuplot PNG");
                let label = ui.label("Gnuplot dir");
                let input = ui.add(
                    TextEdit::singleline(&mut self.config.postprocessing.gnuplot_dir)
                        .desired_width(360.0)
                        .hint_text("auto"),
                );
                let _ = input.labelled_by(label.id);
                ui.end_row();

                ui.checkbox(
                    &mut self.config.postprocessing.terminal_plot,
                    "Terminal plot",
                );
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

    pub(crate) fn request_calculation(&mut self) {
        if self.calculation_state.is_active() {
            self.last_error = Some("A calculation is already running.".to_string());
            return;
        }

        self.invalidate_published_result();

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

        self.invalidate_published_result();

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

        self.invalidate_published_result();

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

    fn invalidate_published_result(&mut self) {
        self.plot_window = None;
        self.task.clear_published_solution();
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

        for substance in &config.substances {
            substance.hash(&mut hasher);
        }
        for reaction in &config.reactions {
            if let Ok(serialized) = serde_json::to_vec(reaction) {
                hasher.write(&serialized);
            } else {
                reaction.eq.hash(&mut hasher);
            }
        }

        for effect in &config.thermal_effects {
            hasher.write_u64(effect.to_bits());
        }

        config.solver_family.hash(&mut hasher);
        if config.solver_family == ReactorIvpSolverFamily::Lsode2 {
            config.method.hash(&mut hasher);
            format!("{:?}", config.native_execution).hash(&mut hasher);
        }
        config.symbolic_backend.hash(&mut hasher);
        config.matrix_backend.hash(&mut hasher);
        config.execution_mode.hash(&mut hasher);
        config.callback_policy.hash(&mut hasher);
        hasher.write_usize(config.callback_min_work);
        hasher.write_usize(config.bdf_max_steps);
        hasher.write_usize(config.radau_max_steps);
        hasher.write_usize(config.radau_max_newton_iterations);
        format!("{:?}", config.telemetry).hash(&mut hasher);
        format!("{:?}", config.logging).hash(&mut hasher);
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
                "H2=>2H".to_string(),
                vec![1.0e3, 0.0, 1.0e4],
                None,
            )],
            None,
        )
        .map_err(|err| IvpError::InvalidConfiguration(err.to_string()))?;
    task.kindata.stecheodata.vec_of_molmasses = Some(vec![2.016, 1.008]);
    task.kindata.substances = vec!["H2".to_string(), "H".to_string()];
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
        ("H2".to_string(), 0.7),
        ("H".to_string(), 0.3),
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
        config.telemetry = SolverTelemetryMode::Detailed;
        config.logging = SolverLogLevel::Info;

        let solver = config.to_solver_config();
        assert_eq!(solver.family(), ReactorIvpSolverFamily::Lsode2);
        assert_eq!(solver.telemetry, SolverTelemetryMode::Detailed);
        assert_eq!(solver.logging, SolverLogLevel::Info);
        assert_eq!(solver.lsode2_controller(), Some(ReactorIvpMethod::Adams));
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

    #[test]
    fn gui_solver_family_is_preserved_in_backend_config() {
        let mut config = IvpGuiConfig::default();
        config.solver_family = ReactorIvpSolverFamily::Radau;

        let solver = config.to_solver_config();
        assert_eq!(solver.family(), ReactorIvpSolverFamily::Radau);
    }

    #[test]
    fn standalone_gui_fingerprint_ignores_inactive_lsode2_controls() {
        let mut first = IvpGuiConfig::default();
        first.solver_family = ReactorIvpSolverFamily::Bdf;
        let mut second = first.clone();
        second.method = ReactorIvpMethod::Adams;
        second.native_execution = Lsode2NativeExecutionConfig::faithful_bdf_solve(17, 19);

        assert_eq!(
            ReactorIvpApp::fingerprint_config(&first),
            ReactorIvpApp::fingerprint_config(&second)
        );
    }
}
