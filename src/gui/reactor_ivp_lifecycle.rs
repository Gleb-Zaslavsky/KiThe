use std::fs;
use std::path::PathBuf;

use RustedSciThe::command_interpreter::task_parser_ivp::{
    IvpMethodSpec, IvpSolverOptionsSpec, IvpSolverSettingsSpec, Lsode2TaskExecutionSpec,
    Lsode2TaskOptionsSpec, SolverSelectionSpec, TaskKindSpec,
};
use RustedSciThe::numerical::LSODE2::config::{
    Lsode2LinearSystemStructure, Lsode2SymbolicAssemblyBackend,
};
use egui::Ui;
use log::{info, warn};
use rfd::FileDialog;

use crate::gui::document_lifecycle::{GuiFileOperationKind, GuiFileOperationResult};

use super::{IvpGuiConfig, ReactorIvpApp, ReactorIvpCalculationState};
use crate::ReactorsIVP::SimpleReactorIVP::IvpError;
use crate::ReactorsIVP::solver_backend::{
    ReactorIvpExecutionBackend, ReactorIvpMatrixBackend, ReactorIvpMethod, ReactorIvpSolverConfig,
    ReactorIvpSymbolicBackend,
};
use crate::ReactorsIVP::task_parser_reactor_IVP::{
    ReactorIvpDocumentSpec, ReactorIvpStopConditionSpec,
};

impl ReactorIvpApp {
    /// Rebuild the canonical document snapshot for the current editor state.
    fn current_document_spec(&self) -> Result<ReactorIvpDocumentSpec, IvpError> {
        let mut preview_task = self.task.clone();
        self.config.apply_to_task(&mut preview_task)?;

        let stop_condition = match preview_task.solver_backend_config.stop_condition {
            Some(stop_condition) => {
                let species = preview_task
                    .kindata
                    .substances
                    .get(stop_condition.species_index)
                    .cloned()
                    .ok_or_else(|| {
                        IvpError::InvalidConfiguration(format!(
                            "stop condition species index {} is out of range",
                            stop_condition.species_index
                        ))
                    })?;
                ReactorIvpStopConditionSpec {
                    enabled: true,
                    species: Some(species),
                    threshold: Some(stop_condition.threshold),
                }
            }
            None => ReactorIvpStopConditionSpec::default(),
        };

        Ok(ReactorIvpDocumentSpec {
            problem_name: preview_task.problem_name.clone(),
            problem_description: preview_task.problem_description.clone(),
            physical: crate::ReactorsIVP::SimpleReactorIVP::ReactorIvpPhysicalConfig {
                ro: preview_task.ro,
                Cp: preview_task.Cp,
                Lambda: preview_task.Lambda,
                m: preview_task.m,
                L: preview_task.L,
            },
            scaling: preview_task.scaling.clone(),
            x0: preview_task.solver_backend_config.x0,
            x_bound: preview_task.solver_backend_config.x_bound,
            initial_conditions: preview_task.initial_conditions.clone(),
            thermal_effects: preview_task.thermal_effects.clone(),
            stop_condition,
            solver_settings: reactor_solver_config_to_parser_settings(
                &preview_task.solver_backend_config,
            ),
        })
    }

    /// Return the canonical serialized document for lifecycle fingerprinting.
    pub(crate) fn current_document_fingerprint(&self) -> Result<String, IvpError> {
        Ok(self.current_document_spec()?.document_to_string())
    }

    /// Synchronize the dirty flag against the current serialized snapshot.
    pub(crate) fn refresh_document_lifecycle_dirty_state(&mut self) {
        match self.current_document_fingerprint() {
            Ok(fingerprint) => self.document_lifecycle.sync_dirty_state(fingerprint),
            Err(err) => {
                warn!("failed to fingerprint reactor IVP document: {err}");
                self.document_lifecycle.mark_dirty();
            }
        }
    }

    /// Update the current document path in both the compatibility field and lifecycle state.
    pub(crate) fn set_document_path(&mut self, path: Option<PathBuf>) {
        self.current_file_path = path.clone();
        self.document_lifecycle.current_path = path;
    }

    /// Render the document lifecycle toolbar.
    pub(crate) fn render_document_toolbar(&mut self, ui: &mut Ui) {
        self.refresh_document_lifecycle_dirty_state();

        ui.horizontal_wrapped(|ui| {
            ui.label("Document");
            ui.monospace(
                self.current_file_path
                    .as_ref()
                    .map(|path| path.display().to_string())
                    .unwrap_or_else(|| "unsaved".to_string()),
            );
            if self.document_lifecycle.dirty {
                ui.colored_label(egui::Color32::from_rgb(185, 110, 35), "unsaved changes");
            } else {
                ui.label("clean");
            }

            if ui.button("Load Document").clicked() {
                if let Some(path) = FileDialog::new().pick_file() {
                    let _ = self.request_document_load(path);
                }
            }

            if ui.button("Save").clicked() {
                if let Some(path) = self.current_file_path.clone() {
                    let _ = self.save_document(path);
                } else if let Some(path) = FileDialog::new().save_file() {
                    let _ = self.save_document(path);
                }
            }

            if ui.button("Save As").clicked() {
                if let Some(path) = FileDialog::new().save_file() {
                    let _ = self.save_document(path);
                }
            }
        });
    }

    /// Persist the current typed document to disk.
    pub fn save_document(&mut self, path: PathBuf) -> GuiFileOperationResult {
        let result = match self.current_document_spec() {
            Ok(spec) => {
                let serialized = spec.document_to_string();
                match fs::write(&path, &serialized) {
                    Ok(_) => {
                        self.set_document_path(Some(path.clone()));
                        self.document_lifecycle.record_clean(
                            GuiFileOperationResult::Saved { path: path.clone() },
                            serialized,
                        );
                        self.pending_document_load = None;
                        self.pending_window_close = false;
                        self.last_status =
                            Some(format!("Saved IVP document to {}.", path.display()));
                        self.last_error = None;
                        info!("Successfully saved reactor IVP document to {:?}", path);
                        GuiFileOperationResult::Saved { path: path.clone() }
                    }
                    Err(error) => {
                        let message = error.to_string();
                        warn!("Error saving reactor IVP document: {}", message);
                        GuiFileOperationResult::Failed {
                            kind: GuiFileOperationKind::Save,
                            path: Some(path.clone()),
                            message,
                        }
                    }
                }
            }
            Err(error) => GuiFileOperationResult::Failed {
                kind: GuiFileOperationKind::Save,
                path: Some(path.clone()),
                message: error.to_string(),
            },
        };

        if result.is_error() {
            self.document_lifecycle.record(result.clone());
            if let GuiFileOperationResult::Failed { message, .. } = &result {
                self.last_error = Some(message.clone());
                self.last_status = Some("Save failed.".to_string());
            }
        }

        result
    }

    /// Load a typed document from disk and resynchronize the GUI state.
    pub(crate) fn load_document_from_path(&mut self, path: PathBuf) -> GuiFileOperationResult {
        let result = match fs::read_to_string(&path) {
            Ok(content) => match crate::ReactorsIVP::task_parser_reactor_IVP::parse_reactor_ivp_document_from_str(&content) {
                Ok(spec) => {
                    let mut task = self.task.clone();
                    if let Err(error) = spec.validate(&task.kindata.substances) {
                        GuiFileOperationResult::Failed {
                            kind: GuiFileOperationKind::Load,
                            path: Some(path.clone()),
                            message: error.to_string(),
                        }
                    } else if let Err(error) = spec.apply_to_task(&mut task) {
                        GuiFileOperationResult::Failed {
                            kind: GuiFileOperationKind::Load,
                            path: Some(path.clone()),
                            message: error.to_string(),
                        }
                    } else {
                        let mut solver_config =
                            parser_solver_settings_to_reactor_config(&spec.solver_settings);
                        solver_config =
                            solver_config.with_integration_domain(spec.x0, spec.x_bound);
                        if spec.stop_condition.enabled {
                            let Some(species_name) = spec.stop_condition.species.as_ref() else {
                                return GuiFileOperationResult::Failed {
                                    kind: GuiFileOperationKind::Load,
                                    path: Some(path.clone()),
                                    message:
                                        "Missing `species` in `stop_condition` section"
                                            .to_string(),
                                };
                            };
                            let Some(species_index) = task
                                .kindata
                                .substances
                                .iter()
                                .position(|name| name == species_name)
                            else {
                                return GuiFileOperationResult::Failed {
                                    kind: GuiFileOperationKind::Load,
                                    path: Some(path.clone()),
                                    message: format!(
                                        "Unknown stop-condition species `{species_name}`"
                                    ),
                                };
                            };
                            let Some(threshold) = spec.stop_condition.threshold else {
                                return GuiFileOperationResult::Failed {
                                    kind: GuiFileOperationKind::Load,
                                    path: Some(path.clone()),
                                    message:
                                        "Missing `threshold` in `stop_condition` section"
                                            .to_string(),
                                };
                            };
                            solver_config =
                                solver_config.with_stop_condition_le(species_index, threshold);
                        } else {
                            solver_config = solver_config.without_stop_condition();
                        }
                        task.set_solver_backend_config(solver_config);

                        self.task = task;
                        self.config = IvpGuiConfig::from_task(&self.task);
                        self.set_document_path(Some(path.clone()));
                        self.document_lifecycle.record_clean(
                            GuiFileOperationResult::Loaded { path: path.clone() },
                            spec.document_to_string(),
                        );
                        self.pending_document_load = None;
                        self.pending_window_close = false;
                        self.plot_window = None;
                        self.last_preview = None;
                        self.calculation_state = ReactorIvpCalculationState::Idle;
                        self.calculation_receiver = None;
                        self.last_status =
                            Some(format!("Loaded IVP document from {}.", path.display()));
                        self.last_error = None;
                        info!("Successfully loaded reactor IVP document from {:?}", path);
                        GuiFileOperationResult::Loaded { path: path.clone() }
                    }
                }
                Err(error) => GuiFileOperationResult::Failed {
                    kind: GuiFileOperationKind::Load,
                    path: Some(path.clone()),
                    message: error.to_string(),
                },
            },
            Err(error) => GuiFileOperationResult::Failed {
                kind: GuiFileOperationKind::Load,
                path: Some(path.clone()),
                message: error.to_string(),
            },
        };

        if result.is_error() {
            self.document_lifecycle.record(result.clone());
            if let GuiFileOperationResult::Failed { message, .. } = &result {
                warn!("Error loading reactor IVP document: {}", message);
                self.last_error = Some(message.clone());
                self.last_status = Some("Load failed.".to_string());
            }
        }

        result
    }

    /// Request loading a new document and confirm before replacing unsaved work.
    pub(crate) fn request_document_load(&mut self, path: PathBuf) -> bool {
        self.refresh_document_lifecycle_dirty_state();
        if self.document_lifecycle.dirty || self.calculation_state.is_active() {
            self.pending_document_load = Some(path.clone());
            self.last_status = Some(format!(
                "Confirm loading {} before discarding the current task state.",
                path.display()
            ));
            false
        } else {
            !self.load_document_from_path(path).is_error()
        }
    }

    /// Confirm a queued document load.
    pub(crate) fn confirm_pending_document_load(&mut self) -> bool {
        let Some(path) = self.pending_document_load.take() else {
            return false;
        };
        !self.load_document_from_path(path).is_error()
    }

    /// Cancel a queued document load without touching the active task.
    pub(crate) fn cancel_pending_document_load(&mut self) {
        if let Some(path) = self.pending_document_load.take() {
            self.last_status = Some(format!(
                "Load cancelled; the current task was preserved instead of loading {}.",
                path.display()
            ));
        }
    }

    /// Request closing the main window and confirm before discarding unsaved work.
    pub(crate) fn request_window_close(&mut self) -> bool {
        self.refresh_document_lifecycle_dirty_state();
        if self.document_lifecycle.dirty || self.calculation_state.is_active() {
            self.pending_window_close = true;
            self.last_status = Some(
                "Confirm closing the window before discarding the current task state.".to_string(),
            );
            false
        } else {
            true
        }
    }

    /// Confirm the queued close request.
    pub(crate) fn confirm_pending_window_close(&mut self) -> bool {
        if !self.pending_window_close {
            return false;
        }
        self.pending_window_close = false;
        true
    }

    /// Cancel the queued close request without touching the task.
    pub(crate) fn cancel_pending_window_close(&mut self) {
        if self.pending_window_close {
            self.pending_window_close = false;
            self.last_status = Some("Close cancelled; the current task was preserved.".to_string());
        }
    }

    /// Render confirmation for loading a document while unsaved work exists.
    pub(crate) fn show_document_load_confirmation(&mut self, ctx: &egui::Context) {
        let Some(path) = self.pending_document_load.clone() else {
            return;
        };

        let mut confirm = false;
        egui::Window::new("Confirm document load")
            .collapsible(false)
            .resizable(false)
            .show(ctx, |ui| {
                ui.label(format!(
                    "Load {} and discard the current editor state?",
                    path.display()
                ));
                ui.horizontal(|ui| {
                    if ui.button("Cancel load").clicked() {
                        self.cancel_pending_document_load();
                    }
                    if ui.button("Load anyway").clicked() {
                        confirm = true;
                    }
                });
            });

        if confirm {
            let _ = self.confirm_pending_document_load();
        }
    }

    /// Render confirmation for closing the window with unsaved work.
    pub(crate) fn show_window_close_confirmation(&mut self, ctx: &egui::Context, open: &mut bool) {
        if !self.pending_window_close {
            return;
        }

        let mut confirm = false;
        egui::Window::new("Confirm window close")
            .collapsible(false)
            .resizable(false)
            .show(ctx, |ui| {
                ui.label("Close the window and discard the current editor state?");
                ui.horizontal(|ui| {
                    if ui.button("Cancel close").clicked() {
                        self.cancel_pending_window_close();
                    }
                    if ui.button("Close anyway").clicked() {
                        confirm = true;
                    }
                });
            });

        if confirm {
            if self.confirm_pending_window_close() {
                *open = false;
            }
        }
    }
}

fn reactor_solver_config_to_parser_settings(
    config: &ReactorIvpSolverConfig,
) -> IvpSolverSettingsSpec {
    let execution_backend = match config.execution_backend {
        ReactorIvpExecutionBackend::Lambdify => Lsode2TaskExecutionSpec::LambdifyExpr,
        ReactorIvpExecutionBackend::Aot { toolchain, profile } => Lsode2TaskExecutionSpec::Aot {
            toolchain,
            profile,
            output_parent_dir: None,
        },
    };

    IvpSolverSettingsSpec {
        solver: SolverSelectionSpec {
            task_kind: TaskKindSpec::Ivp,
            method: IvpMethodSpec::Lsode2,
        },
        solver_options: IvpSolverOptionsSpec {
            step_size: None,
            tolerance: None,
            max_iterations: None,
            rtol: Some(config.rtol),
            atol: Some(config.atol),
            max_step: Some(config.max_step),
            first_step: config.first_step,
            vectorized: None,
            parallel: None,
            neighborhood_check: None,
            lsode2: Some(Lsode2TaskOptionsSpec {
                controller: Some(config.method.to_controller()),
                symbolic_assembly: Some(config.symbolic_backend.to_rusted()),
                symbolic_execution: Some(execution_backend),
                linear_system_structure: Some(match config.matrix_backend {
                    ReactorIvpMatrixBackend::Sparse => Lsode2LinearSystemStructure::Sparse,
                    ReactorIvpMatrixBackend::Banded => {
                        Lsode2LinearSystemStructure::Banded { kl: 0, ku: 0 }
                    }
                }),
                linear_solver_policy: Some(config.to_rusted_linear_policy()),
                native_execution: Some(config.native_execution),
            }),
        },
    }
}

fn parser_solver_settings_to_reactor_config(
    settings: &IvpSolverSettingsSpec,
) -> ReactorIvpSolverConfig {
    let mut config = ReactorIvpSolverConfig::default();

    config.method = match settings
        .solver_options
        .lsode2
        .as_ref()
        .and_then(|lsode2| lsode2.controller)
        .map(|controller| controller.mode)
    {
        Some(mode) => match mode {
            RustedSciThe::numerical::LSODE2::algorithm::Lsode2ControllerMode::AutomaticAdamsBdf => {
                ReactorIvpMethod::Auto
            }
            RustedSciThe::numerical::LSODE2::algorithm::Lsode2ControllerMode::BdfOnly => {
                ReactorIvpMethod::Bdf
            }
            RustedSciThe::numerical::LSODE2::algorithm::Lsode2ControllerMode::AdamsOnly => {
                ReactorIvpMethod::Adams
            }
        },
        None => ReactorIvpMethod::default(),
    };

    config.symbolic_backend = settings
        .solver_options
        .lsode2
        .as_ref()
        .and_then(|lsode2| lsode2.symbolic_assembly)
        .map(|backend| match backend {
            Lsode2SymbolicAssemblyBackend::AtomView => ReactorIvpSymbolicBackend::AtomView,
            Lsode2SymbolicAssemblyBackend::ExprLegacy => ReactorIvpSymbolicBackend::ExprLegacy,
        })
        .unwrap_or_default();

    if let Some(lsode2) = settings.solver_options.lsode2.as_ref() {
        if let Some(execution) = &lsode2.symbolic_execution {
            match execution {
                Lsode2TaskExecutionSpec::LambdifyExpr => {
                    config.execution_backend = ReactorIvpExecutionBackend::Lambdify;
                }
                Lsode2TaskExecutionSpec::Aot {
                    toolchain, profile, ..
                } => {
                    config.execution_backend = ReactorIvpExecutionBackend::Aot {
                        toolchain: *toolchain,
                        profile: *profile,
                    };
                }
            }
        }

        if let Some(structure) = lsode2.linear_system_structure {
            config.matrix_backend = match structure {
                Lsode2LinearSystemStructure::Sparse => ReactorIvpMatrixBackend::Sparse,
                Lsode2LinearSystemStructure::Banded { .. } => ReactorIvpMatrixBackend::Banded,
                Lsode2LinearSystemStructure::Dense => ReactorIvpMatrixBackend::Sparse,
            };
        }

        if let Some(native_execution) = lsode2.native_execution {
            config.native_execution = native_execution;
        }
    }

    config.rtol = settings.solver_options.rtol.unwrap_or(config.rtol);
    config.atol = settings.solver_options.atol.unwrap_or(config.atol);
    config.max_step = settings.solver_options.max_step.unwrap_or(config.max_step);
    config.first_step = settings.solver_options.first_step;
    config
}
