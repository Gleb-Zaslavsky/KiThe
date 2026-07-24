use std::path::PathBuf;
use std::sync::atomic::Ordering;

use RustedSciThe::command_interpreter::task_parser::DocumentParser;
use log::{info, warn};

use super::{
    BvpGuiConfig, CombustionApp, GuiFileOperationKind, GuiFileOperationResult, ProblemsEnum,
};

impl CombustionApp {
    /// Saves the current document to the specified file path.
    ///
    /// The returned typed result is also stored in the explicit document
    /// lifecycle state so the GUI can surface it to the user.
    pub fn save_document(&mut self, path: PathBuf) -> GuiFileOperationResult {
        let content = self.document_to_string();
        let result = match std::fs::write(&path, content) {
            Ok(_) => GuiFileOperationResult::Saved { path: path.clone() },
            Err(e) => GuiFileOperationResult::Failed {
                kind: GuiFileOperationKind::Save,
                path: Some(path.clone()),
                message: e.to_string(),
            },
        };

        match &result {
            GuiFileOperationResult::Saved { path } => {
                info!("Successfully saved to: {:?}", path);
            }
            GuiFileOperationResult::Failed { message, .. } => {
                warn!("Error saving file: {}", message);
            }
            GuiFileOperationResult::Loaded { .. } => {}
        }

        if let GuiFileOperationResult::Saved { path } = &result {
            self.set_document_path(Some(path.clone()));
            self.document_lifecycle
                .record_clean(result.clone(), self.current_document_fingerprint());
        } else {
            self.document_lifecycle.record(result.clone());
        }
        result
    }

    /// Loads a document from disk, parses it, and updates the active GUI state.
    pub(crate) fn load_document_from_path(&mut self, path: PathBuf) -> GuiFileOperationResult {
        let result = match std::fs::read_to_string(&path) {
            Ok(content) => {
                let mut parser = DocumentParser::new(content);
                match parser.parse_document() {
                    Ok(_) => match parser.get_result() {
                        Some(parsed_doc) => {
                            self.document = parsed_doc.clone();
                            if matches!(self.selected_problem, ProblemsEnum::BVPSimple) {
                                super::seed_bvp_document_defaults(&mut self.document);
                            }
                            let (bvp_gui_config, bvp_gui_migration_report) =
                                BvpGuiConfig::from_document(&self.document);
                            self.bvp_gui_config = bvp_gui_config;
                            self.bvp_gui_migration_report = bvp_gui_migration_report;
                            self.validation_fingerprint = None;
                            self.refresh_validation_report();
                            self.set_document_path(Some(path.clone()));
                            self.document_lifecycle.record_clean(
                                GuiFileOperationResult::Loaded { path: path.clone() },
                                self.current_document_fingerprint(),
                            );
                            GuiFileOperationResult::Loaded { path: path.clone() }
                        }
                        None => GuiFileOperationResult::Failed {
                            kind: GuiFileOperationKind::Load,
                            path: Some(path.clone()),
                            message: "Parser returned no result".to_string(),
                        },
                    },
                    Err(error) => GuiFileOperationResult::Failed {
                        kind: GuiFileOperationKind::Load,
                        path: Some(path.clone()),
                        message: error.to_string(),
                    },
                }
            }
            Err(error) => GuiFileOperationResult::Failed {
                kind: GuiFileOperationKind::Load,
                path: Some(path.clone()),
                message: error.to_string(),
            },
        };

        match &result {
            GuiFileOperationResult::Loaded { path } => {
                info!("Successfully loaded and parsed file: {:?}", path);
            }
            GuiFileOperationResult::Failed { message, .. } => {
                warn!("Error loading file: {}", message);
            }
            GuiFileOperationResult::Saved { .. } => {}
        }

        if matches!(result, GuiFileOperationResult::Failed { .. }) {
            self.document_lifecycle.record(result.clone());
        }
        result
    }

    /// Synchronize typed snapshot and lifecycle state after a document edit.
    ///
    /// The GUI mutates the live document from several UI paths, so this helper
    /// keeps the snapshot and dirty flag aligned immediately after each edit.
    pub(crate) fn sync_document_state_after_edit(&mut self) {
        self.refresh_bvp_gui_snapshot();
        self.refresh_document_lifecycle_dirty_state();
    }

    /// Updates the active file path in both the compatibility field and the
    /// explicit lifecycle state.
    fn set_document_path(&mut self, path: Option<PathBuf>) {
        self.current_file_path = path.clone();
        self.document_lifecycle.current_path = path;
    }

    /// Returns a stable fingerprint for the current document and active problem.
    fn current_document_fingerprint(&self) -> String {
        self.current_validation_fingerprint()
    }

    /// Recomputes the dirty flag against the latest clean snapshot.
    pub(crate) fn refresh_document_lifecycle_dirty_state(&mut self) {
        self.document_lifecycle
            .sync_dirty_state(self.current_document_fingerprint());
    }

    /// Returns true when switching problems should ask for confirmation.
    fn problem_switch_requires_confirmation(&mut self) -> bool {
        self.refresh_document_lifecycle_dirty_state();
        self.document_lifecycle.dirty || self.calculation_state.is_active()
    }

    /// Request loading a new document and confirm before replacing unsaved work.
    pub(crate) fn request_document_load(&mut self, path: PathBuf) -> bool {
        if self.problem_switch_requires_confirmation() {
            self.pending_document_load = Some(path.clone());
            self.last_run_message = Some(format!(
                "Confirm loading {} before discarding the current task state.",
                path.display()
            ));
            self.last_run_is_error = false;
            false
        } else {
            self.load_document_from_path(path);
            true
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
            self.last_run_message = Some(format!(
                "Load cancelled; the current task was preserved instead of loading {}.",
                path.display()
            ));
            self.last_run_is_error = false;
        }
    }

    /// Request closing the main window and confirm before discarding unsaved work.
    pub(crate) fn request_window_close(&mut self) -> bool {
        if self.problem_switch_requires_confirmation() {
            self.pending_window_close = true;
            self.last_run_message = Some(
                "Confirm closing the window before discarding the current task state.".to_string(),
            );
            self.last_run_is_error = false;
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
        if self.calculation_state.is_active() {
            self.cancel_calculation();
        }
        true
    }

    /// Cancel the queued close request without touching the task.
    pub(crate) fn cancel_pending_window_close(&mut self) {
        if self.pending_window_close {
            self.pending_window_close = false;
            self.last_run_message =
                Some("Close cancelled; the current task was preserved.".to_string());
            self.last_run_is_error = false;
        }
    }

    /// Requests a problem switch and either applies it immediately or queues
    /// a confirmation dialog if the current task is not safe to discard.
    pub(crate) fn request_problem_switch(&mut self, problem: ProblemsEnum) -> bool {
        if self.selected_problem == problem {
            return true;
        }

        if self.problem_switch_requires_confirmation() {
            self.pending_problem_switch = Some(problem.clone());
            self.last_run_message = Some(format!(
                "Confirm switching to {} before discarding the current task state.",
                problem
            ));
            self.last_run_is_error = false;
            false
        } else {
            self.replace_problem(problem);
            true
        }
    }

    /// Confirms the queued problem switch and replaces the active problem.
    pub(crate) fn confirm_pending_problem_switch(&mut self) -> bool {
        let Some(problem) = self.pending_problem_switch.take() else {
            return false;
        };
        self.replace_problem(problem);
        true
    }

    /// Cancels the queued problem switch without mutating the current problem.
    pub(crate) fn cancel_pending_problem_switch(&mut self) {
        if self.pending_problem_switch.take().is_some() {
            self.last_run_message =
                Some("Problem switch cancelled; the current task was preserved.".to_string());
            self.last_run_is_error = false;
        }
    }

    /// Replace the active problem as one complete document-lifecycle transition.
    ///
    /// Retaining only the old file path while replacing the document allowed a
    /// later Save action to overwrite the previous task with unrelated data.
    /// Reconstructing the complete state also clears stale plots, run messages,
    /// edit buffers, and AOT confirmations.
    pub(crate) fn replace_problem(&mut self, problem: ProblemsEnum) {
        if let Some(cancel) = &self.calculation_cancel {
            cancel.store(true, Ordering::Release);
        }
        *self = Self::new_with_problem(problem);
    }
}
