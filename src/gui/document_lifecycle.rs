use std::path::PathBuf;

/// Typed outcome for one GUI document file operation.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum GuiFileOperationResult {
    Saved {
        path: PathBuf,
    },
    Loaded {
        path: PathBuf,
    },
    Failed {
        kind: GuiFileOperationKind,
        path: Option<PathBuf>,
        message: String,
    },
}

/// File operation kind used by the lifecycle state and status messages.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum GuiFileOperationKind {
    Save,
    Load,
}

impl GuiFileOperationKind {
    fn label(self) -> &'static str {
        match self {
            Self::Save => "Save",
            Self::Load => "Load",
        }
    }
}

impl GuiFileOperationResult {
    pub fn is_error(&self) -> bool {
        matches!(self, Self::Failed { .. })
    }

    pub fn status_text(&self) -> String {
        match self {
            Self::Saved { path } => format!("Saved document to {}.", path.display()),
            Self::Loaded { path } => format!("Loaded document from {}.", path.display()),
            Self::Failed {
                kind,
                path,
                message,
            } => match path {
                Some(path) => format!(
                    "{} failed for {}: {}",
                    kind.label(),
                    path.display(),
                    message
                ),
                None => format!("{} failed: {}", kind.label(), message),
            },
        }
    }
}

/// Explicit lifecycle state for the active GUI document.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct DocumentLifecycleState {
    pub current_path: Option<PathBuf>,
    pub dirty: bool,
    pub last_file_operation: Option<GuiFileOperationResult>,
    clean_fingerprint: Option<String>,
}

impl DocumentLifecycleState {
    /// Marks the current document as modified.
    pub fn mark_dirty(&mut self) {
        self.dirty = true;
    }

    /// Marks the current document as clean after a successful load or save.
    pub fn mark_clean(&mut self) {
        self.dirty = false;
    }

    /// Stores the current document snapshot as the clean baseline.
    pub fn mark_clean_snapshot(&mut self, fingerprint: String) {
        self.clean_fingerprint = Some(fingerprint);
        self.mark_clean();
    }

    /// Updates the dirty flag by comparing the current fingerprint to the last
    /// clean snapshot.
    pub fn sync_dirty_state(&mut self, fingerprint: String) {
        self.dirty = self.clean_fingerprint.as_deref() != Some(fingerprint.as_str());
    }

    /// Records the result of a file operation and updates lifecycle fields.
    pub fn record(&mut self, result: GuiFileOperationResult) {
        match &result {
            GuiFileOperationResult::Saved { path } | GuiFileOperationResult::Loaded { path } => {
                self.current_path = Some(path.clone());
                self.mark_clean();
            }
            GuiFileOperationResult::Failed { .. } => {}
        }

        self.last_file_operation = Some(result);
    }

    /// Records a successful load or save using the provided clean fingerprint.
    pub fn record_clean(&mut self, result: GuiFileOperationResult, fingerprint: String) {
        self.record(result);
        self.clean_fingerprint = Some(fingerprint);
        self.mark_clean();
    }

    pub fn last_status_text(&self) -> Option<String> {
        self.last_file_operation
            .as_ref()
            .map(GuiFileOperationResult::status_text)
    }

    pub fn last_status_is_error(&self) -> bool {
        self.last_file_operation
            .as_ref()
            .is_some_and(GuiFileOperationResult::is_error)
    }
}
