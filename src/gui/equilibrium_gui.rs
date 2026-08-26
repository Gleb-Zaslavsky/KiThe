//! First editor surface for the production chemical-equilibrium workflow.
//!
//! This view prepares and validates requests, then owns a background worker
//! channel so repository lookup and nonlinear solving never block egui.

use crate::Thermodynamics::ChemEquilibrium::equilibrium_diagnostics::EquilibriumDiagnosticEvent;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_execution::{
    EquilibriumExecutionControl, EquilibriumProgressEvent, EquilibriumProgressStage,
};
use crate::Thermodynamics::ChemEquilibrium::prelude::ThermoRepository;
use crate::Thermodynamics::thermo_lib_api::ThermoData;
use crate::gui::equilibrium_gui_candidate::EquilibriumGuiCandidatePreview;
use crate::gui::equilibrium_gui_execution::{
    EquilibriumGuiExecution, EquilibriumGuiPublication, EquilibriumGuiRunState,
    EquilibriumGuiRunTicket, EquilibriumGuiWorkerMessage,
};
use crate::gui::equilibrium_gui_model::{
    ComponentDraft, EquilibriumGuiDocument, EquilibriumGuiDocumentError, EquilibriumInventoryDraft,
    EquilibriumLookupDraft, EquilibriumPhaseModeDraft, EquilibriumProblemDraft,
    EquilibriumSolverDraft, GuiInterpolationSpace, GuiKeqValidationMode, GuiPhaseLifecycleTrace,
    GuiPhaseModel, GuiPhysicalState, GuiPlotScale, GuiPlotTarget, GuiRangeLifecycleTrace,
    GuiResamplingDraft, GuiResultBasis, GuiSolverBackend, GuiSolverCascadeBudgetDraft,
    GuiTraceSeedPolicyDraft, PhTemperatureBoundsDraft, TemperatureDraft,
    ValidatedEquilibriumGuiConfig, ValidationIssue,
};
use crate::gui::equilibrium_gui_plot::{EquilibriumGuiKiThePlotWindow, EquilibriumGuiPlotData};
use crate::gui::equilibrium_gui_request::{
    EquilibriumGuiRequestError, EquilibriumGuiSolveOutcome, EquilibriumGuiSolveRequest,
    build_equilibrium_request, select_equilibrium_candidates,
};
use crate::gui::equilibrium_gui_result::EquilibriumGuiResultSnapshot;
use crate::gui::gui_plot::PlotWindow;
use eframe::egui;
use std::collections::BTreeSet;
use std::collections::hash_map::DefaultHasher;
use std::hash::{Hash, Hasher};
use std::sync::Arc;
use std::sync::mpsc::{self, Receiver};
use std::thread;

enum EquilibriumGuiWorkerEvent {
    Diagnostic {
        ticket: EquilibriumGuiRunTicket,
        event: EquilibriumDiagnosticEvent,
    },
    Completed {
        ticket: EquilibriumGuiRunTicket,
        outcome: EquilibriumGuiSolveOutcome,
    },
    Failed {
        ticket: EquilibriumGuiRunTicket,
        error: String,
    },
}

enum EquilibriumGuiCandidateWorkerEvent {
    Completed {
        fingerprint: u64,
        preview: EquilibriumGuiCandidatePreview,
    },
    Failed {
        fingerprint: u64,
        error: String,
    },
}

enum EquilibriumGuiLibraryWorkerEvent {
    Loaded { libraries: Vec<String> },
    Failed { error: String },
}

/// State owned by the chemical-equilibrium editor window.
pub struct EquilibriumApp {
    pub document: EquilibriumGuiDocument,
    pub open: bool,
    validation_issues: Vec<ValidationIssue>,
    prepared_request: Option<EquilibriumGuiSolveRequest>,
    execution: EquilibriumGuiExecution<EquilibriumGuiResultSnapshot>,
    candidate_preview: Option<EquilibriumGuiCandidatePreview>,
    candidate_preview_fingerprint: Option<u64>,
    candidate_target_phase: String,
    candidate_worker_receiver: Option<Receiver<EquilibriumGuiCandidateWorkerEvent>>,
    library_worker_receiver: Option<Receiver<EquilibriumGuiLibraryWorkerEvent>>,
    available_libraries: Option<Vec<String>>,
    library_picker: String,
    embedded_plot: Option<PlotWindow>,
    kithe_plot: EquilibriumGuiKiThePlotWindow,
    worker_receiver: Option<Receiver<EquilibriumGuiWorkerEvent>>,
    execution_control: Option<EquilibriumExecutionControl>,
    progress_receiver: Option<Receiver<EquilibriumProgressEvent>>,
    last_progress: Option<EquilibriumProgressEvent>,
    /// Bounded transient worker observations. Accepted results replace this
    /// generic live view with their phase-qualified immutable trace.
    live_diagnostic_events: Vec<String>,
    status: String,
    prepared_fingerprint: Option<u64>,
    /// Display-only series selection. It is intentionally outside the
    /// serialized document and solver fingerprint.
    hidden_plot_series: BTreeSet<String>,
}

impl Default for EquilibriumApp {
    fn default() -> Self {
        Self::new()
    }
}

impl EquilibriumApp {
    pub fn new() -> Self {
        Self {
            document: EquilibriumGuiDocument::new(),
            open: true,
            validation_issues: Vec::new(),
            prepared_request: None,
            execution: EquilibriumGuiExecution::default(),
            candidate_preview: None,
            candidate_preview_fingerprint: None,
            candidate_target_phase: "gas".into(),
            candidate_worker_receiver: None,
            library_worker_receiver: None,
            available_libraries: None,
            library_picker: String::new(),
            embedded_plot: None,
            kithe_plot: EquilibriumGuiKiThePlotWindow::default(),
            worker_receiver: None,
            execution_control: None,
            progress_receiver: None,
            last_progress: None,
            live_diagnostic_events: Vec::new(),
            status: "Ready for a canonical P,T request".into(),
            prepared_fingerprint: None,
            hidden_plot_series: BTreeSet::new(),
        }
    }

    /// Serializes only the versioned editable document. Workers, repository
    /// handles, plots, and accepted results are intentionally excluded.
    pub fn save_document_json(&self) -> Result<String, serde_json::Error> {
        self.document.to_json()
    }

    /// Replaces the editable document and starts its runtime state clean.
    ///
    /// A late message from a dropped worker is harmless: the receiver is
    /// removed and the execution gate no longer owns its ticket or result.
    pub fn load_document_json(&mut self, json: &str) -> Result<(), EquilibriumGuiDocumentError> {
        let document = EquilibriumGuiDocument::from_json(json)?;
        if let Some(control) = &self.execution_control {
            control.request_cancel();
        }
        self.worker_receiver = None;
        self.candidate_worker_receiver = None;
        self.library_worker_receiver = None;
        self.execution.reset();
        self.document = document;
        self.validation_issues.clear();
        self.prepared_request = None;
        self.prepared_fingerprint = None;
        self.candidate_preview = None;
        self.candidate_preview_fingerprint = None;
        self.embedded_plot = None;
        self.kithe_plot = EquilibriumGuiKiThePlotWindow::default();
        self.execution_control = None;
        self.progress_receiver = None;
        self.last_progress = None;
        self.live_diagnostic_events.clear();
        self.hidden_plot_series.clear();
        self.status = "Document loaded; runtime state reset to Idle".into();
        Ok(())
    }

    /// Prepares a request without touching the repository or starting a solve.
    pub fn prepare_request(&mut self) -> Result<(), EquilibriumGuiRequestError> {
        self.validation_issues.clear();
        let validated = match self.document.validate_for_run() {
            Ok(config) => config,
            Err(report) => {
                self.validation_issues = report.issues;
                self.prepared_request = None;
                self.prepared_fingerprint = None;
                self.status = "Cannot prepare request: fix validation errors".into();
                return Err(EquilibriumGuiRequestError::FeatureUnavailable(
                    "invalid GUI document",
                ));
            }
        };
        self.prepare_validated_request(validated)
    }

    fn prepare_validated_request(
        &mut self,
        validated: ValidatedEquilibriumGuiConfig,
    ) -> Result<(), EquilibriumGuiRequestError> {
        let request = build_equilibrium_request(validated, None)?;
        self.prepared_fingerprint = Some(self.document_fingerprint());
        self.prepared_request = Some(request);
        self.status = "Canonical request prepared; solver worker is not started".into();
        Ok(())
    }

    pub fn prepared_request(&self) -> Option<&EquilibriumGuiSolveRequest> {
        self.prepared_request.as_ref()
    }

    pub fn prepared_request_is_current(&self) -> bool {
        self.prepared_fingerprint == Some(self.document_fingerprint())
    }

    /// Performs a repository-backed element preview for the current document.
    /// The caller should invoke this from a worker because real catalogs can
    /// contain many thousands of records.
    pub fn preview_candidates(
        &mut self,
        repository: Arc<ThermoRepository>,
    ) -> Result<(), EquilibriumGuiRequestError> {
        self.validation_issues.clear();
        let validated = match self.document.validate() {
            Ok(config) => config,
            Err(report) => {
                self.validation_issues = report.issues;
                self.candidate_preview = None;
                self.candidate_preview_fingerprint = None;
                return Err(EquilibriumGuiRequestError::FeatureUnavailable(
                    "invalid GUI document",
                ));
            }
        };
        let report = select_equilibrium_candidates(&validated, repository)?;
        self.candidate_preview = Some(EquilibriumGuiCandidatePreview::from_report(&report));
        self.candidate_preview_fingerprint = Some(self.candidate_query_fingerprint());
        self.status =
            "Candidate preview accepted; assign candidates to phases before running".into();
        Ok(())
    }

    pub fn candidate_preview(&self) -> Option<&EquilibriumGuiCandidatePreview> {
        (self.candidate_preview_fingerprint == Some(self.candidate_query_fingerprint()))
            .then_some(self.candidate_preview.as_ref())
            .flatten()
    }

    /// Starts local-catalog candidate discovery without blocking egui.
    pub fn start_candidate_preview(&mut self) -> bool {
        if self.candidate_worker_receiver.is_some() {
            return false;
        }
        let validated = match self.document.validate() {
            Ok(config) => config,
            Err(report) => {
                self.validation_issues = report.issues;
                self.status = "Cannot preview candidates: fix validation errors".into();
                return false;
            }
        };
        if !matches!(
            &validated.inventory,
            crate::gui::equilibrium_gui_model::ValidatedInventory::ElementCandidates { .. }
        ) {
            self.status = "Candidate preview requires element inventory mode".into();
            return false;
        }
        let fingerprint = self.candidate_query_fingerprint();
        let (sender, receiver) = mpsc::channel();
        self.candidate_worker_receiver = Some(receiver);
        self.status = "Loading local catalog and selecting candidates".into();
        thread::spawn(move || {
            let event = match ThermoData::try_new_fresh() {
                Ok(data) => match select_equilibrium_candidates(&validated, data.repository) {
                    Ok(report) => EquilibriumGuiCandidateWorkerEvent::Completed {
                        fingerprint,
                        preview: EquilibriumGuiCandidatePreview::from_report(&report),
                    },
                    Err(error) => EquilibriumGuiCandidateWorkerEvent::Failed {
                        fingerprint,
                        error: error.to_string(),
                    },
                },
                Err(error) => EquilibriumGuiCandidateWorkerEvent::Failed {
                    fingerprint,
                    error: format!("thermochemistry repository loading failed: {error}"),
                },
            };
            let _ = sender.send(event);
        });
        true
    }

    fn poll_candidate_worker(&mut self) {
        let Some(receiver) = self.candidate_worker_receiver.take() else {
            return;
        };
        match receiver.try_recv() {
            Ok(EquilibriumGuiCandidateWorkerEvent::Completed {
                fingerprint,
                preview,
            }) if fingerprint == self.candidate_query_fingerprint() => {
                self.candidate_preview = Some(preview);
                self.candidate_preview_fingerprint = Some(fingerprint);
                self.status =
                    "Candidate preview accepted; assign candidates to phases before running".into();
            }
            Ok(EquilibriumGuiCandidateWorkerEvent::Failed { fingerprint, error })
                if fingerprint == self.candidate_query_fingerprint() =>
            {
                self.status = format!("Candidate preview failed: {error}");
            }
            Ok(_) => {
                self.status = "Discarded candidate preview from an edited document".into();
            }
            Err(mpsc::TryRecvError::Empty) => {
                self.candidate_worker_receiver = Some(receiver);
            }
            Err(mpsc::TryRecvError::Disconnected) => {
                self.status = "Candidate preview worker disconnected".into();
            }
        }
    }

    pub fn cancel_candidate_preview(&mut self) -> bool {
        if self.candidate_worker_receiver.take().is_some() {
            self.status = "Candidate preview cancelled; late output will be discarded".into();
            true
        } else {
            false
        }
    }

    /// Starts repository catalog discovery without blocking the egui thread.
    /// The catalog is display metadata only; solving still owns its own
    /// repository snapshot and remains transactional.
    pub fn start_library_catalog_load(&mut self) -> bool {
        if self.library_worker_receiver.is_some() || self.available_libraries.is_some() {
            return false;
        }
        let (sender, receiver) = mpsc::channel();
        self.library_worker_receiver = Some(receiver);
        self.status = "Loading thermochemical library choices".into();
        thread::spawn(move || {
            let event = match ThermoData::try_new_fresh() {
                Ok(data) => EquilibriumGuiLibraryWorkerEvent::Loaded {
                    libraries: normalize_library_choices(data.repository.thermo_libs.as_ref()),
                },
                Err(error) => EquilibriumGuiLibraryWorkerEvent::Failed {
                    error: format!("thermochemistry repository loading failed: {error}"),
                },
            };
            let _ = sender.send(event);
        });
        true
    }

    fn poll_library_worker(&mut self) {
        let Some(receiver) = self.library_worker_receiver.take() else {
            return;
        };
        match receiver.try_recv() {
            Ok(EquilibriumGuiLibraryWorkerEvent::Loaded { libraries }) => {
                if self.library_picker.is_empty() {
                    if let Some(first) = libraries.first() {
                        self.library_picker = first.clone();
                    }
                }
                self.available_libraries = Some(libraries);
                self.status = "Thermochemical library choices loaded".into();
            }
            Ok(EquilibriumGuiLibraryWorkerEvent::Failed { error }) => {
                self.status = format!("Library catalog failed: {error}");
            }
            Err(mpsc::TryRecvError::Empty) => {
                self.library_worker_receiver = Some(receiver);
            }
            Err(mpsc::TryRecvError::Disconnected) => {
                self.status = "Library catalog worker disconnected".into();
            }
        }
    }

    /// Assigns one selected exact record to the currently chosen phase row.
    /// The initial amount is visible and editable in the phase editor; it is
    /// never inferred from catalog metadata.
    pub fn assign_candidate_to_target_phase(&mut self, record_key: &str) -> bool {
        let source_library = self.candidate_preview().and_then(|preview| {
            preview
                .rows()
                .iter()
                .find(|row| row.record_key() == Some(record_key))
                .map(|row| row.library().to_string())
        });
        let EquilibriumInventoryDraft::ElementCandidates { assignments, .. } =
            &mut self.document.config.inventory
        else {
            return false;
        };
        let Some(phase) = assignments
            .iter_mut()
            .find(|phase| phase.id == self.candidate_target_phase)
        else {
            return false;
        };
        if phase
            .components
            .iter()
            .any(|component| component.substance == record_key)
        {
            return false;
        }
        phase.components.push(ComponentDraft {
            substance: record_key.to_string(),
            initial_moles: "1.0".into(),
            source_library,
        });
        self.invalidate_prepared_request();
        true
    }

    pub fn unassign_candidate(&mut self, record_key: &str) -> bool {
        let EquilibriumInventoryDraft::ElementCandidates { assignments, .. } =
            &mut self.document.config.inventory
        else {
            return false;
        };
        let mut removed = false;
        for phase in assignments {
            let before = phase.components.len();
            phase
                .components
                .retain(|component| component.substance != record_key);
            removed |= before != phase.components.len();
        }
        if removed {
            self.invalidate_prepared_request();
        }
        removed
    }

    /// Returns the editable initial amount for an assigned candidate.
    pub fn candidate_initial_moles(&self, record_key: &str) -> Option<&str> {
        let EquilibriumInventoryDraft::ElementCandidates { assignments, .. } =
            &self.document.config.inventory
        else {
            return None;
        };
        assignments
            .iter()
            .flat_map(|phase| phase.components.iter())
            .find(|component| component.substance == record_key)
            .map(|component| component.initial_moles.as_str())
    }

    /// Updates an assigned candidate's initial amount without parsing it in
    /// the view. Numeric validation remains centralized in the document
    /// validator, so partially typed values can be edited normally.
    pub fn set_candidate_initial_moles(
        &mut self,
        record_key: &str,
        initial_moles: impl Into<String>,
    ) -> bool {
        let initial_moles = initial_moles.into();
        let EquilibriumInventoryDraft::ElementCandidates { assignments, .. } =
            &mut self.document.config.inventory
        else {
            return false;
        };
        let Some(component) = assignments
            .iter_mut()
            .flat_map(|phase| phase.components.iter_mut())
            .find(|component| component.substance == record_key)
        else {
            return false;
        };
        if component.initial_moles == initial_moles {
            return false;
        }
        component.initial_moles = initial_moles;
        self.invalidate_prepared_request();
        true
    }

    fn candidate_query_fingerprint(&self) -> u64 {
        let mut hasher = DefaultHasher::new();
        match &self.document.config.inventory {
            EquilibriumInventoryDraft::ElementCandidates {
                elements,
                candidate_policy,
                ..
            } => {
                if let Ok(json) = serde_json::to_string(&(
                    elements,
                    candidate_policy,
                    &self.document.config.lookup,
                    &self.document.config.problem,
                )) {
                    json.hash(&mut hasher);
                }
            }
            EquilibriumInventoryDraft::ExplicitPhases { .. } => {
                "explicit-inventory".hash(&mut hasher);
            }
        }
        hasher.finish()
    }

    fn is_candidate_assigned(&self, record_key: Option<&str>) -> bool {
        let Some(record_key) = record_key else {
            return false;
        };
        match &self.document.config.inventory {
            EquilibriumInventoryDraft::ElementCandidates { assignments, .. } => assignments
                .iter()
                .flat_map(|phase| phase.components.iter())
                .any(|component| component.substance == record_key),
            EquilibriumInventoryDraft::ExplicitPhases { .. } => false,
        }
    }

    /// Builds an embedded plot from the accepted immutable snapshot only.
    /// Plot creation is display work and never invokes the equilibrium worker.
    pub fn open_embedded_plot(&mut self) -> Result<(), String> {
        let data = self.plot_data()?;
        self.embedded_plot = Some(data.to_embedded_plot_window());
        Ok(())
    }

    pub fn open_kithe_plot(&mut self) -> Result<(), String> {
        let data = self.plot_data()?;
        self.kithe_plot.open_from_data(&data)
    }

    /// Returns whether a result series is currently visible in plot windows.
    pub fn plot_series_visible(&self, label: &str) -> bool {
        !self.hidden_plot_series.contains(label)
    }

    /// Changes only display state; accepted solver results remain current.
    pub fn set_plot_series_visible(&mut self, label: impl Into<String>, visible: bool) {
        let label = label.into();
        if visible {
            self.hidden_plot_series.remove(&label);
        } else {
            self.hidden_plot_series.insert(label);
        }
    }

    fn plot_series_labels(&self) -> Vec<String> {
        let Some(snapshot) = self.result_snapshot() else {
            return Vec::new();
        };
        match self.document.config.postprocessing.result_basis {
            GuiResultBasis::ComponentMoles | GuiResultBasis::MoleFractions => {
                snapshot.component_labels().to_vec()
            }
            GuiResultBasis::PhaseTotals => snapshot.phase_labels().to_vec(),
        }
    }

    fn plot_data(&self) -> Result<EquilibriumGuiPlotData, String> {
        let Some(snapshot) = self.result_snapshot() else {
            return Err("no accepted equilibrium result to plot".into());
        };
        if !self.accepted_result_is_current() {
            return Err("accepted equilibrium result belongs to a previous document".into());
        }
        let data = EquilibriumGuiPlotData::from_snapshot(
            snapshot,
            self.document.config.postprocessing.result_basis,
        )?;
        let visible = data
            .column_labels()
            .filter(|label| self.plot_series_visible(label))
            .collect::<Vec<_>>();
        let data = data.retain_columns(visible)?;
        let data = match &self.document.config.postprocessing.resampling {
            GuiResamplingDraft::None => Ok(data),
            GuiResamplingDraft::Pchip {
                output_points,
                interpolation_space,
                clamp,
            } => {
                let output_points = output_points
                    .trim()
                    .parse::<usize>()
                    .map_err(|_| "PCHIP display point count must be an integer".to_string())?;
                let space = match interpolation_space {
                    GuiInterpolationSpace::Linear => RustedSciThe::numerical::optimization::inter_n_extrapolate::InterpolationSpace::Linear,
                    GuiInterpolationSpace::Log => RustedSciThe::numerical::optimization::inter_n_extrapolate::InterpolationSpace::Log,
                };
                data.resample_pchip(output_points, space, *clamp)
            }
        }?;
        match self.document.config.postprocessing.y_scale {
            GuiPlotScale::Linear => Ok(data),
            GuiPlotScale::Log10 => data.with_log10_y_scale(),
        }
    }

    /// Opens a lifecycle ticket for a previously prepared current request.
    /// The caller owns the actual worker/executor and must return its outcome
    /// through [`Self::publish_outcome`] or [`Self::publish_failure`].
    pub fn begin_prepared_run(&mut self) -> Option<EquilibriumGuiRunTicket> {
        if !self.prepared_request_is_current() {
            self.status = "Prepared request is stale; prepare it again".into();
            return None;
        }
        let ticket = self.execution.begin(self.document_fingerprint());
        self.execution.mark_solving(ticket);
        self.status = "Solver worker is running".into();
        Some(ticket)
    }

    /// Starts one worker for the prepared canonical request.
    ///
    /// The request is moved into the thread. The worker sends exactly one
    /// terminal event; the egui thread performs the stale-fingerprint check
    /// and immutable snapshot conversion before publishing anything.
    pub fn start_prepared_run(&mut self) -> bool {
        if self.worker_receiver.is_some() {
            self.status = "A solver worker is already running".into();
            return false;
        }
        let Some(ticket) = self.begin_prepared_run() else {
            return false;
        };
        let Some(request) = self.prepared_request.take() else {
            self.publish_failure(ticket, "prepared request disappeared before worker start");
            return false;
        };
        let (progress_sender, progress_receiver) = mpsc::channel();
        let control = EquilibriumExecutionControl::new().with_progress_sink(move |event| {
            let _ = progress_sender.send(event);
        });
        let request = request.with_execution_control(control.clone());
        self.execution_control = Some(control);
        self.progress_receiver = Some(progress_receiver);
        self.last_progress = None;
        self.live_diagnostic_events.clear();
        let (sender, receiver) = mpsc::channel();
        self.worker_receiver = Some(receiver);
        let diagnostic_sender = sender.clone();
        let request = request.with_diagnostic_sink(move |event| {
            let _ = diagnostic_sender.send(EquilibriumGuiWorkerEvent::Diagnostic { ticket, event });
        });
        thread::spawn(move || {
            let event = match request.solve() {
                Ok(outcome) => EquilibriumGuiWorkerEvent::Completed { ticket, outcome },
                Err(error) => EquilibriumGuiWorkerEvent::Failed {
                    ticket,
                    error: error.to_string(),
                },
            };
            let _ = sender.send(event);
        });
        true
    }

    /// Drains terminal worker messages on the egui thread.
    pub fn poll_worker(&mut self) {
        let Some(receiver) = self.worker_receiver.take() else {
            return;
        };
        let mut terminal = false;
        let mut disconnected = false;
        loop {
            match receiver.try_recv() {
                Ok(EquilibriumGuiWorkerEvent::Diagnostic { ticket, event })
                    if self.execution.active_ticket() == Some(ticket) =>
                {
                    self.push_live_diagnostic(event);
                }
                Ok(EquilibriumGuiWorkerEvent::Diagnostic { .. }) => {}
                Ok(EquilibriumGuiWorkerEvent::Completed { ticket, outcome }) => {
                    self.publish_outcome(ticket, outcome);
                    self.execution_control = None;
                    self.progress_receiver = None;
                    terminal = true;
                    break;
                }
                Ok(EquilibriumGuiWorkerEvent::Failed { ticket, error }) => {
                    self.publish_failure(ticket, error);
                    self.execution_control = None;
                    self.progress_receiver = None;
                    terminal = true;
                    break;
                }
                Err(mpsc::TryRecvError::Empty) => break,
                Err(mpsc::TryRecvError::Disconnected) => {
                    disconnected = true;
                    break;
                }
            }
        }
        if !terminal {
            if disconnected {
                if let Some(ticket) = self.execution.active_ticket() {
                    self.publish_failure(ticket, "equilibrium worker disconnected");
                }
                self.execution_control = None;
                self.progress_receiver = None;
            } else {
                self.worker_receiver = Some(receiver);
            }
        }
    }

    fn push_live_diagnostic(&mut self, event: EquilibriumDiagnosticEvent) {
        const MAX_LIVE_EVENTS: usize = 64;
        if self.live_diagnostic_events.len() == MAX_LIVE_EVENTS {
            self.live_diagnostic_events.remove(0);
        }
        self.live_diagnostic_events
            .push(format_live_diagnostic_event(&event));
    }

    /// Drains non-terminal progress without blocking the egui thread.
    pub fn poll_progress(&mut self) {
        let Some(receiver) = self.progress_receiver.as_ref() else {
            return;
        };
        while let Ok(event) = receiver.try_recv() {
            self.last_progress = Some(event);
        }
    }

    pub fn last_progress(&self) -> Option<EquilibriumProgressEvent> {
        self.last_progress
    }

    /// Publishes an accepted worker outcome through the stale-result gate.
    pub fn publish_outcome(
        &mut self,
        ticket: EquilibriumGuiRunTicket,
        outcome: EquilibriumGuiSolveOutcome,
    ) -> EquilibriumGuiPublication {
        if ticket.fingerprint() != self.document_fingerprint() {
            self.execution.cancel(ticket);
            self.status = "Discarded result from an edited document".into();
            return EquilibriumGuiPublication::Stale;
        }
        let snapshot = match EquilibriumGuiResultSnapshot::from_outcome(outcome) {
            Ok(snapshot) => snapshot,
            Err(error) => return self.publish_failure(ticket, error),
        };
        let publication = self
            .execution
            .publish(EquilibriumGuiWorkerMessage::Completed {
                ticket,
                result: snapshot,
            });
        if publication == EquilibriumGuiPublication::Accepted {
            self.status = "Equilibrium result accepted".into();
        }
        publication
    }

    /// Publishes a worker failure without destroying the previous accepted
    /// snapshot.
    pub fn publish_failure(
        &mut self,
        ticket: EquilibriumGuiRunTicket,
        error: impl Into<String>,
    ) -> EquilibriumGuiPublication {
        if ticket.fingerprint() != self.document_fingerprint() {
            self.execution.cancel(ticket);
            self.status = "Discarded failure from an edited document".into();
            return EquilibriumGuiPublication::Stale;
        }
        let publication = self.execution.publish(EquilibriumGuiWorkerMessage::Failed {
            ticket,
            error: error.into(),
        });
        if publication == EquilibriumGuiPublication::Accepted {
            self.status = "Equilibrium worker failed".into();
        }
        publication
    }

    pub fn result_snapshot(&self) -> Option<&EquilibriumGuiResultSnapshot> {
        self.execution.accepted_result()
    }

    pub fn run_state(&self) -> EquilibriumGuiRunState {
        self.execution.state()
    }

    pub fn accepted_result_is_current(&self) -> bool {
        self.execution.result_matches(self.document_fingerprint())
    }

    pub fn last_error(&self) -> Option<&str> {
        self.execution.last_error()
    }

    pub fn cancel_run(&mut self) -> bool {
        let Some(ticket) = self.execution.active_ticket() else {
            return false;
        };
        let cancelled = self.execution.cancel(ticket);
        if cancelled {
            if let Some(control) = &self.execution_control {
                control.request_cancel();
            }
            self.status = "Cancellation requested; worker is winding down".into();
        }
        cancelled
    }

    pub fn document_fingerprint(&self) -> u64 {
        // Postprocessing is deliberately excluded: basis, PCHIP, visibility,
        // and plot target consume an accepted snapshot and never alter the
        // equilibrium request itself.
        let mut request_document = self.document.clone();
        request_document.config.postprocessing = Default::default();
        let mut hasher = DefaultHasher::new();
        if let Ok(json) = request_document.to_json() {
            json.hash(&mut hasher);
        }
        hasher.finish()
    }

    /// Renders the editor window. It is deliberately a normal egui view so it
    /// can be exercised by `egui_kittest` without launching a native window.
    pub fn show(&mut self, ctx: &egui::Context, open: &mut bool) {
        let request_fingerprint_before_frame = self.document_fingerprint();
        self.poll_candidate_worker();
        self.poll_library_worker();
        self.poll_progress();
        self.poll_worker();
        egui::Window::new("Chemical equilibrium")
            .open(open)
            .resizable(true)
            .default_size([760.0, 720.0])
            .show(ctx, |ui| {
                ui.heading("Chemical equilibrium");
                ui.label(
                    if matches!(
                        self.document.config.problem,
                        EquilibriumProblemDraft::FixedPh { .. }
                    ) {
                        "Production P,H workflow"
                    } else {
                        "Production P,T workflow"
                    },
                );
                ui.separator();

                self.render_problem(ui);
                self.render_inventory(ui);
                self.render_phase_policy(ui);
                self.render_lookup(ui);
                self.render_solver(ui);
                self.render_diagnostics(ui);
                self.render_postprocessing(ui);
                self.render_results(ui);
                if matches!(
                    self.document.config.postprocessing.plot_target,
                    GuiPlotTarget::Embedded | GuiPlotTarget::Both
                ) {
                    if ui.button("Open embedded plot").clicked() {
                        if let Err(error) = self.open_embedded_plot() {
                            self.status = format!("Cannot open equilibrium plot: {error}");
                        }
                    }
                }
                if matches!(
                    self.document.config.postprocessing.plot_target,
                    GuiPlotTarget::KiThePlot | GuiPlotTarget::Both
                ) {
                    if ui.button("Open KiThePlot editor").clicked() {
                        if let Err(error) = self.open_kithe_plot() {
                            self.status = format!("Cannot open KiThePlot editor: {error}");
                        }
                    }
                }

                ui.separator();
                ui.horizontal(|ui| {
                    if ui.button("Validate document").clicked() {
                        self.validation_issues = match self.document.validate_for_run() {
                            Ok(_) => {
                                self.status = "Document is valid for a production request".into();
                                Vec::new()
                            }
                            Err(report) => {
                                self.status = "Document needs attention".into();
                                report.issues
                            }
                        };
                    }
                    if ui.button("Prepare canonical request").clicked() {
                        let _ = self.prepare_request();
                    }
                    if self.prepared_request_is_current()
                        && self.worker_receiver.is_none()
                        && ui.button("Run prepared request").clicked()
                    {
                        self.start_prepared_run();
                    }
                });
                ui.label(&self.status);
                ui.label(format!("Run state: {:?}", self.run_state()));
                if let Some(progress) = self.last_progress {
                    ui.label(format_progress(progress));
                }
                if self.execution.active_ticket().is_some()
                    && !self.live_diagnostic_events.is_empty()
                {
                    ui.collapsing("Live phase lifecycle", |ui| {
                        egui::ScrollArea::vertical()
                            .max_height(180.0)
                            .show(ui, |ui| {
                                for event in &self.live_diagnostic_events {
                                    ui.label(event);
                                }
                            });
                    });
                }
                if self.execution.active_ticket().is_some() && ui.button("Cancel run").clicked() {
                    self.cancel_run();
                }
                if let Some(snapshot) = self.result_snapshot() {
                    ui.label(format!(
                        "Accepted result{}: {} point(s), {} component(s)",
                        if self.accepted_result_is_current() {
                            ""
                        } else {
                            " (previous document)"
                        },
                        snapshot.points().len(),
                        snapshot.component_labels().len()
                    ));
                }
                if let Some(fingerprint) = self.prepared_fingerprint {
                    let current = self.document_fingerprint() == fingerprint;
                    ui.label(if current {
                        "Prepared request matches current document"
                    } else {
                        "Prepared request is stale after an editor change"
                    });
                }
                for issue in &self.validation_issues {
                    ui.colored_label(
                        egui::Color32::from_rgb(210, 80, 70),
                        format!("{}: {}", issue.field, issue.message),
                    );
                }
            });
        if let Some(plot) = self.embedded_plot.as_mut() {
            plot.show(ctx);
            if !plot.visible {
                self.embedded_plot = None;
            }
        }
        self.kithe_plot.show(ctx);
        // Text fields and display-independent editor controls can change the
        // document without going through one dedicated button callback. A
        // frame-level comparison keeps prepared requests transactional.
        if self.document_fingerprint() != request_fingerprint_before_frame {
            self.invalidate_prepared_request();
        }
    }

    fn render_problem(&mut self, ui: &mut egui::Ui) {
        ui.collapsing("Problem", |ui| {
            let is_pt = matches!(
                self.document.config.problem,
                EquilibriumProblemDraft::FixedPt { .. }
            );
            ui.horizontal(|ui| {
                if ui.selectable_label(is_pt, "P,T = const").clicked() && !is_pt {
                    self.document.config.problem = EquilibriumProblemDraft::default();
                    self.invalidate_prepared_request();
                }
                if ui.selectable_label(!is_pt, "P,H = const").clicked() && is_pt {
                    self.document.config.problem = EquilibriumProblemDraft::FixedPh {
                        pressure_pa: "101325".into(),
                        reference_pressure_pa: "101325".into(),
                        target_enthalpy_j: "0".into(),
                        temperature_bounds: PhTemperatureBoundsDraft::default(),
                    };
                    self.invalidate_prepared_request();
                }
            });
            match &mut self.document.config.problem {
                EquilibriumProblemDraft::FixedPt {
                    pressure_pa,
                    reference_pressure_pa,
                    temperature,
                } => {
                    labeled_text(ui, "Pressure [Pa]", pressure_pa);
                    labeled_text(ui, "Reference pressure [Pa]", reference_pressure_pa);
                    render_temperature(ui, temperature);
                }
                EquilibriumProblemDraft::FixedPh {
                    pressure_pa,
                    reference_pressure_pa,
                    target_enthalpy_j,
                    temperature_bounds,
                } => {
                    labeled_text(ui, "Pressure [Pa]", pressure_pa);
                    labeled_text(ui, "Reference pressure [Pa]", reference_pressure_pa);
                    labeled_text(ui, "Target total enthalpy [J]", target_enthalpy_j);
                    render_ph_temperature_bounds(ui, temperature_bounds);
                    ui.label(
                        "The P,H solve uses the common selected thermochemistry interval; "
                            .to_string(),
                    );
                }
            }
        });
    }

    fn render_inventory(&mut self, ui: &mut egui::Ui) {
        ui.collapsing("Components and phases", |ui| {
            let mut changed = false;
            let explicit = matches!(
                self.document.config.inventory,
                EquilibriumInventoryDraft::ExplicitPhases { .. }
            );
            ui.horizontal(|ui| {
                if ui.selectable_label(explicit, "Explicit species").clicked() && !explicit {
                    self.document.config.inventory = EquilibriumInventoryDraft::default();
                    self.invalidate_prepared_request();
                }
                if ui
                    .selectable_label(!explicit, "Search by elements")
                    .clicked()
                    && explicit
                {
                    self.document.config.inventory = EquilibriumInventoryDraft::ElementCandidates {
                        elements: vec!["C".into(), "H".into(), "O".into()],
                        candidate_policy: Default::default(),
                        assignments: Vec::new(),
                    };
                    self.invalidate_prepared_request();
                }
                if ui.button("Simple ideal-gas preset").clicked() {
                    self.document = EquilibriumGuiDocument::simple_ideal_gas(["H2O"]);
                    self.invalidate_prepared_request();
                    self.status = "Simple ideal-gas preset applied".into();
                }
            });
            match &mut self.document.config.inventory {
                EquilibriumInventoryDraft::ExplicitPhases { phases }
                | EquilibriumInventoryDraft::ElementCandidates {
                    assignments: phases,
                    ..
                } => {
                    let mut remove_phase = None;
                    let mut move_phase = None;
                    let phase_count = phases.len();
                    for (phase_index, phase) in phases.iter_mut().enumerate() {
                        // Phase IDs are validated as unique before a request is
                        // built, so they are a better widget scope than the
                        // current vector index. Reordering therefore keeps text
                        // edit state attached to the same phase.
                        let phase_key = editor_phase_key(phase, phase_index);
                        ui.push_id(("equilibrium-phase", phase_key.as_str()), |ui| {
                            ui.group(|ui| {
                                ui.horizontal(|ui| {
                                    ui.label(format!("Phase '{}'", phase.id.trim()));
                                    if phase_index > 0 && ui.button("Up").clicked() {
                                        move_phase = Some((phase_index, phase_index - 1));
                                    }
                                    if phase_index + 1 < phase_count && ui.button("Down").clicked()
                                    {
                                        move_phase = Some((phase_index, phase_index + 1));
                                    }
                                    if ui.button("Remove").clicked() {
                                        remove_phase = Some(phase_index);
                                    }
                                });
                                labeled_text(ui, "Phase id", &mut phase.id);
                                ui.horizontal(|ui| {
                                    ui.label("Physical state");
                                    for (label, value) in [
                                        ("Gas", GuiPhysicalState::Gas),
                                        ("Liquid", GuiPhysicalState::Liquid),
                                        ("Solid", GuiPhysicalState::Solid),
                                    ] {
                                        ui.selectable_value(
                                            &mut phase.physical_state,
                                            value,
                                            label,
                                        );
                                    }
                                });
                                ui.horizontal(|ui| {
                                    ui.label("Model");
                                    ui.selectable_value(
                                        &mut phase.model,
                                        GuiPhaseModel::IdealGas,
                                        "Ideal gas",
                                    );
                                    ui.selectable_value(
                                        &mut phase.model,
                                        GuiPhaseModel::IdealSolution,
                                        "Ideal solution",
                                    );
                                    ui.selectable_value(
                                        &mut phase.model,
                                        GuiPhaseModel::PureCondensed,
                                        "Pure condensed",
                                    );
                                });
                                let mut remove_component = None;
                                let mut move_component = None;
                                let component_count = phase.components.len();
                                for (component_index, component) in
                                    phase.components.iter_mut().enumerate()
                                {
                                    // Substance keys are unique inside a
                                    // validated phase. The index is only a
                                    // deterministic fallback while a row is
                                    // still empty or temporarily duplicated.
                                    let component_key = editor_component_key(
                                        &phase_key,
                                        component,
                                        component_index,
                                    );
                                    ui.push_id(
                                        ("equilibrium-component", component_key.as_str()),
                                        |ui| {
                                            ui.horizontal(|ui| {
                                                labeled_text(
                                                    ui,
                                                    "Substance",
                                                    &mut component.substance,
                                                );
                                                labeled_text(
                                                    ui,
                                                    "Moles",
                                                    &mut component.initial_moles,
                                                );
                                                if let Some(source_library) =
                                                    &component.source_library
                                                {
                                                    ui.label(format!("pinned: {source_library}"));
                                                }
                                                if component_index > 0 && ui.button("Up").clicked()
                                                {
                                                    move_component = Some((
                                                        component_index,
                                                        component_index - 1,
                                                    ));
                                                }
                                                if component_index + 1 < component_count
                                                    && ui.button("Down").clicked()
                                                {
                                                    move_component = Some((
                                                        component_index,
                                                        component_index + 1,
                                                    ));
                                                }
                                                if ui.button("Remove").clicked() {
                                                    remove_component = Some(component_index);
                                                }
                                            });
                                        },
                                    );
                                }
                                if ui.button("Add component").clicked() {
                                    phase.components.push(Default::default());
                                    changed = true;
                                }
                                if let Some(index) = remove_component {
                                    phase.components.remove(index);
                                    changed = true;
                                } else if let Some((from, to)) = move_component {
                                    phase.components.swap(from, to);
                                    changed = true;
                                }
                            });
                        });
                    }
                    if ui.button("Add phase").clicked() {
                        phases.push(Default::default());
                        changed = true;
                    }
                    if let Some(index) = remove_phase {
                        phases.remove(index);
                        changed = true;
                    } else if let Some((from, to)) = move_phase {
                        phases.swap(from, to);
                        changed = true;
                    }
                    if let EquilibriumInventoryDraft::ElementCandidates {
                        elements,
                        candidate_policy,
                        assignments,
                    } = &mut self.document.config.inventory
                    {
                        ui.label("Candidate preview is required before phase assignment");
                        let element_count = elements.len();
                        let mut remove_element = None;
                        for (index, element) in elements.iter_mut().enumerate() {
                            // Element symbols are unique after validation; the
                            // index keeps the editor deterministic for an
                            // incomplete/duplicated draft.
                            let element_key = editor_element_key(element, index);
                            ui.push_id(("equilibrium-element", element_key.as_str()), |ui| {
                                ui.horizontal(|ui| {
                                    labeled_text(ui, "Element", element);
                                    if element_count > 1 && ui.button("Remove").clicked() {
                                        remove_element = Some(index);
                                    }
                                });
                            });
                        }
                        if ui.button("Add element").clicked() {
                            elements.push(String::new());
                            changed = true;
                        }
                        if let Some(index) = remove_element {
                            elements.remove(index);
                            changed = true;
                        }
                        ui.horizontal(|ui| {
                            ui.label("Element matching");
                            ui.selectable_value(
                                &mut candidate_policy.element_mode,
                                crate::gui::equilibrium_gui_model::GuiElementSearchMode::SubsetOf,
                                "Subset of",
                            );
                            ui.selectable_value(
                                &mut candidate_policy.element_mode,
                                crate::gui::equilibrium_gui_model::GuiElementSearchMode::Exact,
                                "Exact set",
                            );
                        });
                        ui.horizontal(|ui| {
                            ui.label("Allowed physical states");
                            for (label, state) in [
                                ("Gas", GuiPhysicalState::Gas),
                                ("Liquid", GuiPhysicalState::Liquid),
                                ("Solid", GuiPhysicalState::Solid),
                                ("Condensed", GuiPhysicalState::Condensed),
                            ] {
                                let mut enabled = candidate_policy.physical_states.contains(&state);
                                if ui.checkbox(&mut enabled, label).changed() {
                                    if enabled {
                                        if !candidate_policy.physical_states.contains(&state) {
                                            candidate_policy.physical_states.push(state);
                                        }
                                    } else {
                                        candidate_policy
                                            .physical_states
                                            .retain(|item| *item != state);
                                    }
                                    changed = true;
                                }
                            }
                        });
                        labeled_text(
                            ui,
                            "Candidate temperature lower [K]",
                            &mut candidate_policy.temperature_lower_k,
                        );
                        labeled_text(
                            ui,
                            "Candidate temperature upper [K]",
                            &mut candidate_policy.temperature_upper_k,
                        );
                        labeled_text(ui, "Max candidates", &mut candidate_policy.max_candidates);
                        ui.label(format!(
                            "Confirmed phase assignments: {}",
                            assignments.len()
                        ));
                    }
                }
            }
            self.render_candidate_preview(ui);
            if matches!(
                self.document.config.inventory,
                EquilibriumInventoryDraft::ElementCandidates { .. }
            ) {
                ui.horizontal(|ui| {
                    if self.candidate_worker_receiver.is_some() {
                        ui.label("Candidate worker is running");
                        if ui.button("Cancel candidate preview").clicked() {
                            self.cancel_candidate_preview();
                        }
                    } else if self.candidate_preview().is_none()
                        && ui.button("Preview candidates").clicked()
                    {
                        self.start_candidate_preview();
                    }
                });
            }
            if changed {
                self.invalidate_prepared_request();
            }
        });
    }

    fn render_candidate_preview(&mut self, ui: &mut egui::Ui) {
        let Some(preview) = self.candidate_preview() else {
            return;
        };
        let requested_elements = preview.requested_elements().join(", ");
        let selected_count = preview.selected_count();
        let rejected_count = preview.rejected_count();
        let rows = preview.rows().to_vec();
        let phase_ids = match &self.document.config.inventory {
            EquilibriumInventoryDraft::ElementCandidates { assignments, .. } => assignments
                .iter()
                .map(|phase| phase.id.clone())
                .collect::<Vec<_>>(),
            EquilibriumInventoryDraft::ExplicitPhases { .. } => Vec::new(),
        };
        let mut target_phase = self.candidate_target_phase.clone();
        if !phase_ids.iter().any(|phase_id| phase_id == &target_phase) {
            if let Some(first_phase) = phase_ids.first() {
                target_phase = first_phase.clone();
            }
        }
        let mut pending_assignment = None;
        let mut pending_unassignment = None;
        let mut pending_amount = None;
        ui.separator();
        ui.label(format!(
            "Candidate preview for {}: {} selected, {} rejected",
            requested_elements, selected_count, rejected_count
        ));
        if phase_ids.is_empty() {
            ui.label("Add a phase assignment before assigning candidates");
        } else {
            egui::ComboBox::from_label("Target phase")
                .selected_text(&target_phase)
                .show_ui(ui, |ui| {
                    for phase_id in &phase_ids {
                        ui.selectable_value(&mut target_phase, phase_id.clone(), phase_id);
                    }
                });
        }
        egui::ScrollArea::vertical()
            .max_height(220.0)
            .show(ui, |ui| {
                egui::Grid::new("equilibrium-candidate-preview")
                    .striped(true)
                    .show(ui, |ui| {
                        for heading in [
                            "Use",
                            "Substance",
                            "Record key",
                            "Library",
                            "Physical state",
                            "Temperature",
                            "Assignment",
                            "Initial moles",
                            "Action",
                            "Decision",
                        ] {
                            ui.label(heading);
                        }
                        ui.end_row();
                        for row in &rows {
                            let assigned = self.is_candidate_assigned(row.record_key());
                            ui.label(if row.included() { "yes" } else { "no" });
                            ui.label(row.substance());
                            ui.label(row.record_key().unwrap_or("<not reported>"));
                            ui.label(row.library());
                            ui.label(format!("{:?}", row.physical_state()));
                            ui.label(format!("{:?}", row.temperature_support()));
                            ui.label(if assigned { "assigned" } else { "unassigned" });
                            if assigned {
                                let mut amount = self
                                    .candidate_initial_moles(row.record_key().unwrap_or_default())
                                    .unwrap_or_default()
                                    .to_string();
                                if ui.text_edit_singleline(&mut amount).changed() {
                                    pending_amount = row
                                        .record_key()
                                        .map(|record_key| (record_key.to_string(), amount));
                                }
                            } else {
                                ui.label("-");
                            }
                            if assigned && row.record_key().is_some() {
                                if ui.button("Unassign").clicked() {
                                    pending_unassignment = row.record_key().map(str::to_string);
                                }
                            } else if row.included()
                                && row.record_key().is_some()
                                && !phase_ids.is_empty()
                                && ui.button("Assign").clicked()
                            {
                                pending_assignment = row.record_key().map(str::to_string);
                            } else {
                                ui.label("");
                            }
                            ui.label(
                                row.rejection()
                                    .map(|reason| format!("rejected: {reason:?}"))
                                    .unwrap_or_else(|| "selected".into()),
                            );
                            ui.end_row();
                        }
                    });
            });
        if target_phase != self.candidate_target_phase {
            self.candidate_target_phase = target_phase;
        }
        if let Some(record_key) = pending_assignment {
            self.assign_candidate_to_target_phase(&record_key);
        }
        if let Some(record_key) = pending_unassignment {
            self.unassign_candidate(&record_key);
        }
        if let Some((record_key, amount)) = pending_amount {
            self.set_candidate_initial_moles(&record_key, amount);
        }
    }

    fn render_lookup(&mut self, ui: &mut egui::Ui) {
        ui.collapsing("Library lookup", |ui| {
            let mut changed = false;
            let available_libraries = self.available_libraries.clone();
            if let Some(libraries) = available_libraries.as_ref() {
                ui.label(format!(
                    "Available thermochemical libraries: {}",
                    libraries.len()
                ));
                ui.horizontal(|ui| {
                    egui::ComboBox::from_label("Catalog library")
                        .selected_text(if self.library_picker.is_empty() {
                            "Select library"
                        } else {
                            self.library_picker.as_str()
                        })
                        .show_ui(ui, |ui| {
                            for library in libraries {
                                ui.selectable_value(
                                    &mut self.library_picker,
                                    library.clone(),
                                    library,
                                );
                            }
                        });
                });
            } else if self.library_worker_receiver.is_some() {
                ui.label("Loading library choices...");
            } else if ui.button("Load local library choices").clicked() {
                self.start_library_catalog_load();
            }
            let default_policy =
                matches!(self.document.config.lookup, EquilibriumLookupDraft::Default);
            ui.horizontal(|ui| {
                if ui
                    .selectable_label(default_policy, "Engine default")
                    .clicked()
                    && !default_policy
                {
                    self.document.config.lookup = EquilibriumLookupDraft::Default;
                    self.invalidate_prepared_request();
                }
                if ui
                    .selectable_label(!default_policy, "Explicit policy")
                    .clicked()
                    && default_policy
                {
                    self.document.config.lookup = EquilibriumLookupDraft::Explicit {
                        priority_libraries: Vec::new(),
                        permitted_libraries: Vec::new(),
                        explicit_search_instructions: Default::default(),
                        search_in_nist: false,
                    };
                    self.invalidate_prepared_request();
                }
            });
            if let EquilibriumLookupDraft::Explicit {
                priority_libraries,
                permitted_libraries,
                search_in_nist,
                ..
            } = &mut self.document.config.lookup
            {
                ui.label("Priority libraries (ordered)");
                let mut remove_priority = None;
                for (index, library) in priority_libraries.iter_mut().enumerate() {
                    ui.horizontal(|ui| {
                        labeled_text(ui, &format!("Priority {}", index + 1), library);
                        if ui.button("Remove").clicked() {
                            remove_priority = Some(index);
                        }
                    });
                }
                if ui.button("Add priority library").clicked() {
                    priority_libraries.push(if self.library_picker.is_empty() {
                        String::new()
                    } else {
                        self.library_picker.clone()
                    });
                    changed = true;
                }
                if let Some(index) = remove_priority {
                    priority_libraries.remove(index);
                    changed = true;
                }

                ui.separator();
                ui.label("Permitted libraries (closed candidate set)");
                let mut remove_permitted = None;
                for (index, library) in permitted_libraries.iter_mut().enumerate() {
                    ui.horizontal(|ui| {
                        labeled_text(ui, &format!("Permitted {}", index + 1), library);
                        if ui.button("Remove").clicked() {
                            remove_permitted = Some(index);
                        }
                    });
                }
                if ui.button("Add permitted library").clicked() {
                    permitted_libraries.push(if self.library_picker.is_empty() {
                        String::new()
                    } else {
                        self.library_picker.clone()
                    });
                    changed = true;
                }
                if let Some(index) = remove_permitted {
                    permitted_libraries.remove(index);
                    changed = true;
                }
                if ui.checkbox(search_in_nist, "Allow NIST fallback").changed() {
                    changed = true;
                }
                ui.label("Library names are canonicalized and validated at request build time");
            }
            if changed {
                self.invalidate_prepared_request();
            }
        });
    }

    fn render_phase_policy(&mut self, ui: &mut egui::Ui) {
        ui.collapsing("Phase policy", |ui| {
            let fixed = matches!(
                self.document.config.phase_mode,
                EquilibriumPhaseModeDraft::FixedDeclared
            );
            ui.horizontal(|ui| {
                if ui
                    .selectable_label(fixed, "Fixed declared phases")
                    .clicked()
                    && !fixed
                {
                    self.document.config.phase_mode = EquilibriumPhaseModeDraft::FixedDeclared;
                    self.invalidate_prepared_request();
                }
                if ui
                    .selectable_label(!fixed, "Bounded phase control")
                    .clicked()
                    && fixed
                {
                    self.document.config.phase_mode = EquilibriumPhaseModeDraft::Bounded {
                        phase_epsilon: "1e-12".into(),
                        dg_create: "-1e-6".into(),
                        dg_keep: "1e-8".into(),
                        max_phase_iterations: "20".into(),
                    };
                    self.invalidate_prepared_request();
                }
            });
            if let EquilibriumPhaseModeDraft::Bounded {
                phase_epsilon,
                dg_create,
                dg_keep,
                max_phase_iterations,
            } = &mut self.document.config.phase_mode
            {
                labeled_text(ui, "Phase epsilon", phase_epsilon);
                labeled_text(ui, "Creation driving force", dg_create);
                labeled_text(ui, "Keep driving force", dg_keep);
                labeled_text(ui, "Maximum phase iterations", max_phase_iterations);
                ui.label(
                    "A phase is retained inside the hysteresis band; transitions are reported after solve",
                );
            }
        });
    }

    fn render_solver(&mut self, ui: &mut egui::Ui) {
        let mut invalidate = false;
        ui.collapsing("Solver", |ui| {
            let selection = &mut self.document.config.solver.selection;
            let production = matches!(selection, EquilibriumSolverDraft::ProductionDefault);
            if ui
                .selectable_label(production, "Production default cascade")
                .clicked()
                && !production
            {
                *selection = EquilibriumSolverDraft::ProductionDefault;
            }
            let selected_backend = match selection {
                EquilibriumSolverDraft::ProductionDefault => "Production default cascade".into(),
                EquilibriumSolverDraft::SingleBackend { backend } => backend.label().to_string(),
                EquilibriumSolverDraft::CustomCascade { backends } => format!(
                    "Custom cascade ({} backend{})",
                    backends.len(),
                    if backends.len() == 1 { "" } else { "s" }
                ),
            };
            egui::ComboBox::from_label("Concrete backend")
                .selected_text(selected_backend)
                .show_ui(ui, |ui| {
                    for backend in GuiSolverBackend::ALL {
                        ui.selectable_value(
                            selection,
                            EquilibriumSolverDraft::SingleBackend { backend },
                            backend.label(),
                        );
                    }
                });
            if ui.button("Use custom cascade").clicked()
                && !matches!(selection, EquilibriumSolverDraft::CustomCascade { .. })
            {
                *selection = EquilibriumSolverDraft::CustomCascade {
                    backends: vec![GuiSolverBackend::LegacyNr, GuiSolverBackend::LegacyLm],
                };
                invalidate = true;
            }
            if let EquilibriumSolverDraft::CustomCascade { backends } = selection {
                ui.label("Ordered fallback sequence");
                let mut remove_index = None;
                let mut move_up = None;
                let mut move_down = None;
                for index in 0..backends.len() {
                    let backend = backends[index];
                    ui.horizontal(|ui| {
                        ui.label(format!("{}.", index + 1));
                        egui::ComboBox::from_id_salt(("equilibrium-cascade", index))
                            .selected_text(backend.label())
                            .show_ui(ui, |ui| {
                                for candidate in GuiSolverBackend::ALL {
                                    ui.selectable_value(
                                        &mut backends[index],
                                        candidate,
                                        candidate.label(),
                                    );
                                }
                            });
                        if ui.small_button("Up").clicked() && index > 0 {
                            move_up = Some(index);
                        }
                        if ui.small_button("Down").clicked() && index + 1 < backends.len() {
                            move_down = Some(index);
                        }
                        if ui.small_button("Remove").clicked() {
                            remove_index = Some(index);
                        }
                    });
                }
                if ui.button("Add backend").clicked() {
                    backends.push(GuiSolverBackend::LegacyTr);
                    invalidate = true;
                }
                if let Some(index) = move_up {
                    backends.swap(index - 1, index);
                    invalidate = true;
                }
                if let Some(index) = move_down {
                    backends.swap(index, index + 1);
                    invalidate = true;
                }
                if let Some(index) = remove_index {
                    backends.remove(index);
                    invalidate = true;
                }
            }
            labeled_text(
                ui,
                "Tolerance override",
                &mut self.document.config.solver.overrides.tolerance,
            );
            labeled_text(
                ui,
                "Max iterations override",
                &mut self.document.config.solver.overrides.max_iterations,
            );
            let mut budget_override = self
                .document
                .config
                .solver
                .overrides
                .cascade_budget
                .is_some();
            if ui
                .checkbox(&mut budget_override, "Override cascade budget")
                .changed()
            {
                self.document.config.solver.overrides.cascade_budget = if budget_override {
                    Some(GuiSolverCascadeBudgetDraft {
                        max_attempts: "2".into(),
                        max_iterations_per_attempt: "100".into(),
                        max_total_iterations: "200".into(),
                    })
                } else {
                    None
                };
            }
            if let Some(budget) = self
                .document
                .config
                .solver
                .overrides
                .cascade_budget
                .as_mut()
            {
                labeled_text(ui, "Cascade max attempts", &mut budget.max_attempts);
                labeled_text(
                    ui,
                    "Cascade iterations per attempt",
                    &mut budget.max_iterations_per_attempt,
                );
                labeled_text(
                    ui,
                    "Cascade total iterations",
                    &mut budget.max_total_iterations,
                );
            }
            ui.checkbox(
                &mut self.document.config.solver.overrides.scaling_enabled,
                "Enable scaling",
            );
            let mut trace_override = self
                .document
                .config
                .solver
                .overrides
                .trace_seed_policy
                .is_some();
            if ui
                .checkbox(&mut trace_override, "Override trace-species seed policy")
                .changed()
            {
                self.document.config.solver.overrides.trace_seed_policy = if trace_override {
                    Some(GuiTraceSeedPolicyDraft::Absolute {
                        floor: "1e-30".into(),
                    })
                } else {
                    None
                };
                invalidate = true;
            }
            if self
                .document
                .config
                .solver
                .overrides
                .trace_seed_policy
                .is_some()
            {
                let absolute = matches!(
                    self.document.config.solver.overrides.trace_seed_policy,
                    Some(GuiTraceSeedPolicyDraft::Absolute { .. })
                );
                ui.horizontal(|ui| {
                    ui.label("Trace seed strategy");
                    if ui.selectable_label(absolute, "Absolute floor").clicked() && !absolute {
                        self.document.config.solver.overrides.trace_seed_policy =
                            Some(GuiTraceSeedPolicyDraft::Absolute {
                                floor: "1e-30".into(),
                            });
                        invalidate = true;
                    }
                    if ui
                        .selectable_label(!absolute, "Relative to largest mole")
                        .clicked()
                        && absolute
                    {
                        self.document.config.solver.overrides.trace_seed_policy =
                            Some(GuiTraceSeedPolicyDraft::RelativeToLargestInitialMole {
                                fraction: "1e-12".into(),
                                minimum_floor: "1e-30".into(),
                            });
                        invalidate = true;
                    }
                });
            }
            if let Some(policy) = self
                .document
                .config
                .solver
                .overrides
                .trace_seed_policy
                .as_mut()
            {
                match policy {
                    GuiTraceSeedPolicyDraft::Absolute { floor } => {
                        labeled_text(ui, "Trace floor", floor);
                    }
                    GuiTraceSeedPolicyDraft::RelativeToLargestInitialMole {
                        fraction,
                        minimum_floor,
                    } => {
                        labeled_text(ui, "Trace fraction", fraction);
                        labeled_text(ui, "Minimum trace floor", minimum_floor);
                    }
                }
            }
        });
        if invalidate {
            self.invalidate_prepared_request();
        }
    }

    fn render_diagnostics(&mut self, ui: &mut egui::Ui) {
        ui.collapsing("Diagnostics", |ui| {
            let diagnostics = &mut self.document.config.diagnostics;
            ui.checkbox(&mut diagnostics.collect_timing, "Collect timing");
            ui.checkbox(
                &mut diagnostics.retain_backend_attempts,
                "Retain backend attempts",
            );
            ui.checkbox(
                &mut diagnostics.retain_conservation_report,
                "Retain conservation report",
            );
            ui.checkbox(
                &mut diagnostics.retain_phase_transitions,
                "Retain phase transitions",
            );
            egui::ComboBox::from_label("Phase lifecycle trace")
                .selected_text(match diagnostics.phase_lifecycle_trace {
                    GuiPhaseLifecycleTrace::Off => "Off",
                    GuiPhaseLifecycleTrace::Summary => "Summary",
                    GuiPhaseLifecycleTrace::PhaseLifecycle => "Lifecycle",
                    GuiPhaseLifecycleTrace::Detailed => "Detailed",
                })
                .show_ui(ui, |ui| {
                    ui.selectable_value(
                        &mut diagnostics.phase_lifecycle_trace,
                        GuiPhaseLifecycleTrace::Off,
                        "Off",
                    );
                    ui.selectable_value(
                        &mut diagnostics.phase_lifecycle_trace,
                        GuiPhaseLifecycleTrace::Summary,
                        "Summary",
                    );
                    ui.selectable_value(
                        &mut diagnostics.phase_lifecycle_trace,
                        GuiPhaseLifecycleTrace::PhaseLifecycle,
                        "Lifecycle",
                    );
                    ui.selectable_value(
                        &mut diagnostics.phase_lifecycle_trace,
                        GuiPhaseLifecycleTrace::Detailed,
                        "Detailed",
                    );
                });
            egui::ComboBox::from_label("Range lifecycle trace")
                .selected_text(match diagnostics.range_lifecycle_trace {
                    GuiRangeLifecycleTrace::Endpoints => "Endpoints",
                    GuiRangeLifecycleTrace::TransitionsOnly => "Transitions only",
                    GuiRangeLifecycleTrace::EveryPoint => "Every point",
                })
                .show_ui(ui, |ui| {
                    ui.selectable_value(
                        &mut diagnostics.range_lifecycle_trace,
                        GuiRangeLifecycleTrace::Endpoints,
                        "Endpoints",
                    );
                    ui.selectable_value(
                        &mut diagnostics.range_lifecycle_trace,
                        GuiRangeLifecycleTrace::TransitionsOnly,
                        "Transitions only",
                    );
                    ui.selectable_value(
                        &mut diagnostics.range_lifecycle_trace,
                        GuiRangeLifecycleTrace::EveryPoint,
                        "Every point",
                    );
                });
            egui::ComboBox::from_label("Equilibrium-constant validation")
                .selected_text(match diagnostics.keq_validation {
                    GuiKeqValidationMode::Off => "Off",
                    GuiKeqValidationMode::WhenApplicable => "When applicable",
                    GuiKeqValidationMode::Required => "Required",
                })
                .show_ui(ui, |ui| {
                    ui.selectable_value(
                        &mut diagnostics.keq_validation,
                        GuiKeqValidationMode::Off,
                        "Off",
                    );
                    ui.selectable_value(
                        &mut diagnostics.keq_validation,
                        GuiKeqValidationMode::WhenApplicable,
                        "When applicable",
                    );
                    ui.selectable_value(
                        &mut diagnostics.keq_validation,
                        GuiKeqValidationMode::Required,
                        "Required",
                    );
                });
        });
    }

    fn render_postprocessing(&mut self, ui: &mut egui::Ui) {
        ui.collapsing("Postprocessing and plots", |ui| {
            ui.horizontal(|ui| {
                ui.label("Result basis");
                ui.selectable_value(
                    &mut self.document.config.postprocessing.result_basis,
                    GuiResultBasis::ComponentMoles,
                    "Moles",
                );
                ui.selectable_value(
                    &mut self.document.config.postprocessing.result_basis,
                    GuiResultBasis::MoleFractions,
                    "Mole fractions",
                );
                ui.selectable_value(
                    &mut self.document.config.postprocessing.result_basis,
                    GuiResultBasis::PhaseTotals,
                    "Phase totals",
                );
            });
            ui.horizontal(|ui| {
                ui.label("Plot target");
                ui.selectable_value(
                    &mut self.document.config.postprocessing.plot_target,
                    GuiPlotTarget::None,
                    "None",
                );
                ui.selectable_value(
                    &mut self.document.config.postprocessing.plot_target,
                    GuiPlotTarget::Embedded,
                    "Embedded",
                );
                ui.selectable_value(
                    &mut self.document.config.postprocessing.plot_target,
                    GuiPlotTarget::KiThePlot,
                    "KiThePlot",
                );
                ui.selectable_value(
                    &mut self.document.config.postprocessing.plot_target,
                    GuiPlotTarget::Both,
                    "Both",
                );
            });
            ui.horizontal(|ui| {
                ui.label("Y scale");
                ui.selectable_value(
                    &mut self.document.config.postprocessing.y_scale,
                    GuiPlotScale::Linear,
                    "Linear",
                );
                ui.selectable_value(
                    &mut self.document.config.postprocessing.y_scale,
                    GuiPlotScale::Log10,
                    "Log10",
                );
            });
            let plot_labels = self.plot_series_labels();
            if !plot_labels.is_empty() {
                ui.collapsing("Visible plot series", |ui| {
                    ui.label("Display filters do not change the accepted solver result");
                    for label in plot_labels {
                        let mut visible = self.plot_series_visible(&label);
                        if ui.checkbox(&mut visible, &label).changed() {
                            self.set_plot_series_visible(label, visible);
                        }
                    }
                });
            }
            let mut pchip = matches!(
                self.document.config.postprocessing.resampling,
                GuiResamplingDraft::Pchip { .. }
            );
            if ui
                .checkbox(&mut pchip, "PCHIP display resampling")
                .changed()
            {
                self.document.config.postprocessing.resampling = if pchip {
                    GuiResamplingDraft::Pchip {
                        output_points: "200".into(),
                        interpolation_space: GuiInterpolationSpace::Linear,
                        clamp: false,
                    }
                } else {
                    GuiResamplingDraft::None
                };
            }
            if let GuiResamplingDraft::Pchip {
                output_points,
                interpolation_space,
                clamp,
            } = &mut self.document.config.postprocessing.resampling
            {
                labeled_text(ui, "Display points", output_points);
                ui.horizontal(|ui| {
                    ui.label("Interpolation space");
                    ui.selectable_value(
                        interpolation_space,
                        GuiInterpolationSpace::Linear,
                        "Linear",
                    );
                    ui.selectable_value(interpolation_space, GuiInterpolationSpace::Log, "Log");
                });
                ui.checkbox(clamp, "Clamp display values to solved range");
            }
        });
    }

    /// Presents accepted solver data without exposing mutable engine objects
    /// or editable numeric result fields. The same snapshot is later consumed
    /// by both plot adapters.
    fn render_results(&self, ui: &mut egui::Ui) {
        let has_result = self.result_snapshot().is_some();
        egui::CollapsingHeader::new("Results")
            .default_open(has_result)
            .show(ui, |ui| {
                let Some(snapshot) = self.result_snapshot() else {
                    ui.label("No accepted equilibrium result yet");
                    return;
                };

                ui.label(format!(
                    "{} accepted point(s), {} phase-qualified component(s)",
                    snapshot.points().len(),
                    snapshot.component_labels().len()
                ));
                if let Some(enthalpy) = snapshot.enthalpy() {
                    ui.collapsing("P,H energy contract", |ui| {
                        egui::Grid::new("equilibrium-ph-result")
                            .striped(true)
                            .show(ui, |ui| {
                                ui.label("Target total enthalpy [J]");
                                ui.label(format!("{:.8e}", enthalpy.target_enthalpy_j()));
                                ui.end_row();
                                ui.label("Calculated total enthalpy [J]");
                                ui.label(format!("{:.8e}", enthalpy.calculated_enthalpy_j()));
                                ui.end_row();
                                ui.label("Absolute error [J]");
                                ui.label(format!("{:.8e}", enthalpy.enthalpy_error_j()));
                                ui.end_row();
                                ui.label("Relative error");
                                ui.label(format!("{:.8e}", enthalpy.relative_enthalpy_error()));
                                ui.end_row();
                            });
                    });
                }
                if let Some(diagnostics) = snapshot.ph_diagnostics() {
                    ui.collapsing("P,H solve diagnostics", |ui| {
                        egui::Grid::new("equilibrium-ph-diagnostics")
                            .striped(true)
                            .show(ui, |ui| {
                                ui.label("Solved temperature [K]");
                                ui.label(format!("{:.8}", diagnostics.solved_temperature_k()));
                                ui.end_row();
                                ui.label("Outer solve path");
                                ui.label(diagnostics.solve_path());
                                ui.end_row();
                                ui.label("Outer iterations");
                                ui.label(diagnostics.iterations().to_string());
                                ui.end_row();
                                if diagnostics.monolithic().is_none() {
                                    ui.label("Temperature trials");
                                    ui.label(diagnostics.trial_count().to_string());
                                    ui.end_row();
                                }
                                ui.label("Inner backend attempts");
                                ui.label(diagnostics.inner_backend_attempts().to_string());
                                ui.end_row();
                                ui.label("Inner nonlinear iterations");
                                ui.label(diagnostics.inner_nonlinear_iterations().to_string());
                                ui.end_row();
                                ui.label("Phase-control transitions");
                                ui.label(diagnostics.phase_control_transitions().to_string());
                                ui.end_row();
                                ui.label("Formulation builds");
                                ui.label(diagnostics.fixed_formulation_builds().to_string());
                                ui.end_row();
                                ui.label("Formulation reuses");
                                ui.label(diagnostics.fixed_formulation_reuses().to_string());
                                ui.end_row();
                                ui.label("Accepted error limit [J]");
                                ui.label(format!(
                                    "{:.8e}",
                                    diagnostics.accepted_enthalpy_error_limit_j()
                                ));
                                ui.end_row();
                                if diagnostics.timing_enabled() {
                                    ui.label("Outer timing [ms]");
                                    ui.label(format!("{:.3}", diagnostics.timing_total_ms()));
                                    ui.end_row();
                                }
                                if let Some(reason) = diagnostics.fallback_reason() {
                                    ui.label("Fallback reason");
                                    ui.label(reason);
                                    ui.end_row();
                                }
                            });
                        if let Some(monolithic) = diagnostics.monolithic() {
                            ui.collapsing("P,H monolithic evidence", |ui| {
                                ui.label(format!(
                                    "{} backend attempt(s)",
                                    monolithic.backend_attempts().len()
                                ));
                                egui::Grid::new("equilibrium-ph-monolithic-summary")
                                    .striped(true)
                                    .show(ui, |ui| {
                                        ui.label("Accepted backend");
                                        ui.label(monolithic.accepted_backend());
                                        ui.end_row();
                                        ui.label("Residual evaluations");
                                        ui.label(monolithic.residual_evaluations().to_string());
                                        ui.end_row();
                                        ui.label("Jacobian evaluations");
                                        ui.label(monolithic.jacobian_evaluations().to_string());
                                        ui.end_row();
                                        ui.label("Inner timing [ms]");
                                        ui.label(format!("{:.3}", monolithic.inner_timing_ms()));
                                        ui.end_row();
                                    });
                                ui.collapsing("Backend attempts", |ui| {
                                    for attempt in monolithic.backend_attempts() {
                                        ui.label(attempt);
                                    }
                                });
                                if !monolithic.acceptance_rows().is_empty() {
                                    ui.collapsing("Acceptance", |ui| {
                                        egui::Grid::new("equilibrium-ph-monolithic-acceptance")
                                            .striped(true)
                                            .show(ui, |ui| {
                                                for (label, value) in monolithic.acceptance_rows() {
                                                    ui.label(label);
                                                    ui.label(value);
                                                    ui.end_row();
                                                }
                                            });
                                    });
                                }
                                if !monolithic.phase_control_rows().is_empty() {
                                    ui.collapsing("Phase control", |ui| {
                                        egui::Grid::new("equilibrium-ph-monolithic-phases")
                                            .striped(true)
                                            .show(ui, |ui| {
                                                for (label, value) in
                                                    monolithic.phase_control_rows()
                                                {
                                                    ui.label(label);
                                                    ui.label(value);
                                                    ui.end_row();
                                                }
                                            });
                                    });
                                }
                            });
                        }
                        if !diagnostics.route_decisions().is_empty() {
                            ui.collapsing("P,H route decisions", |ui| {
                                egui::Grid::new("equilibrium-ph-route-decisions")
                                    .striped(true)
                                    .show(ui, |ui| {
                                        ui.label("From");
                                        ui.label("To");
                                        ui.label("Reason");
                                        ui.end_row();
                                        for decision in diagnostics.route_decisions() {
                                            ui.label(decision.from_route());
                                            ui.label(decision.to_route());
                                            ui.label(decision.reason());
                                            ui.end_row();
                                        }
                                    });
                            });
                        }
                        if diagnostics.monolithic().is_none() && !diagnostics.trials().is_empty() {
                            ui.collapsing("P,H temperature trials", |ui| {
                                ui.label(format!("{} trial(s)", diagnostics.trials().len()));
                                egui::ScrollArea::vertical()
                                    .max_height(240.0)
                                    .show(ui, |ui| {
                                        egui::Grid::new("equilibrium-ph-trials")
                                            .striped(true)
                                            .show(ui, |ui| {
                                                for heading in [
                                                    "Trial",
                                                    "Temperature [K]",
                                                    "Step",
                                                    "H(T) [J]",
                                                    "Error [J]",
                                                    "Scaled error",
                                                    "Inner attempts",
                                                    "Transitions",
                                                    "Time [ms]",
                                                ] {
                                                    ui.label(heading);
                                                }
                                                ui.end_row();
                                                for (index, trial) in
                                                    diagnostics.trials().iter().enumerate()
                                                {
                                                    ui.label((index + 1).to_string());
                                                    ui.label(format!(
                                                        "{:.6}",
                                                        trial.temperature_k()
                                                    ));
                                                    ui.label(trial.step_kind());
                                                    ui.label(format!(
                                                        "{:.6e}",
                                                        trial.total_enthalpy_j()
                                                    ));
                                                    ui.label(format!(
                                                        "{:.6e}",
                                                        trial.enthalpy_error_j()
                                                    ));
                                                    ui.label(format!(
                                                        "{:.6e}",
                                                        trial.scaled_error()
                                                    ));
                                                    ui.label(
                                                        trial.inner_backend_attempts().to_string(),
                                                    );
                                                    ui.label(
                                                        trial
                                                            .phase_control_transitions()
                                                            .to_string(),
                                                    );
                                                    ui.label(format!("{:.3}", trial.total_ms()));
                                                    ui.end_row();
                                                }
                                            });
                                    });
                            });
                        }
                    });
                }
                if let Some(range_report) = snapshot.range_report() {
                    ui.label(format!(
                        "Range: {:?}, builds={}, reuses={}, transitions={}",
                        range_report.direction(),
                        range_report.formulation_builds(),
                        range_report.formulation_reuses(),
                        range_report.phase_control_transitions()
                    ));
                }

                for (point_index, point) in snapshot.points().iter().enumerate() {
                    ui.collapsing(
                        format!("Point {}: {:.6} K", point_index + 1, point.temperature_k()),
                        |ui| {
                            egui::Grid::new(("equilibrium-result", point_index))
                                .striped(true)
                                .show(ui, |ui| {
                                    ui.label("Component");
                                    ui.label("Moles [mol]");
                                    ui.label("Mole fraction");
                                    ui.end_row();
                                    for (index, label) in
                                        snapshot.component_labels().iter().enumerate()
                                    {
                                        ui.label(label);
                                        ui.label(format!("{:.8e}", point.component_moles()[index]));
                                        ui.label(format!("{:.8e}", point.mole_fractions()[index]));
                                        ui.end_row();
                                    }
                                });
                            ui.label(format!(
                                "Phase status: {}",
                                point.phase_statuses().join(", ")
                            ));
                            ui.label(format!(
                                "Solver attempts: {}",
                                point.source().solve_report().attempt_count()
                            ));
                            ui.collapsing("Solve diagnostics", |ui| {
                                // Keep these as separate evidence sections. A
                                // single mixed table made it too easy to miss
                                // conservation or fallback failures while
                                // scanning a successful-looking result.
                                let validation = point.source().accepted_solution().validation();
                                ui.collapsing("Conservation", |ui| {
                                    egui::Grid::new(("equilibrium-conservation", point_index))
                                        .striped(true)
                                        .show(ui, |ui| {
                                            ui.label("Maximum elemental-balance error");
                                            ui.label(format!(
                                                "{:.6e}",
                                                validation.max_abs_element_balance_error
                                            ));
                                            ui.end_row();
                                            ui.label("Minimum reconstructed moles");
                                            ui.label(format!("{:.6e}", validation.min_moles));
                                            ui.end_row();
                                        });
                                });
                                ui.collapsing("Residuals", |ui| {
                                    egui::Grid::new(("equilibrium-residuals", point_index))
                                        .striped(true)
                                        .show(ui, |ui| {
                                            ui.label("Scaled residual L2 norm");
                                            ui.label(format!(
                                                "{:.6e}",
                                                validation.residual_l2_norm
                                            ));
                                            ui.end_row();
                                            ui.label("Raw residual L2 norm");
                                            ui.label(format!(
                                                "{:.6e}",
                                                validation.raw_residual_l2_norm
                                            ));
                                            ui.end_row();
                                            ui.label("Maximum residual");
                                            ui.label(format!(
                                                "{:.6e}",
                                                validation.max_abs_residual
                                            ));
                                            ui.end_row();
                                        });
                                });
                                ui.collapsing("Fallback attempts", |ui| {
                                    let report = point.source().solve_report();
                                    egui::Grid::new(("equilibrium-fallback", point_index))
                                        .striped(true)
                                        .show(ui, |ui| {
                                            ui.label("Backend attempts");
                                            ui.label(report.attempt_count().to_string());
                                            ui.end_row();
                                            ui.label("Started attempts");
                                            ui.label(report.started_attempt_count().to_string());
                                            ui.end_row();
                                            ui.label("Fallback attempts");
                                            ui.label(report.fallback_attempt_count().to_string());
                                            ui.end_row();
                                            ui.label("Accepted after fallback");
                                            ui.label(report.accepted_after_fallback().to_string());
                                            ui.end_row();
                                        });
                                    for index in 0..report.attempt_count() {
                                        if let Some(attempt) = report.attempt(index) {
                                            ui.label(attempt.summary());
                                        }
                                    }
                                });
                                ui.collapsing("K_eq validation", |ui| {
                                    let status = point
                                        .source()
                                        .keq_validation_status()
                                        .map(|status| format!("{status:?}"))
                                        .unwrap_or_else(|| "not requested".to_string());
                                    ui.label(status);
                                });
                                ui.collapsing("Phase and acceptance evidence", |ui| {
                                    egui::Grid::new(("equilibrium-phase-evidence", point_index))
                                        .striped(true)
                                        .show(ui, |ui| {
                                            for row in
                                                point.source().summary_rows().into_iter().filter(
                                                    |row| {
                                                        matches!(
                                                            row.section,
                                                            "backend"
                                                                | "phase"
                                                                | "phase_control"
                                                                | "acceptance"
                                                                | "complementarity"
                                                                | "phase_stability"
                                                        )
                                                    },
                                                )
                                            {
                                                ui.label(row.section);
                                                ui.label(row.label);
                                                ui.label(row.value);
                                                ui.end_row();
                                            }
                                        });
                                });
                            });
                            if let Some(trace) = point.lifecycle_trace() {
                                ui.collapsing("Phase lifecycle", |ui| {
                                    ui.label(format!(
                                        "{} trace: {} retained event(s)",
                                        trace.mode(),
                                        trace.events().len()
                                    ));
                                    egui::ScrollArea::vertical()
                                        .max_height(260.0)
                                        .show(ui, |ui| {
                                            egui::Grid::new((
                                                "equilibrium-lifecycle-trace",
                                                point_index,
                                            ))
                                            .striped(true)
                                            .show(ui, |ui| {
                                                ui.label("Decision");
                                                ui.label("Evidence");
                                                ui.end_row();
                                                for event in trace.events() {
                                                    ui.label(event.title());
                                                    ui.label(event.detail());
                                                    ui.end_row();
                                                }
                                            });
                                        });
                                    if trace.dropped_events() > 0 {
                                        ui.label(format!(
                                            "{} additional event(s) were omitted by the retention limit",
                                            trace.dropped_events()
                                        ));
                                    }
                                });
                            }
                            ui.collapsing("Lookup provenance", |ui| {
                                egui::Grid::new(("equilibrium-provenance", point_index))
                                    .striped(true)
                                    .show(ui, |ui| {
                                        ui.label("Component");
                                        ui.label("Library");
                                        ui.label("Record");
                                        ui.label("State");
                                        ui.end_row();
                                        for component in point.source().build_report().components()
                                        {
                                            let source = component.thermo_source();
                                            ui.label(component.component().label());
                                            ui.label(source.library());
                                            ui.label(source.record_key());
                                            ui.label(source.state());
                                            ui.end_row();
                                        }
                                    });
                            });
                        },
                    );
                }
            });
    }

    fn invalidate_prepared_request(&mut self) {
        if let Some(ticket) = self.execution.active_ticket() {
            if let Some(control) = &self.execution_control {
                control.request_cancel();
            }
            self.execution.cancel(ticket);
            self.status = "Editor changed; active equilibrium run is being cancelled".into();
        }
        self.prepared_request = None;
        self.prepared_fingerprint = None;
    }
}

fn format_progress(event: EquilibriumProgressEvent) -> String {
    match event.stage() {
        EquilibriumProgressStage::RepositoryLookup => "Progress: repository lookup complete".into(),
        EquilibriumProgressStage::FormulationPreparation => {
            if let Some(count) = event.point_count() {
                format!("Progress: formulation ready for {count} point(s)")
            } else {
                "Progress: formulation preparation complete".into()
            }
        }
        EquilibriumProgressStage::PointStarted => format_point_progress("started", event),
        EquilibriumProgressStage::PointAccepted => format_point_progress("accepted", event),
        EquilibriumProgressStage::TemperatureTrialStarted => {
            format_point_progress("temperature trial started", event)
        }
        EquilibriumProgressStage::TemperatureTrialAccepted => {
            format_point_progress("temperature trial accepted", event)
        }
        EquilibriumProgressStage::TemperatureTrialRejected => {
            format_point_progress("temperature trial rejected", event)
        }
        EquilibriumProgressStage::InnerSolveStarted => {
            format_point_progress("inner P,T solve started", event)
        }
        EquilibriumProgressStage::InnerSolveCompleted => {
            format_point_progress("inner P,T solve completed", event)
        }
        EquilibriumProgressStage::PhaseTransitionAccepted => {
            format_point_progress("phase transition accepted", event)
        }
        EquilibriumProgressStage::PublicationStarted => {
            "Progress: P,H publication validation started".into()
        }
        EquilibriumProgressStage::PublicationCompleted => "Progress: P,H result published".into(),
        EquilibriumProgressStage::InnerBackendAttemptStarted => {
            "Progress: nonlinear backend attempt started".into()
        }
        EquilibriumProgressStage::InnerBackendAttemptFinished => {
            "Progress: nonlinear backend attempt finished".into()
        }
    }
}

/// Produces a short transient worker message without duplicating the
/// phase-qualified renderer used by accepted snapshots. During a solve the
/// worker may not yet have an accepted layout; after publication the GUI uses
/// `EquilibriumGuiLifecycleTraceSnapshot` instead.
fn format_live_diagnostic_event(event: &EquilibriumDiagnosticEvent) -> String {
    match event {
        EquilibriumDiagnosticEvent::SolveStarted { .. } => "Solve started".into(),
        EquilibriumDiagnosticEvent::OuterIterationStarted { iteration, .. } => {
            format!("Outer iteration {iteration} started")
        }
        EquilibriumDiagnosticEvent::ActiveSetCandidateAccepted { iteration, .. } => {
            format!("Fixed active-set candidate accepted at iteration {iteration}")
        }
        EquilibriumDiagnosticEvent::ActiveSetCandidateRejected {
            iteration, message, ..
        } => format!("Candidate rejected at iteration {iteration}: {message}"),
        EquilibriumDiagnosticEvent::StabilityEvaluated { iteration, .. } => {
            format!("TPD stability evaluated at iteration {iteration}")
        }
        EquilibriumDiagnosticEvent::TransitionAccepted {
            iteration, reason, ..
        } => {
            format!("Phase transition accepted at iteration {iteration}: {reason:?}")
        }
        EquilibriumDiagnosticEvent::TransitionHeldByHysteresis {
            iteration,
            phase_index,
        } => format!("Hysteresis retained phase {phase_index} at iteration {iteration}"),
        EquilibriumDiagnosticEvent::RecoveryProbeStarted {
            iteration,
            removed_phase_index,
            ..
        } => format!(
            "Boundary recovery started at iteration {iteration}: remove phase {removed_phase_index}"
        ),
        EquilibriumDiagnosticEvent::RecoveryProbeAccepted {
            iteration,
            phase_index,
            ..
        } => format!(
            "Boundary recovery accepted at iteration {iteration}: removed phase {phase_index}"
        ),
        EquilibriumDiagnosticEvent::RecoveryProbeRejected {
            iteration,
            removed_phase_index,
            message,
        } => format!(
            "Boundary recovery rejected at iteration {iteration}: remove phase {removed_phase_index}; {message}"
        ),
        EquilibriumDiagnosticEvent::ContinuationRestored { retained_seed, .. } => {
            format!("Continuation restored; seed_retained={retained_seed}")
        }
        EquilibriumDiagnosticEvent::PhaseControlBudgetExhausted {
            max_outer_iterations,
        } => format!("Phase-control budget exhausted after {max_outer_iterations} iteration(s)"),
        EquilibriumDiagnosticEvent::PhaseControlCycleDetected { iteration, .. } => {
            format!("Phase-control cycle detected at iteration {iteration}")
        }
        EquilibriumDiagnosticEvent::PhRouteFallback {
            from_route,
            to_route,
            error_kind,
            ..
        } => format!("P,H route fallback: {from_route:?} -> {to_route:?} ({error_kind:?})"),
        EquilibriumDiagnosticEvent::PhRouteFailed {
            route, error_kind, ..
        } => format!("P,H route failed: {route:?} ({error_kind:?})"),
        EquilibriumDiagnosticEvent::SolveFailed {
            continuation_restored,
            message,
        } => format!("Solve failed; continuation_restored={continuation_restored}; {message}"),
        EquilibriumDiagnosticEvent::SolveAccepted {
            outer_iterations,
            transition_count,
            ..
        } => format!(
            "Solve accepted: outer_iterations={outer_iterations}, transitions={transition_count}"
        ),
    }
}

fn format_point_progress(verb: &str, event: EquilibriumProgressEvent) -> String {
    match (
        event.point_index(),
        event.point_count(),
        event.temperature(),
    ) {
        (Some(index), Some(count), Some(temperature)) => format!(
            "Progress: point {}/{} {verb} at {temperature:.6} K",
            index + 1,
            count
        ),
        _ => format!("Progress: point {verb}"),
    }
}

fn labeled_text(ui: &mut egui::Ui, label: &str, value: &mut String) {
    ui.horizontal(|ui| {
        ui.label(label);
        ui.text_edit_singleline(value);
    });
}

/// Returns a stable egui scope for a phase editor row.
///
/// The semantic phase ID survives reorder operations. An empty ID is still
/// common while the user is typing, so the index is used only as a temporary
/// deterministic fallback until validation can assign the semantic identity.
fn editor_phase_key(phase: &crate::gui::equilibrium_gui_model::PhaseDraft, index: usize) -> String {
    let id = phase.id.trim();
    if id.is_empty() {
        format!("anonymous-{index}")
    } else {
        format!("id-{id}")
    }
}

/// Returns a stable egui scope for a component row within one phase.
fn editor_component_key(phase_key: &str, component: &ComponentDraft, index: usize) -> String {
    let substance = component.substance.trim();
    if substance.is_empty() {
        format!("{phase_key}-anonymous-{index}")
    } else {
        format!("{phase_key}-substance-{substance}")
    }
}

/// Returns a stable egui scope for an editable element row.
fn editor_element_key(element: &str, index: usize) -> String {
    let element = element.trim();
    if element.is_empty() {
        format!("anonymous-{index}")
    } else {
        format!("symbol-{element}")
    }
}

fn normalize_library_choices(libraries: &[String]) -> Vec<String> {
    let mut choices = libraries
        .iter()
        .map(|library| library.trim().to_string())
        .filter(|library| !library.is_empty())
        .collect::<Vec<_>>();
    choices.sort();
    choices.dedup();
    choices
}

fn render_temperature(ui: &mut egui::Ui, temperature: &mut TemperatureDraft) {
    let is_point = matches!(temperature, TemperatureDraft::Point { .. });
    ui.horizontal(|ui| {
        if ui.selectable_label(is_point, "Point").clicked() && !is_point {
            *temperature = TemperatureDraft::Point {
                temperature_k: "1000".into(),
            };
        }
        if ui.selectable_label(!is_point, "Range").clicked() && is_point {
            *temperature = TemperatureDraft::Range {
                start_k: "300".into(),
                end_k: "1500".into(),
                point_count: "20".into(),
            };
        }
    });
    match temperature {
        TemperatureDraft::Point { temperature_k } => {
            labeled_text(ui, "Temperature [K]", temperature_k);
        }
        TemperatureDraft::Range {
            start_k,
            end_k,
            point_count,
        } => {
            labeled_text(ui, "Start [K]", start_k);
            labeled_text(ui, "End [K]", end_k);
            labeled_text(ui, "Solved points", point_count);
        }
    }
}

fn render_ph_temperature_bounds(ui: &mut egui::Ui, bounds: &mut PhTemperatureBoundsDraft) {
    ui.label("P,H temperature solve");
    labeled_text(ui, "Lower bound [K]", &mut bounds.lower_k);
    labeled_text(ui, "Upper bound [K]", &mut bounds.upper_k);
    labeled_text(ui, "Initial seed [K]", &mut bounds.seed_k);
}

#[cfg(test)]
mod tests {
    use super::normalize_library_choices;

    #[test]
    fn library_choices_are_trimmed_sorted_and_deduplicated() {
        let choices = normalize_library_choices(&[
            " NASA_gas ".into(),
            "NIST".into(),
            "NASA_gas".into(),
            "".into(),
        ]);
        assert_eq!(choices, vec!["NASA_gas", "NIST"]);
    }
}
