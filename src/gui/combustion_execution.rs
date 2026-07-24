use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Arc, mpsc};
use std::time::Duration;

use RustedSciThe::command_interpreter::task_parser::DocumentMap;
use eframe::egui;
use log::warn;

use super::{
    BvpWorkerMessage, BvpWorkerOutcome, CalculationState, CombustionApp, PlotWindow,
    aot_toolchain_request, build_bvp_task_preview_snapshot, execute_bvp_calculation,
    print_bvp_task_preview, validate_aot_toolchain_fields, validate_bvp_gui_document,
};

impl CombustionApp {
    /// Executes the appropriate calculation based on the selected problem type.
    ///
    /// This method dispatches to the correct solver based on the current problem
    /// type selection. It passes the current document state to the solver,
    /// allowing users to run calculations with their configured parameters.
    ///
    /// # Problem Type Dispatch
    ///
    /// * `BVPSimple` - Creates a SimpleReactorTask and calls solve_from_parsed()
    /// * `None` - Shows a message that no calculation is available
    ///
    /// # Error Handling
    ///
    /// Solver errors are handled internally by the respective solver implementations.
    /// This method focuses on dispatch logic rather than error propagation.
    #[allow(dead_code)]
    pub(crate) fn run_calculation(&mut self) {
        let document = match self.validated_document_for_run() {
            Ok(document) => document,
            Err(message) => {
                self.calculation_state = CalculationState::Failed;
                self.last_run_message = Some(message);
                self.last_run_is_error = true;
                return;
            }
        };
        self.plot_window = None;
        self.calculation_state = CalculationState::Running { run_id: 0 };
        self.last_run_message = Some("Running calculation...".to_string());
        self.last_run_is_error = false;
        let outcome = execute_bvp_calculation(document, self.selected_problem.clone());
        self.apply_worker_outcome(outcome);
    }

    /// Prints a read-only BVP preview without mutating the loaded task data.
    ///
    /// The preview uses the same normalized snapshot that the solver handoff
    /// would see, but it stops short of starting any calculation.
    pub(crate) fn preview_bvp_task(&mut self) {
        let snapshot = build_bvp_task_preview_snapshot(&self.document, &self.selected_problem);
        print_bvp_task_preview(&snapshot);
        self.last_run_message = Some("Task preview printed to console.".to_string());
        self.last_run_is_error = false;
    }

    fn apply_worker_outcome(&mut self, outcome: BvpWorkerOutcome) {
        match outcome {
            BvpWorkerOutcome::Completed { plot } => {
                self.plot_window = plot
                    .map(|plot| PlotWindow::new(plot.arg, plot.values, plot.x_mesh, plot.solution));
                self.calculation_state = CalculationState::Completed;
                self.last_run_message = Some("Calculation completed successfully.".to_string());
                self.last_run_is_error = false;
            }
            BvpWorkerOutcome::Failed(message) => {
                warn!("BVP GUI calculation failed: {}", message);
                self.plot_window = None;
                self.calculation_state = CalculationState::Failed;
                self.last_run_message = Some(message);
                self.last_run_is_error = true;
            }
            BvpWorkerOutcome::Cancelled => {
                self.plot_window = None;
                self.calculation_state = CalculationState::Idle;
                self.last_run_message = Some("Calculation cancelled.".to_string());
                self.last_run_is_error = false;
            }
        }
    }

    /// Start one owned calculation worker while keeping egui responsive.
    fn start_calculation_worker(&mut self) {
        let document = match self.validated_document_for_run() {
            Ok(document) => document,
            Err(message) => {
                self.calculation_state = CalculationState::Failed;
                self.last_run_message = Some(message);
                self.last_run_is_error = true;
                return;
            }
        };
        let problem = self.selected_problem.clone();
        self.start_worker_job(move |worker_cancel| {
            if worker_cancel.load(Ordering::Acquire) {
                return BvpWorkerOutcome::Cancelled;
            }
            let calculated = execute_bvp_calculation(document, problem);
            if worker_cancel.load(Ordering::Acquire) {
                BvpWorkerOutcome::Cancelled
            } else {
                calculated
            }
        });
    }

    pub(crate) fn current_validation_fingerprint(&self) -> String {
        format!("{}\n{}", self.selected_problem, self.document_to_string())
    }

    /// Refresh cached validation only after the problem or document changes.
    pub(crate) fn refresh_validation_report(&mut self) {
        self.refresh_bvp_gui_snapshot();
        let fingerprint = self.current_validation_fingerprint();
        if self.validation_fingerprint.as_deref() == Some(fingerprint.as_str()) {
            return;
        }
        let (_, report) = validate_bvp_gui_document(&self.document, &self.selected_problem);
        self.validation_fingerprint = Some(fingerprint);
        self.validation_report = report;
    }

    /// Force preflight immediately before a synchronous or worker solve request.
    pub(crate) fn validated_document_for_run(&mut self) -> Result<DocumentMap, String> {
        self.refresh_bvp_gui_snapshot();
        let fingerprint = self.current_validation_fingerprint();
        let (document, report) = validate_bvp_gui_document(&self.document, &self.selected_problem);
        self.validation_fingerprint = Some(fingerprint);
        self.validation_report = report;
        if self.validation_report.is_valid() {
            Ok(document)
        } else {
            Err(self.validation_report.failure_message())
        }
    }

    /// Rebuild the typed BVP GUI snapshot without mutating the raw document.
    pub(crate) fn refresh_bvp_gui_snapshot(&mut self) {
        let (config, report) = super::BvpGuiConfig::from_document(&self.document);
        if self.bvp_gui_config != config {
            self.bvp_gui_config = config;
        }
        self.bvp_gui_migration_report = report;
    }

    fn start_worker_job<F>(&mut self, job: F)
    where
        F: FnOnce(Arc<AtomicBool>) -> BvpWorkerOutcome + Send + 'static,
    {
        if self.calculation_state.is_active() {
            self.last_run_message = Some("A calculation is already running.".to_string());
            self.last_run_is_error = true;
            return;
        }

        self.plot_window = None;
        self.last_run_message = Some("Running calculation...".to_string());
        self.last_run_is_error = false;
        let run_id = self.next_run_id;
        self.next_run_id = self.next_run_id.wrapping_add(1);
        let cancel = Arc::new(AtomicBool::new(false));
        let worker_cancel = Arc::clone(&cancel);
        let (sender, receiver) = mpsc::channel();

        match std::thread::Builder::new()
            .name(format!("kithe-bvp-{}", run_id))
            .spawn(move || {
                let outcome = job(worker_cancel);
                let _ = sender.send(BvpWorkerMessage { run_id, outcome });
            }) {
            Ok(_) => {
                self.calculation_receiver = Some(receiver);
                self.calculation_cancel = Some(cancel);
                self.calculation_state = CalculationState::Running { run_id };
            }
            Err(error) => {
                self.calculation_receiver = None;
                self.calculation_cancel = None;
                self.calculation_state = CalculationState::Failed;
                self.last_run_message =
                    Some(format!("Failed to start calculation worker: {}", error));
                self.last_run_is_error = true;
            }
        }
    }

    /// Start a deterministic worker used by lifecycle story tests.
    #[cfg(test)]
    pub(crate) fn start_test_calculation_worker(
        &mut self,
        delay: Duration,
        result: Result<(), String>,
    ) {
        self.start_worker_job(move |cancel| {
            let deadline = std::time::Instant::now() + delay;
            while std::time::Instant::now() < deadline {
                if cancel.load(Ordering::Acquire) {
                    return BvpWorkerOutcome::Cancelled;
                }
                std::thread::sleep(Duration::from_millis(2));
            }
            if cancel.load(Ordering::Acquire) {
                return BvpWorkerOutcome::Cancelled;
            }
            match result {
                Ok(()) => BvpWorkerOutcome::Completed { plot: None },
                Err(message) => BvpWorkerOutcome::Failed(message),
            }
        });
    }

    /// Request cancellation and discard the result once the backend returns.
    ///
    /// RustedSciThe does not currently expose an interrupt hook for an active
    /// damped solve. The flag therefore prevents pending work and discards a
    /// completed result, while the GUI remains responsive in `Cancelling`.
    pub(crate) fn cancel_calculation(&mut self) {
        let Some(run_id) = self.calculation_state.active_run_id() else {
            return;
        };
        if let Some(cancel) = &self.calculation_cancel {
            cancel.store(true, Ordering::Release);
        }
        self.calculation_state = CalculationState::Cancelling { run_id };
        self.last_run_message =
            Some("Cancellation requested; waiting for the solver backend to return.".to_string());
        self.last_run_is_error = false;
    }

    /// Poll the worker channel without blocking the egui thread.
    pub(crate) fn poll_calculation_worker(&mut self, ctx: &egui::Context) {
        let Some(receiver) = &self.calculation_receiver else {
            return;
        };

        match receiver.try_recv() {
            Ok(message) => {
                let active_run_id = self.calculation_state.active_run_id();
                self.calculation_receiver = None;
                self.calculation_cancel = None;
                if active_run_id == Some(message.run_id) {
                    self.apply_worker_outcome(message.outcome);
                }
                ctx.request_repaint();
            }
            Err(mpsc::TryRecvError::Empty) => {
                ctx.request_repaint_after(Duration::from_millis(100));
            }
            Err(mpsc::TryRecvError::Disconnected) => {
                self.calculation_receiver = None;
                self.calculation_cancel = None;
                self.calculation_state = CalculationState::Failed;
                self.last_run_message =
                    Some("Calculation worker disconnected without a result.".to_string());
                self.last_run_is_error = true;
                ctx.request_repaint();
            }
        }
    }

    /// Validate a run request and require explicit consent before invoking AOT tools.
    ///
    /// Lambdify requests proceed immediately. A new or changed AOT configuration
    /// is parked in `pending_aot_confirmation` and executed only after the user
    /// confirms the exact toolchain summary shown by the GUI.
    pub(crate) fn request_calculation(&mut self) {
        if let Err(message) = validate_aot_toolchain_fields(&self.document) {
            warn!("Rejected unsafe AOT GUI configuration: {}", message);
            self.pending_aot_confirmation = None;
            self.calculation_state = CalculationState::Failed;
            self.last_run_message = Some(format!("Calculation blocked: {}", message));
            self.last_run_is_error = true;
            return;
        }

        if let Some(request) = aot_toolchain_request(&self.document) {
            if self.confirmed_aot_request.as_ref() != Some(&request) {
                self.last_run_message =
                    Some("AOT toolchain confirmation is required before calculation.".to_string());
                self.last_run_is_error = false;
                self.pending_aot_confirmation = Some(request);
                return;
            }
        }

        self.pending_aot_confirmation = None;
        self.start_calculation_worker();
    }

    /// Render the modal trust boundary for external AOT toolchain execution.
    pub(crate) fn show_aot_confirmation(&mut self, ctx: &egui::Context) {
        let Some(request) = self.pending_aot_confirmation.clone() else {
            return;
        };

        let mut approve = false;
        let mut cancel = false;
        egui::Window::new("Confirm AOT toolchain execution")
            .collapsible(false)
            .resizable(false)
            .show(ctx, |ui| {
                ui.label("This calculation will invoke an external compiler toolchain.");
                ui.label(request.summary());
                ui.horizontal(|ui| {
                    if ui.button("Cancel AOT run").clicked() {
                        cancel = true;
                    }
                    if ui.button("Run AOT calculation").clicked() {
                        approve = true;
                    }
                });
            });

        if cancel {
            self.pending_aot_confirmation = None;
            self.last_run_message = Some("AOT calculation cancelled before execution.".to_string());
            self.last_run_is_error = false;
        } else if approve {
            self.confirmed_aot_request = Some(request);
            self.pending_aot_confirmation = None;
            self.start_calculation_worker();
        }
    }
}
