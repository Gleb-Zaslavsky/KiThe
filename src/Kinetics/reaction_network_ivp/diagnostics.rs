use std::time::Duration;

use RustedSciThe::numerical::SolverTelemetryMode;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NetworkIvpStage {
    EquationAssembly,
    BackendPreparation,
    Integration,
    OutputPostprocessing,
}

/// Optional stage values follow the requested mode: counts are absent in
/// Timings, durations are absent in Counters, and the whole snapshot is absent
/// in Off. Durations are cumulative for a prepared session.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct NetworkIvpStageMeasurement {
    pub calls: Option<u64>,
    pub elapsed: Option<Duration>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct NetworkIvpStageDiagnostics {
    pub mode: SolverTelemetryMode,
    pub equation_assembly: NetworkIvpStageMeasurement,
    pub backend_preparation: NetworkIvpStageMeasurement,
    pub integration: NetworkIvpStageMeasurement,
    pub output_postprocessing: NetworkIvpStageMeasurement,
}

impl NetworkIvpStageDiagnostics {
    pub(crate) fn new(mode: SolverTelemetryMode) -> Option<Self> {
        (mode != SolverTelemetryMode::Off).then_some(Self {
            mode,
            equation_assembly: NetworkIvpStageMeasurement::default(),
            backend_preparation: NetworkIvpStageMeasurement::default(),
            integration: NetworkIvpStageMeasurement::default(),
            output_postprocessing: NetworkIvpStageMeasurement::default(),
        })
    }

    pub(crate) fn records_timings(&self) -> bool {
        matches!(
            self.mode,
            SolverTelemetryMode::Timings | SolverTelemetryMode::Detailed
        )
    }

    pub(crate) fn record(&mut self, stage: NetworkIvpStage, elapsed: Option<Duration>) {
        let measurement = match stage {
            NetworkIvpStage::EquationAssembly => &mut self.equation_assembly,
            NetworkIvpStage::BackendPreparation => &mut self.backend_preparation,
            NetworkIvpStage::Integration => &mut self.integration,
            NetworkIvpStage::OutputPostprocessing => &mut self.output_postprocessing,
        };
        if matches!(
            self.mode,
            SolverTelemetryMode::Counters | SolverTelemetryMode::Detailed
        ) {
            *measurement.calls.get_or_insert(0) += 1;
        }
        if let Some(elapsed) = elapsed {
            *measurement.elapsed.get_or_insert(Duration::ZERO) += elapsed;
        }
    }
}
