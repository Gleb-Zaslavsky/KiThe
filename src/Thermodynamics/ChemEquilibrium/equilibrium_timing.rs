//! Optional stage timing for the canonical equilibrium workflow.
//!
//! Timing is deliberately a value object rather than logging side effects.
//! A caller can enable it for characterization or diagnostics, inspect the
//! immutable report attached to a solution, and leave the default production
//! path free from per-stage clock reads.

use std::time::{Duration, Instant};

/// Controls whether the canonical workflow records stage timings.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum EquilibriumTimingMode {
    /// Do not read the clock or populate timing fields.
    #[default]
    Disabled,
    /// Record the requested stage durations in the immutable result report.
    Enabled,
}

/// Named stages used by the canonical equilibrium pipeline.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum EquilibriumTimingStage {
    RepositoryLookup,
    ThermochemistryPreparation,
    NumericClosureConstruction,
    SymbolicConstruction,
    EquationConstruction,
    NumericalProblemPreparation,
    ProjectionBuild,
    NonlinearSolve,
    PhaseControl,
    Validation,
    Postprocessing,
}

/// Immutable timing evidence attached to an accepted equilibrium result.
///
/// Durations are kept as `Duration` values so callers can choose their own
/// display unit. The `total` value is wall-clock time for the measured public
/// operation; stage values are inclusive where a stage contains nested solver
/// work and therefore should not be summed as an independent wall-clock total.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct EquilibriumTimingReport {
    enabled: bool,
    total: Duration,
    repository_lookup: Duration,
    thermochemistry_preparation: Duration,
    numeric_closure_construction: Duration,
    symbolic_construction: Duration,
    equation_construction: Duration,
    numerical_problem_preparation: Duration,
    projection_build: Duration,
    nonlinear_solve: Duration,
    phase_control: Duration,
    validation: Duration,
    postprocessing: Duration,
}

impl EquilibriumTimingReport {
    /// Whether this report contains measured values rather than zero defaults.
    pub fn enabled(&self) -> bool {
        self.enabled
    }

    /// Wall-clock time for the measured public operation.
    pub fn total(&self) -> Duration {
        self.total
    }

    pub fn repository_lookup(&self) -> Duration {
        self.repository_lookup
    }

    pub fn thermochemistry_preparation(&self) -> Duration {
        self.thermochemistry_preparation
    }

    pub fn numeric_closure_construction(&self) -> Duration {
        self.numeric_closure_construction
    }

    pub fn symbolic_construction(&self) -> Duration {
        self.symbolic_construction
    }

    pub fn equation_construction(&self) -> Duration {
        self.equation_construction
    }

    pub fn numerical_problem_preparation(&self) -> Duration {
        self.numerical_problem_preparation
    }

    pub fn projection_build(&self) -> Duration {
        self.projection_build
    }

    pub fn nonlinear_solve(&self) -> Duration {
        self.nonlinear_solve
    }

    pub fn phase_control(&self) -> Duration {
        self.phase_control
    }

    pub fn validation(&self) -> Duration {
        self.validation
    }

    pub fn postprocessing(&self) -> Duration {
        self.postprocessing
    }

    /// Adds one nested timing report to an outer aggregate.
    ///
    /// Stage values are inclusive, so callers must interpret the aggregate as
    /// accumulated evidence rather than summing it with the outer wall clock.
    pub(crate) fn accumulate(&mut self, nested: Self) {
        if !nested.enabled {
            return;
        }
        if !self.enabled {
            *self = nested;
            return;
        }
        self.total += nested.total;
        self.repository_lookup += nested.repository_lookup;
        self.thermochemistry_preparation += nested.thermochemistry_preparation;
        self.numeric_closure_construction += nested.numeric_closure_construction;
        self.symbolic_construction += nested.symbolic_construction;
        self.equation_construction += nested.equation_construction;
        self.numerical_problem_preparation += nested.numerical_problem_preparation;
        self.projection_build += nested.projection_build;
        self.nonlinear_solve += nested.nonlinear_solve;
        self.phase_control += nested.phase_control;
        self.validation += nested.validation;
        self.postprocessing += nested.postprocessing;
    }
}

/// Internal mutable collector used only while an immutable result is built.
#[derive(Debug)]
pub(crate) struct EquilibriumTimingCollector {
    report: EquilibriumTimingReport,
}

impl EquilibriumTimingCollector {
    pub(crate) fn new(mode: EquilibriumTimingMode) -> Self {
        Self {
            report: EquilibriumTimingReport {
                enabled: matches!(mode, EquilibriumTimingMode::Enabled),
                ..EquilibriumTimingReport::default()
            },
        }
    }

    pub(crate) fn from_report(report: EquilibriumTimingReport) -> Self {
        Self { report }
    }

    pub(crate) fn measure<T, F>(&mut self, stage: EquilibriumTimingStage, operation: F) -> T
    where
        F: FnOnce() -> T,
    {
        if !self.report.enabled {
            return operation();
        }

        let started = Instant::now();
        let result = operation();
        self.record(stage, started.elapsed());
        result
    }

    pub(crate) fn record(&mut self, stage: EquilibriumTimingStage, duration: Duration) {
        if !self.report.enabled {
            return;
        }

        let target = match stage {
            EquilibriumTimingStage::RepositoryLookup => &mut self.report.repository_lookup,
            EquilibriumTimingStage::ThermochemistryPreparation => {
                &mut self.report.thermochemistry_preparation
            }
            EquilibriumTimingStage::NumericClosureConstruction => {
                &mut self.report.numeric_closure_construction
            }
            EquilibriumTimingStage::SymbolicConstruction => &mut self.report.symbolic_construction,
            EquilibriumTimingStage::EquationConstruction => &mut self.report.equation_construction,
            EquilibriumTimingStage::NumericalProblemPreparation => {
                &mut self.report.numerical_problem_preparation
            }
            EquilibriumTimingStage::ProjectionBuild => &mut self.report.projection_build,
            EquilibriumTimingStage::NonlinearSolve => &mut self.report.nonlinear_solve,
            EquilibriumTimingStage::PhaseControl => &mut self.report.phase_control,
            EquilibriumTimingStage::Validation => &mut self.report.validation,
            EquilibriumTimingStage::Postprocessing => &mut self.report.postprocessing,
        };
        *target += duration;
    }

    pub(crate) fn set_total(&mut self, duration: Duration) {
        if self.report.enabled {
            self.report.total = duration;
        }
    }

    pub(crate) fn finish(self) -> EquilibriumTimingReport {
        self.report
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::thread;

    #[test]
    fn disabled_collection_does_not_publish_measurements() {
        let mut collector = EquilibriumTimingCollector::new(EquilibriumTimingMode::Disabled);
        collector.measure(EquilibriumTimingStage::NonlinearSolve, || {
            thread::yield_now()
        });
        let report = collector.finish();

        assert!(!report.enabled());
        assert_eq!(report.nonlinear_solve(), Duration::ZERO);
    }

    #[test]
    fn enabled_collection_keeps_stage_and_total_evidence() {
        let mut collector = EquilibriumTimingCollector::new(EquilibriumTimingMode::Enabled);
        collector.measure(EquilibriumTimingStage::SymbolicConstruction, || {
            thread::yield_now()
        });
        collector.set_total(Duration::from_nanos(7));
        let report = collector.finish();

        assert!(report.enabled());
        assert!(report.symbolic_construction() >= Duration::ZERO);
        assert_eq!(report.total(), Duration::from_nanos(7));
    }
}
