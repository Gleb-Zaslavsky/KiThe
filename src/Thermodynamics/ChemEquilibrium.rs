//! Chemical equilibrium and thermodynamics module index.
//!
//! The modern typed equilibrium stack lives alongside a small legacy surface
//! that is kept only while older call sites are being retired.

/// Legacy handwritten nonlinear solver implementation.
pub(crate) mod NR_Legacy;
/// Compatibility/experimental single-reaction helper retained during the
/// transition. It is not the production multiphase equilibrium API; use the
/// typed resolved-phase facade from [`prelude`] for new code.
pub(crate) mod easy_equilibrium;

/// Immutable global/local projection for one fixed active phase set.
pub(crate) mod equilibrium_active_set;
/// Shared activity-model contract for numerical and symbolic equilibrium paths.
pub mod equilibrium_activity;
/// Internal bridge between solver policy and concrete backend execution.
pub(crate) mod equilibrium_backend_adapter;
/// Deterministic element-to-record candidate selection and provenance reports.
pub mod equilibrium_candidate_selection;
/// Canonical phase-qualified component identity shared by bridge and solver.
pub mod equilibrium_component;
/// Cross-validation reports comparing canonical and independent K_eq results.
pub mod equilibrium_constant_cross_validation;
/// Independent reaction-extent/equilibrium-constant domain model.
pub mod equilibrium_constant_problem;
/// Independent small-system extent solver for equilibrium-constant validation.
pub mod equilibrium_constant_solver;
#[cfg(test)]
mod equilibrium_constant_solver_tests;
#[cfg(test)]
mod equilibrium_constant_tests;
/// Policies and reports for equilibrium-constant validation.
pub mod equilibrium_constant_validation;
#[cfg(test)]
mod equilibrium_golden_fixtures_tests;
/// Typed, validated input data for the canonical equilibrium formulation.
pub mod equilibrium_ids;
/// Temporary adapter around the historical hand-written nonlinear methods.
pub(crate) mod equilibrium_legacy_backend;
#[cfg(test)]
mod equilibrium_live_data_tests;
/// Canonical chemical-equilibrium solver using logarithmic species moles.
///
/// This is an implementation module. New consumers should use [`prelude`] or
/// the explicitly transitional [`legacy`] namespace instead of depending on
/// the mutable orchestration types directly.
pub(crate) mod equilibrium_log_moles;
#[cfg(test)]
mod equilibrium_log_moles_tests;
#[cfg(test)]
mod equilibrium_log_moles_tests2;
/// Canonical typed layout and physical initial composition for multiphase problems.
pub mod equilibrium_multiphase_domain;
#[cfg(test)]
mod equilibrium_multiphase_domain_tests;
#[cfg(test)]
mod equilibrium_multiphase_story_tests;
/// Reaction-basis construction, equilibrium errors, and temporary legacy backends.
pub mod equilibrium_nonlinear;
#[cfg(test)]
mod equilibrium_offline_regression_matrix_tests;
#[cfg(test)]
mod equilibrium_phase_bridge_tests;
/// Immutable numerical runner for one prepared fixed-`P,T` problem.
pub(crate) mod equilibrium_prepared_runner;
/// Typed, validated equilibrium problem construction and previews.
pub mod equilibrium_problem;
#[cfg(test)]
mod equilibrium_problem_tests;
#[cfg(test)]
mod equilibrium_public_api_tests;
/// Typed reaction-basis contract for independent validation.
pub mod equilibrium_reaction_basis;
/// RustedSciThe symbolic nonlinear-solver adapters for equilibrium.
pub mod equilibrium_rst_backend;
#[cfg(test)]
mod equilibrium_rst_backend_tests;
#[cfg(test)]
mod equilibrium_rst_matrix_tests;
/// Explicit nonlinear-backend policies and fallback diagnostics.
pub mod equilibrium_solver_policy;
/// Postprocessing for range-temperature equilibrium sweeps.
pub mod equilibrium_temperature_postprocessing;
/// Typed transactional fixed-pressure temperature continuation.
pub mod equilibrium_temperature_range;
/// Optional stage timing for the canonical equilibrium workflow.
pub mod equilibrium_timing;
/// Backend-independent validation for reconstructed equilibrium candidates.
pub mod equilibrium_validation;
#[cfg(test)]
mod equilibrium_validation_tests;
#[cfg(test)]
mod equilibrium_workflow_tests;
/// Problem construction, phase-control, and gas-equilibrium workflows.
///
/// The historical mutable orchestration remains inside the crate while
/// callers migrate to the typed resolved-phase workflow.
pub(crate) mod equilibrium_workflows;
/// Typed boundary from resolved phase data to equilibrium problem construction.
pub mod phase_equilibrium_problem;
#[cfg(test)]
mod phase_equilibrium_problem_tests;
/// Immutable phase-aware result assembled from one accepted bridge solve.
pub mod phase_equilibrium_solution;
/// Narrow fixed-pressure, fixed-temperature public workflow facade.
pub mod phase_equilibrium_workflow;
pub(crate) mod prepared_phase_control_runner;

/// Explicit compatibility namespace for the retained mutable equilibrium API.
///
/// The legacy implementations are still useful as policy-selected numerical
/// fallbacks and for old applications, but they are no longer presented as
/// peer production orchestration modules. New code should prefer [`prelude`].
pub mod legacy {
    /// Historical standalone Newton/LM/TR implementation.
    pub mod nr {
        pub use super::super::NR_Legacy::*;
    }

    /// Historical one-reaction convenience calculator.
    pub mod single_reaction {
        pub use super::super::easy_equilibrium::*;
    }

    pub use super::equilibrium_log_moles::*;
    pub use super::equilibrium_workflows::*;
}

/// Narrow re-export set for the production equilibrium path.
pub mod prelude {
    pub use crate::Thermodynamics::ChemEquilibrium::equilibrium_candidate_selection::{
        CandidateRejection, CandidateRejectionReason, CandidateSelectionError,
        CandidateTemperatureRange, CandidateTemperatureSupport, EquilibriumCandidate,
        EquilibriumCandidatePhaseAssignment, EquilibriumCandidatePhasePlan,
        EquilibriumCandidatePolicy, EquilibriumCandidateSelectionReport,
        EquilibriumCandidateSelector,
    };
    pub use crate::Thermodynamics::ChemEquilibrium::equilibrium_multiphase_domain::{
        MultiphaseEquilibriumLayout, MultiphaseInitialComposition,
    };
    pub use crate::Thermodynamics::ChemEquilibrium::equilibrium_nonlinear::{
        ReactionExtentError, ReactionExtentErrorKind,
    };
    pub use crate::Thermodynamics::ChemEquilibrium::equilibrium_problem::{
        EquilibriumConditions, TraceSpeciesSeedPolicy,
    };
    pub use crate::Thermodynamics::ChemEquilibrium::equilibrium_rst_backend::RustedSciTheSolver;
    pub use crate::Thermodynamics::ChemEquilibrium::equilibrium_solver_policy::{
        EquilibriumSolveReport, SolverAttemptMetrics, SolverAttemptOutcome, SolverAttemptReport,
        SolverBackend, SolverCascadeBudget, SolverPolicy,
    };
    pub use crate::Thermodynamics::ChemEquilibrium::equilibrium_temperature_postprocessing::{
        TemperatureInterpolationPolicy, TemperatureInterpolationSpace,
        TemperaturePostprocessingPolicy, TemperaturePostprocessingResult,
        TemperaturePostprocessingRow, TemperatureResamplingGrid, TemperatureSweepSeries,
        postprocess_temperature_range_solution, postprocess_temperature_series,
    };
    pub use crate::Thermodynamics::ChemEquilibrium::equilibrium_temperature_range::{
        TemperatureGrid, TemperatureRangeDirection, TemperatureRangeDurationSummary,
        TemperatureRangePoint, TemperatureRangePointPreparation, TemperatureRangePointReport,
        TemperatureRangeRequest, TemperatureRangeSolution, TemperatureRangeSolveReport,
    };
    pub use crate::Thermodynamics::ChemEquilibrium::equilibrium_timing::{
        EquilibriumTimingMode, EquilibriumTimingReport,
    };
    pub use crate::Thermodynamics::ChemEquilibrium::legacy::{
        InitialPhaseSet, Solvers as LegacyEquilibriumSolver,
    };
    pub use crate::Thermodynamics::ChemEquilibrium::phase_equilibrium_problem::{
        PhaseEquilibriumBuildReport, PhaseEquilibriumBuildRequest, PhaseEquilibriumProblemBundle,
        PhaseEquilibriumSolutionBundle, SupportedPhaseModelPolicy,
    };
    pub use crate::Thermodynamics::ChemEquilibrium::phase_equilibrium_solution::{
        MultiphaseEquilibriumSolution, MultiphaseEquilibriumSummaryRow,
    };
    pub use crate::Thermodynamics::ChemEquilibrium::phase_equilibrium_workflow::{
        EquilibriumSolveOptions, PhaseControlPolicy, PhaseEquilibriumPipelineError,
        PhaseEquilibriumPipelineRequest, PhaseEquilibriumSolveMode,
        ResolvedPhaseEquilibriumOutcome, ResolvedPhaseEquilibriumRequest, solve_resolved_pt,
    };
    pub use crate::Thermodynamics::User_PhaseOrSolution::{
        PhaseModel, PhaseSpec, ResolvedPhaseSystem, ResolvedPhaseSystemReport,
        SubstanceSystemFactory, SubstanceSystemFactoryError, SubstanceSystemSpec,
        SubstanceSystemSpecBuilder, SubstancesContainer,
    };
    pub use crate::Thermodynamics::phase_layout::{PhaseComponentId, PhaseId};
    pub use crate::Thermodynamics::physical_state::PhysicalState;
    pub use crate::Thermodynamics::thermo_lib_api::{ElementSearchMode, ThermoRepository};
}
