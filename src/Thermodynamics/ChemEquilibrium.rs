//! Chemical equilibrium and thermodynamics module index.
//!
//! The modern typed equilibrium stack lives alongside a small legacy surface
//! that is kept only while older call sites are being retired.
/// Legacy handwritten nonlinear solver implementation.
#[path = "ChemEquilibrium/nonlinear_solvers/NR_Legacy.rs"]
pub(crate) mod NR_Legacy;
/// Compatibility/experimental single-reaction helper retained during the
/// transition. It is not the production multiphase equilibrium API; use the
/// typed resolved-phase facade from [`prelude`] for new code.
pub(crate) mod easy_equilibrium;
/// Read-only typed infrastructure for frozen external reference evidence.
///
/// This is test-only: production thermochemistry continues to enter through
/// the repository and format-handler stack.
#[cfg(test)]
#[path = "ChemEquilibrium/frozen_reference/mod.rs"]
mod frozen_reference;
#[cfg(test)]
#[path = "ChemEquilibrium/frozen_reference/ar_water_co2_sublimation_gauge.rs"]
mod frozen_reference_ar_water_co2_sublimation_gauge;
#[cfg(test)]
#[path = "ChemEquilibrium/frozen_reference/ar_water_co2_three_phase_preflight.rs"]
mod frozen_reference_ar_water_co2_three_phase_preflight;
/// First general multicomponent fixed-`P,T` I5 characterization against a
/// frozen Argonne/STANJAN CHON equilibrium table.
#[cfg(test)]
#[path = "ChemEquilibrium/frozen_reference/argonne_stanjan_chon.rs"]
mod frozen_reference_argonne_stanjan_chon;
#[cfg(test)]
#[path = "ChemEquilibrium/frozen_reference/argonne_stanjan_chon_tests.rs"]
mod frozen_reference_argonne_stanjan_chon_tests;
/// Test-only feasibility preflight for competing liquid/ice water candidates.
#[cfg(test)]
#[path = "ChemEquilibrium/frozen_reference/iapws_competing_water_preflight.rs"]
mod frozen_reference_iapws_competing_water_preflight;
/// Test-only competing liquid/ice candidate gauge and phase-order matrix.
#[cfg(test)]
#[path = "ChemEquilibrium/frozen_reference/iapws_exclusive_competing_candidates.rs"]
mod frozen_reference_iapws_exclusive_competing_candidates;
/// Test-only log-moles dynamic-range and vanishing-species evidence.
#[cfg(test)]
#[path = "ChemEquilibrium/frozen_reference/extreme_dynamic_range.rs"]
mod frozen_reference_extreme_dynamic_range;
#[cfg(test)]
#[path = "ChemEquilibrium/frozen_reference/iapws_ice_sublimation.rs"]
mod frozen_reference_iapws_ice_sublimation;
#[cfg(test)]
#[path = "ChemEquilibrium/frozen_reference/iapws_ice_tests.rs"]
mod frozen_reference_iapws_ice_tests;
#[cfg(test)]
#[path = "ChemEquilibrium/frozen_reference/iapws_tests.rs"]
mod frozen_reference_iapws_tests;
#[cfg(test)]
#[path = "ChemEquilibrium/frozen_reference/iapws_water_saturation.rs"]
mod frozen_reference_iapws_water_saturation;
#[cfg(test)]
#[path = "ChemEquilibrium/frozen_reference/janaf_boudouard_boundary.rs"]
mod frozen_reference_janaf_boudouard_boundary;
#[cfg(test)]
#[path = "ChemEquilibrium/frozen_reference/janaf_boudouard_boundary_tests.rs"]
mod frozen_reference_janaf_boudouard_boundary_tests;
#[cfg(test)]
#[path = "ChemEquilibrium/frozen_reference/janaf_boudouard_tests.rs"]
mod frozen_reference_janaf_boudouard_tests;
#[cfg(test)]
#[path = "ChemEquilibrium/frozen_reference/janaf_boudouard_thermochemistry.rs"]
mod frozen_reference_janaf_boudouard_thermochemistry;
#[cfg(test)]
#[path = "ChemEquilibrium/frozen_reference/nasa_cea_h2_o2_hp.rs"]
mod frozen_reference_nasa_cea_h2_o2_hp;
#[cfg(test)]
#[path = "ChemEquilibrium/frozen_reference/nasa_cea_h2_o2_hp_tests.rs"]
mod frozen_reference_nasa_cea_h2_o2_hp_tests;
#[cfg(test)]
#[path = "ChemEquilibrium/frozen_reference/nasa_cea_rp1311_example3.rs"]
mod frozen_reference_nasa_cea_rp1311_example3;
/// NASA TP-1906 specific-enthalpy evidence joined semantically to the existing
/// TP-1907 CHON + graphite composition fixture for the external P,H story.
#[cfg(test)]
#[path = "ChemEquilibrium/frozen_reference/nasa_tp1906_chon_graphite_ph.rs"]
mod frozen_reference_nasa_tp1906_chon_graphite_ph;
#[cfg(test)]
#[path = "ChemEquilibrium/frozen_reference/nasa_tp1906_chon_graphite_ph_tests.rs"]
mod frozen_reference_nasa_tp1906_chon_graphite_ph_tests;
/// Real-data feasibility matrix for a continuation failure after acceptance.
#[cfg(test)]
#[path = "ChemEquilibrium/frozen_reference/nasa_tp1906_chon_graphite_rollback_tests.rs"]
mod frozen_reference_nasa_tp1906_chon_graphite_rollback_tests;
/// NASA TP-1907 heterogeneous CHON + graphite I5 fixture and diagnostics.
#[cfg(test)]
#[path = "ChemEquilibrium/frozen_reference/nasa_tp1907_chon_graphite.rs"]
mod frozen_reference_nasa_tp1907_chon_graphite;
#[cfg(test)]
#[path = "ChemEquilibrium/frozen_reference/nasa_tp1907_chon_graphite_tests.rs"]
mod frozen_reference_nasa_tp1907_chon_graphite_tests;
#[cfg(test)]
#[path = "ChemEquilibrium/frozen_reference/nist_tashkun_harvey_co2_tests.rs"]
mod frozen_reference_nist_tashkun_harvey_co2_tests;

/// Test-only element-conserving interior seed diagnostics for fresh-basin
/// investigations. This is deliberately not production recovery policy.
#[cfg(test)]
#[path = "ChemEquilibrium/frozen_reference/interior_seed.rs"]
mod frozen_reference_interior_seed;

/// Repository-level catalog checks for frozen metadata/rows pairs. Typed
/// physical adapters remain separate from this test-only integrity layer.
#[cfg(test)]
#[path = "ChemEquilibrium/frozen_reference/catalog_tests.rs"]
mod frozen_reference_catalog_tests;
#[cfg(test)]
#[path = "ChemEquilibrium/frozen_reference/nist_benzene_toluene_vle.rs"]
mod frozen_reference_nist_benzene_toluene_vle;
#[cfg(test)]
#[path = "ChemEquilibrium/frozen_reference/nist_benzene_toluene_vle_tests.rs"]
mod frozen_reference_nist_benzene_toluene_vle_tests;
/// Test-only NIST ethylbenzene anchors and independent pure-component oracles.
#[cfg(test)]
#[path = "ChemEquilibrium/frozen_reference/nist_ethylbenzene_pure_component.rs"]
mod frozen_reference_nist_ethylbenzene_pure_component;
#[cfg(test)]
#[path = "ChemEquilibrium/frozen_reference/nist_ethylbenzene_pure_component_tests.rs"]
mod frozen_reference_nist_ethylbenzene_pure_component_tests;
/// Test-only three-component Antoine gauge for the ternary `P,T` VLE story.
#[cfg(test)]
#[path = "ChemEquilibrium/frozen_reference/nist_ternary_vle_antoine_gauge.rs"]
mod frozen_reference_nist_ternary_vle_antoine_gauge;
#[cfg(test)]
#[path = "ChemEquilibrium/frozen_reference/nist_ternary_vle_antoine_gauge_tests.rs"]
mod frozen_reference_nist_ternary_vle_antoine_gauge_tests;
/// Seven-point fixed-inventory forward/reverse lifecycle over the ternary
/// Antoine gauge. Kept separate from the focused activation smoke stories.
#[cfg(test)]
#[path = "ChemEquilibrium/frozen_reference/nist_ternary_vle_full_lifecycle_tests.rs"]
mod frozen_reference_nist_ternary_vle_full_lifecycle_tests;
/// Full bounded phase-control lifecycle for the ternary Antoine-gauge story.
#[cfg(test)]
#[path = "ChemEquilibrium/frozen_reference/nist_ternary_vle_lifecycle_tests.rs"]
mod frozen_reference_nist_ternary_vle_lifecycle_tests;
/// Capability/source audit for the proposed ternary NIST ThermoML VLE fixture.
#[cfg(test)]
#[path = "ChemEquilibrium/frozen_reference/nist_toluene_ethylbenzene_chlorobenzene_preflight.rs"]
mod frozen_reference_nist_toluene_ethylbenzene_chlorobenzene_preflight;
#[cfg(test)]
#[path = "ChemEquilibrium/frozen_reference/nist_toluene_ethylbenzene_chlorobenzene_preflight_tests.rs"]
mod frozen_reference_nist_toluene_ethylbenzene_chlorobenzene_preflight_tests;
#[cfg(test)]
#[path = "ChemEquilibrium/frozen_reference/loader_tests.rs"]
mod frozen_reference_tests;
// Cross-validation source files are grouped under `cross_validation/`, while
// these stable module names avoid unnecessary churn for internal consumers and
// existing test filters.
/// Cross-validation stories for pure-phase boundary mathematics and lifecycle.
#[cfg(test)]
#[path = "ChemEquilibrium/cross_validation/phase_boundary_cross_validation_tests.rs"]
mod phase_boundary_cross_validation_tests;
/// Adapters that derive pure-phase cross-validation evidence from immutable
/// production phase-control results.
#[path = "ChemEquilibrium/cross_validation/phase_boundary_production_adapter.rs"]
pub mod phase_boundary_production_adapter;
#[path = "ChemEquilibrium/cross_validation/phase_boundary_validation.rs"]
pub mod phase_boundary_validation;
#[cfg(test)]
#[path = "ChemEquilibrium/cross_validation/phase_boundary_validation_tests.rs"]
mod phase_boundary_validation_tests;
#[cfg(test)]
#[path = "ChemEquilibrium/cross_validation/pure_phase_ph_lifecycle_tests.rs"]
mod pure_phase_ph_lifecycle_tests;
/// Offline real-data I4 validation for the independent pure-phase `P,H` route.
#[cfg(test)]
#[path = "ChemEquilibrium/cross_validation/pure_phase_ph_live_data_tests.rs"]
mod pure_phase_ph_live_data_tests;
/// Independent scalar `P,H` validation for one gas phase and one pure
/// condensed candidate. This remains crate-private until its I1-I3 evidence
/// ladder is complete; it is not a second production P,H workflow.
#[cfg(test)]
#[path = "ChemEquilibrium/cross_validation/pure_phase_ph_validation.rs"]
pub(crate) mod pure_phase_ph_validation;
#[cfg(test)]
#[path = "ChemEquilibrium/cross_validation/pure_phase_ph_validation_tests.rs"]
mod pure_phase_ph_validation_tests;
#[cfg(test)]
#[path = "ChemEquilibrium/cross_validation/pure_phase_pt_live_data_tests.rs"]
mod pure_phase_pt_live_data_tests;
/// Real Boudouard `P,H` independent/canonical validation stories.
#[cfg(test)]
#[path = "ChemEquilibrium/cross_validation/real_boudouard_ph_tests.rs"]
mod real_boudouard_ph_tests;
/// Test-side real Boudouard `P,H` reference-state builder.
#[cfg(test)]
#[path = "ChemEquilibrium/cross_validation/real_boudouard_ph_validation.rs"]
mod real_boudouard_ph_validation;
/// Shared, offline real-data chemistry definitions for pure-phase validation.
/// This remains test-only until its inventory contract is established.
#[cfg(test)]
#[path = "ChemEquilibrium/cross_validation/real_pure_phase_fixtures.rs"]
mod real_pure_phase_fixtures;

/// Immutable global/local projection for one fixed active phase set.
#[path = "ChemEquilibrium/phase_control/equilibrium_active_set.rs"]
pub(crate) mod equilibrium_active_set;
/// Shared activity-model contract for numerical and symbolic equilibrium paths.
pub mod equilibrium_activity;
/// Internal bridge between solver policy and concrete backend execution.
#[path = "ChemEquilibrium/nonlinear_solvers/equilibrium_backend_adapter.rs"]
pub(crate) mod equilibrium_backend_adapter;
/// Deterministic element-to-record candidate selection and provenance reports.
pub mod equilibrium_candidate_selection;
/// Strict read-only comparison of two accepted equilibrium solutions.
pub mod equilibrium_comparison;
/// Canonical phase-qualified component identity shared by bridge and solver.
pub mod equilibrium_component;
/// Cross-validation reports comparing canonical and independent K_eq results.
#[path = "ChemEquilibrium/equilibrium_constants/equilibrium_constant_cross_validation.rs"]
pub mod equilibrium_constant_cross_validation;
/// Independent reaction-extent/equilibrium-constant domain model.
#[path = "ChemEquilibrium/equilibrium_constants/equilibrium_constant_problem.rs"]
pub mod equilibrium_constant_problem;
/// Independent small-system extent solver for equilibrium-constant validation.
#[path = "ChemEquilibrium/equilibrium_constants/equilibrium_constant_solver.rs"]
pub mod equilibrium_constant_solver;
#[cfg(test)]
#[path = "ChemEquilibrium/equilibrium_constants/equilibrium_constant_solver_tests.rs"]
mod equilibrium_constant_solver_tests;
#[cfg(test)]
#[path = "ChemEquilibrium/equilibrium_constants/equilibrium_constant_tests.rs"]
mod equilibrium_constant_tests;
/// Policies and reports for equilibrium-constant validation.
#[path = "ChemEquilibrium/equilibrium_constants/equilibrium_constant_validation.rs"]
pub mod equilibrium_constant_validation;
/// Typed PT/PH constraints and pure enthalpy-domain helpers.
#[path = "ChemEquilibrium/ph/equilibrium_constraints.rs"]
pub mod equilibrium_constraints;
/// Optional structured lifecycle diagnostics and live diagnostic sinks.
#[path = "ChemEquilibrium/postprocessing_and_logging/equilibrium_diagnostics.rs"]
pub mod equilibrium_diagnostics;
/// Human-readable formatting and explicit log-facade output for diagnostics.
#[path = "ChemEquilibrium/postprocessing_and_logging/equilibrium_diagnostics_display.rs"]
pub mod equilibrium_diagnostics_display;
/// Presentation-only display filtering, units, and numeric formatting.
#[path = "ChemEquilibrium/postprocessing_and_logging/equilibrium_display.rs"]
pub mod equilibrium_display;
/// Cooperative cancellation and progress events for typed workflows.
pub mod equilibrium_execution;
/// Explicit physical-to-normalized representation mapping for extensive data.
///
/// This module defines units and conversions only; it never changes the
/// canonical solver route or retry policy by itself.
pub mod equilibrium_extensive_normalization;
#[cfg(test)]
mod equilibrium_golden_fixtures_tests;
/// Typed, validated input data for the canonical equilibrium formulation.
pub mod equilibrium_ids;
/// Temporary adapter around the historical hand-written nonlinear methods.
#[path = "ChemEquilibrium/nonlinear_solvers/equilibrium_legacy_backend.rs"]
pub(crate) mod equilibrium_legacy_backend;
#[cfg(test)]
#[path = "ChemEquilibrium/test_suites/live_data/equilibrium_live_data_tests.rs"]
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
#[path = "ChemEquilibrium/phase_bridge/equilibrium_multiphase_story_tests.rs"]
mod equilibrium_multiphase_story_tests;
/// Reaction-basis construction, equilibrium errors, and temporary legacy backends.
#[path = "ChemEquilibrium/nonlinear_solvers/equilibrium_nonlinear.rs"]
pub mod equilibrium_nonlinear;
#[cfg(test)]
mod equilibrium_offline_regression_matrix_tests;
/// Structural contracts for the monolithic fixed-P,H formulation.
#[path = "ChemEquilibrium/ph/equilibrium_ph_formulation.rs"]
pub(crate) mod equilibrium_ph_formulation;
/// Fixed-active-set numerical runner for the monolithic P,H formulation.
#[path = "ChemEquilibrium/ph/equilibrium_ph_monolithic.rs"]
pub(crate) mod equilibrium_ph_monolithic;
/// Nested scalar P,H reports and the stateless bracket engine. The stable
/// workflow facade re-exports the report types; callers normally do not need
/// to address this module directly.
#[path = "ChemEquilibrium/ph/equilibrium_ph_nested.rs"]
pub mod equilibrium_ph_nested;
/// Shared and route-specific typed controls for fixed-P,H solving.
#[path = "ChemEquilibrium/ph/equilibrium_ph_options.rs"]
pub mod equilibrium_ph_options;
/// Transactional fixed-pressure target-enthalpy continuation.
#[path = "ChemEquilibrium/ph/equilibrium_ph_range.rs"]
pub mod equilibrium_ph_range;
/// Read-only raw series and diagnostics for fixed-pressure enthalpy ranges.
#[path = "ChemEquilibrium/ph/equilibrium_ph_range_presentation.rs"]
pub mod equilibrium_ph_range_presentation;
/// Typed thermochemistry capabilities shared by P,H formulations and runners.
#[path = "ChemEquilibrium/ph/equilibrium_ph_thermochemistry.rs"]
pub(crate) mod equilibrium_ph_thermochemistry;
/// Outer fixed-pressure, fixed-total-enthalpy workflow over the canonical PT solver.
#[path = "ChemEquilibrium/ph/equilibrium_ph_workflow.rs"]
pub mod equilibrium_ph_workflow;
#[cfg(test)]
#[path = "ChemEquilibrium/phase_bridge/equilibrium_phase_bridge_tests.rs"]
mod equilibrium_phase_bridge_tests;
/// Pure canonical-state and element-potential primitives for phase stability.
#[path = "ChemEquilibrium/phase_control/equilibrium_phase_stability.rs"]
pub(crate) mod equilibrium_phase_stability;
/// Immutable numerical runner for one prepared fixed-`P,T` problem.
pub(crate) mod equilibrium_prepared_runner;
/// Read-only rows and compact rendering for accepted equilibrium evidence.
#[path = "ChemEquilibrium/postprocessing_and_logging/equilibrium_presentation.rs"]
pub mod equilibrium_presentation;
/// Typed, validated equilibrium problem construction and previews.
pub mod equilibrium_problem;
#[cfg(test)]
mod equilibrium_problem_tests;
#[cfg(test)]
mod equilibrium_public_api_tests;
/// Read-only series and transition boundaries for fixed-pressure temperature ranges.
#[path = "ChemEquilibrium/postprocessing_and_logging/equilibrium_range_presentation.rs"]
pub mod equilibrium_range_presentation;
/// Typed reaction-basis contract for independent validation.
pub mod equilibrium_reaction_basis;
/// Exportable immutable provenance and effective-policy capsule for accepted runs.
pub mod equilibrium_reproducibility;
/// RustedSciThe symbolic nonlinear-solver adapters for equilibrium.
#[path = "ChemEquilibrium/nonlinear_solvers/equilibrium_rst_backend.rs"]
pub mod equilibrium_rst_backend;
#[cfg(test)]
#[path = "ChemEquilibrium/nonlinear_solvers/equilibrium_rst_backend_tests.rs"]
mod equilibrium_rst_backend_tests;
#[cfg(test)]
#[path = "ChemEquilibrium/nonlinear_solvers/equilibrium_rst_matrix_tests.rs"]
mod equilibrium_rst_matrix_tests;
/// Explicit nonlinear-backend policies and fallback diagnostics.
#[path = "ChemEquilibrium/nonlinear_solvers/equilibrium_solver_policy.rs"]
pub mod equilibrium_solver_policy;
/// Postprocessing for range-temperature equilibrium sweeps.
#[path = "ChemEquilibrium/postprocessing_and_logging/equilibrium_temperature_postprocessing.rs"]
pub mod equilibrium_temperature_postprocessing;
/// Typed transactional fixed-pressure temperature continuation.
pub mod equilibrium_temperature_range;
/// Optional stage timing for the canonical equilibrium workflow.
#[path = "ChemEquilibrium/postprocessing_and_logging/equilibrium_timing.rs"]
pub mod equilibrium_timing;
/// Backend-independent validation for reconstructed equilibrium candidates.
#[path = "ChemEquilibrium/nonlinear_solvers/equilibrium_validation.rs"]
pub mod equilibrium_validation;
#[cfg(test)]
#[path = "ChemEquilibrium/nonlinear_solvers/equilibrium_validation_tests.rs"]
mod equilibrium_validation_tests;
#[cfg(test)]
#[path = "ChemEquilibrium/phase_control/equilibrium_workflow_tests.rs"]
mod equilibrium_workflow_tests;
/// Problem construction, phase-control, and gas-equilibrium workflows.
///
/// The historical mutable orchestration remains inside the crate while
/// callers migrate to the typed resolved-phase workflow.
#[path = "ChemEquilibrium/phase_control/equilibrium_workflows.rs"]
pub(crate) mod equilibrium_workflows;
/// Typed boundary from resolved phase data to equilibrium problem construction.
#[path = "ChemEquilibrium/phase_bridge/phase_equilibrium_problem.rs"]
pub mod phase_equilibrium_problem;
#[cfg(test)]
#[path = "ChemEquilibrium/phase_bridge/phase_equilibrium_problem_tests.rs"]
mod phase_equilibrium_problem_tests;
/// Immutable phase-aware result assembled from one accepted bridge solve.
#[path = "ChemEquilibrium/phase_bridge/phase_equilibrium_solution.rs"]
pub mod phase_equilibrium_solution;
/// Narrow fixed-pressure, fixed-temperature public workflow facade.
#[path = "ChemEquilibrium/phase_bridge/phase_equilibrium_workflow.rs"]
pub mod phase_equilibrium_workflow;
#[path = "ChemEquilibrium/phase_control/prepared_phase_control_runner.rs"]
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
    pub use crate::Thermodynamics::ChemEquilibrium::equilibrium_comparison::{
        EquilibriumComparisonReport, EquilibriumComparisonSummary,
        EquilibriumComponentComparisonRow, EquilibriumPhaseComparisonRow,
    };
    pub use crate::Thermodynamics::ChemEquilibrium::equilibrium_constant_validation::EquilibriumConstantValidationMode;
    pub use crate::Thermodynamics::ChemEquilibrium::equilibrium_constraints::{
        EnthalpyScale, EquilibriumConstraint, TemperatureBounds, TotalEnthalpyJoules,
        additive_total_enthalpy,
    };
    pub use crate::Thermodynamics::ChemEquilibrium::equilibrium_diagnostics::{
        EquilibriumDiagnosticEvent, EquilibriumDiagnosticSink, EquilibriumDiagnosticsMode,
        EquilibriumDiagnosticsOptions, EquilibriumDiagnosticsReport,
        EquilibriumRangeDiagnosticsPolicy, PhDiagnosticRoute, PhaseStabilityDiagnostic,
    };
    pub use crate::Thermodynamics::ChemEquilibrium::equilibrium_diagnostics_display::{
        format_diagnostics, format_ph_solution_execution_summary,
        format_solution_diagnostics, format_solution_execution_summary,
        log_ph_solution_execution_summary, log_solution_diagnostics,
        log_solution_execution_summary,
    };
    pub use crate::Thermodynamics::ChemEquilibrium::equilibrium_display::{
        DisplayPolicyError, EquilibriumDisplayPolicy, EquilibriumFormattedComponentRow,
        EquilibriumFormattedPhaseRow, EquilibriumNumberStyle, MoleFractionDisplay,
    };
    pub use crate::Thermodynamics::ChemEquilibrium::equilibrium_execution::{
        EquilibriumExecutionControl, EquilibriumProgressEvent, EquilibriumProgressStage,
    };
    pub use crate::Thermodynamics::ChemEquilibrium::equilibrium_extensive_normalization::ExtensiveNormalization;
    pub use crate::Thermodynamics::ChemEquilibrium::equilibrium_ids::{
        ElementId, PhaseIndex, ReactionId, SpeciesId,
    };
    pub use crate::Thermodynamics::ChemEquilibrium::equilibrium_multiphase_domain::{
        MultiphaseEquilibriumLayout, MultiphaseInitialComposition,
    };
    pub use crate::Thermodynamics::ChemEquilibrium::equilibrium_nonlinear::{
        PhMonolithicSeedFailure, ReactionExtentError, ReactionExtentErrorKind,
    };
    pub use crate::Thermodynamics::ChemEquilibrium::equilibrium_ph_options::{
        PhAcceptanceOptions, PhMonolithicOptions, PhMonolithicSeedAttemptReport,
        PhMonolithicSeedPolicy, PhMonolithicSeedRecoveryReport, PhNestedOptions,
    };
    pub use crate::Thermodynamics::ChemEquilibrium::equilibrium_ph_range::{
        PhEnthalpyGrid, PhRangeDirection, PhRangeDurationSummary, PhRangeError, PhRangePoint,
        PhRangePointError, PhRangePointPreparation, PhRangePointReport, PhRangeRequest,
        PhRangeSolution, PhRangeSolveReport,
    };
    pub use crate::Thermodynamics::ChemEquilibrium::equilibrium_ph_range_presentation::{
        EnthalpySweepSeries, PhRangePresentationPointRow, PhRangePresentationReport,
    };
    pub use crate::Thermodynamics::ChemEquilibrium::equilibrium_ph_workflow::{
        EnthalpyEvaluation, EnthalpyModel, FixedPressureEnthalpySolution, MolarEnthalpyFunction,
        MolarThermoFunction, PhFallbackReason, PhMonolithicEvidence, PhMonotonicityPolicy,
        PhRouteDecision, PhSolveMode, PhSolvePath, PhTemperatureSolveOptions,
        PhTemperatureSolveReport, PhTemperatureStepKind, PhTemperatureTimingReport,
        PhTemperatureTrial, PhTrialInnerEvidence, PhTrialPhaseState, PhTrialPreparation,
        PhTrialTimingReport, ResolvedPhaseEnthalpyRequest, ResolvedThermochemistry,
        ThermochemistryProvenance, ThermochemistryStandardStatePressure, solve_resolved_ph,
    };
    pub use crate::Thermodynamics::ChemEquilibrium::equilibrium_presentation::{
        EquilibriumBackendAttemptPresentationRow, EquilibriumComponentPresentationRow,
        EquilibriumPhasePresentationRow, EquilibriumPresentationReport,
        EquilibriumPresentationSummary, EquilibriumTimingPresentationRow, backend_attempt_rows,
        backend_attempt_rows_from_error,
    };
    pub use crate::Thermodynamics::ChemEquilibrium::equilibrium_problem::{
        EquilibriumConditions, LogMolesInitialGuess, TraceSpeciesSeedPolicy,
    };
    pub use crate::Thermodynamics::ChemEquilibrium::equilibrium_range_presentation::{
        TemperatureRangePresentationPointRow, TemperatureRangePresentationReport,
        TemperatureRangeTransitionBoundary,
    };
    pub use crate::Thermodynamics::ChemEquilibrium::equilibrium_reproducibility::{
        EQUILIBRIUM_REPRODUCIBILITY_SCHEMA_VERSION, EquilibriumCandidateRecordSnapshot,
        EquilibriumCandidateSelectionSnapshot, EquilibriumPhaseSpecSnapshot,
        EquilibriumRecordIdentity, EquilibriumReproducibilityCapsule, PhaseStabilitySemantics,
        ReproducibilityCapsuleError, ThermoCatalogSnapshot,
    };
    pub use crate::Thermodynamics::ChemEquilibrium::equilibrium_rst_backend::RustedSciTheSolver;
    pub use crate::Thermodynamics::ChemEquilibrium::equilibrium_solver_policy::{
        EquilibriumSolveReport, MultiStartAttemptReport, MultiStartSolveReport,
        SolverAttemptMetrics, SolverAttemptOutcome, SolverAttemptReport, SolverBackend,
        SolverCascadeBudget, SolverPolicy,
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
    pub use crate::Thermodynamics::ChemEquilibrium::equilibrium_validation::EquilibriumCandidateReport;
    /// Immutable phase-control and TPD evidence emitted by accepted solves.
    ///
    /// These types are re-exported here instead of exposing the mutable
    /// orchestration module that constructs them. Consumers can inspect the
    /// accepted phase lifecycle without acquiring a second solve API.
    pub use crate::Thermodynamics::ChemEquilibrium::equilibrium_workflows::{
        ElementPotentialReport, ElementalFeasibilityReport, MultiphaseAcceptanceReport,
        MultiphaseAcceptanceRow, MultiphaseComplementarityReport, PhaseControlledSolveReport,
        PhaseControlledSolveRow, PhaseSet, PhaseStabilityConditions, PhaseStabilityLayout,
        PhaseStabilityReport, PhaseStabilityStatus, PhaseStatus, PhaseTransitionReason,
        PhaseTransitionRecord, TpdMinimizerReport,
    };
    pub use crate::Thermodynamics::ChemEquilibrium::legacy::{
        InitialPhaseSet, Solvers as LegacyEquilibriumSolver,
    };
    pub use crate::Thermodynamics::ChemEquilibrium::phase_equilibrium_problem::{
        PhaseEquilibriumBuildReport, PhaseEquilibriumBuildRequest, PhaseEquilibriumProblemBundle,
        PhaseEquilibriumSolutionBundle, SupportedPhaseModelPolicy,
    };
    pub use crate::Thermodynamics::ChemEquilibrium::phase_equilibrium_solution::{
        ExtensiveNormalizationRecoveryEvidence, MultiphaseEquilibriumSolution,
        MultiphaseEquilibriumSummaryRow,
    };
    pub use crate::Thermodynamics::ChemEquilibrium::phase_equilibrium_workflow::{
        EquilibriumSolveOptions, EquilibriumSolveOptionsSnapshot, EquilibriumSolverBudgetSnapshot,
        ExtensiveNormalizationPolicy, PhaseControlPolicy, PhaseEquilibriumPipelineError,
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
    pub use crate::Thermodynamics::thermo_lib_api::{
        ElementSearchMode, ThermoCatalogConsistencyReport, ThermoRepository,
    };
}
