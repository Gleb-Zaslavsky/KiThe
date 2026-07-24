//! Chemical equilibrium and thermodynamics module index.
//!
//! The modern typed equilibrium stack lives alongside a small legacy surface
//! that is kept only while older call sites are being retired.

/// Legacy handwritten nonlinear solver.
pub mod NR_Legacy;
/// Small helper examples retained for compatibility during the transition.
pub mod easy_equilibrium;

/// Immutable global/local projection for one fixed active phase set.
pub(crate) mod equilibrium_active_set;
/// Shared activity-model contract for numerical and symbolic equilibrium paths.
pub mod equilibrium_activity;
/// Internal bridge between solver policy and concrete backend execution.
pub(crate) mod equilibrium_backend_adapter;
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
/// Canonical chemical-equilibrium solver using logarithmic species moles.
pub mod equilibrium_log_moles;
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
/// Typed, validated equilibrium problem construction and previews.
pub mod equilibrium_problem;
#[cfg(test)]
mod equilibrium_problem_tests;
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
/// Backend-independent validation for reconstructed equilibrium candidates.
pub mod equilibrium_validation;
#[cfg(test)]
mod equilibrium_validation_tests;
#[cfg(test)]
mod equilibrium_workflow_tests;
/// Problem construction, phase-control, and gas-equilibrium workflows.
pub mod equilibrium_workflows;
/// Typed boundary from resolved phase data to equilibrium problem construction.
pub mod phase_equilibrium_problem;
#[cfg(test)]
mod phase_equilibrium_problem_tests;
/// Immutable phase-aware result assembled from one accepted bridge solve.
pub mod phase_equilibrium_solution;
/// Narrow fixed-pressure, fixed-temperature public workflow facade.
pub mod phase_equilibrium_workflow;
