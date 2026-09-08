#![allow(deprecated)]
//! Phase-control algorithms, compatibility workflows, and fallback helpers.
//!
//! # Purpose
//!
//! This module contains the phase-control mechanics and retained mutable
//! compatibility workflows. The production entry point is
//! [`phase_equilibrium_workflow`](super::phase_equilibrium_workflow), which
//! prepares immutable active-set problems and publishes typed results. The
//! algorithms here are still used by the compatibility layer and by the
//! explicit fallback implementation, but this module is not the public
//! orchestration boundary for new code.
//!
//! The module also provides convenience functions (`gas_solver`,
//! `gas_solver_from_elements`, `gas_solver_for_T_range`) for characterization
//! and migration. They construct the historical mutable
//! [`EquilibriumLogMoles`] host and should not be used as the production
//! phase-resolved entry point.
//!
//! # Physical and Mathematical Background
//!
//! In multiphase equilibrium, the Gibbs free energy minimum may involve phases that
//! are not present in the initial guess. The phase-control algorithm:
//!
//! 1. **Solves** equilibrium with the current set of active phases.
//! 2. **Computes phase stability** — for each inactive phase, evaluates whether
//!    introducing it would lower the Gibbs free energy (using the tangent-plane
//!    criterion: `ΔG = Σ n_i(μ_i - μ_i^ref) < 0`).
//! 3. **Activates** unstable inactive phases (seeding them with trace amounts).
//! 4. **Deactivates** phases whose mole numbers fall below a destruction threshold.
//! 5. **Repeats** until no phase transitions occur or a maximum iteration count is reached.
//!
//! ## Hysteresis
//!
//! To prevent oscillatory phase transitions, the module uses **hysteresis**:
//! a phase must have `ΔG < dg_create` to be activated, but once active it is
//! kept until `ΔG > dg_keep` (with `dg_create < dg_keep`). This creates a
//! dead zone that dampens chattering.
//!
//! # Dataflow
//!
//! ```text
//!                         ┌──────────────────────────────┐
//!                         │   EquilibriumLogMoles         │
//!                         │   (inner solver, fixed set)   │
//!                         └──────────┬───────────────────┘
//!                                    │ solution + reports
//!                                    v
//!   ┌─────────────────────────────────────────────────────────┐
//!   │  PhaseManager                                            │
//!   │   ├── detect_phase_destruction()  ──> deactivation list │
//!   │   ├── classify_phases()           ──> transition plan   │
//!   │   └── thresholds_at(T)            ──> (dg_create, dg_keep)│
//!   └─────────────────────────────────────────────────────────┘
//!                                    │
//!                                    v
//!   ┌─────────────────────────────────────────────────────────┐
//!   │  PhaseSet                                                │
//!   │   ├── from_policy(policy, initial_active)               │
//!   │   ├── active_mask() ──> Vec<bool>                       │
//!   │   ├── active_phases() ──> Vec<PhaseIndex>               │
//!   │   ├── activate(phase) / deactivate(phase)               │
//!   │   └── settle_transitions()                              │
//!   └─────────────────────────────────────────────────────────┘
//!                                    │
//!                                    v
//!   ┌─────────────────────────────────────────────────────────┐
//!   │  PhaseTransitionRecord / PhaseControlledSolveReport     │
//!   │   (diagnostics: iterations, transitions, final state)   │
//!   └─────────────────────────────────────────────────────────┘
//! ```
//!
//! # Key Structures
//!
//! | Structure | Role |
//! |-----------|------|
//! | [`PhaseManager`] | Owns hysteresis policy and phase classification logic |
//! | [`PhaseSet`] | Tracks active/inactive status of all phases |
//! | [`PhaseControlledSolveReport`] | Diagnostic report of the outer loop |
//! | [`PhaseStabilityReport`] | Per-phase stability evidence |
//! | [`MultiphaseAcceptanceReport`] | Human-readable summary of the final state |
//! | [`PhaseTransitionRecord`] | One phase transition event |
//!
//! # Key Functions
//!
//! | Function | Purpose |
//! |----------|---------|
//! | [`compute_phase_stability_reports`] | Evaluates tangent-plane stability for all inactive phases |
//! | [`build_multiphase_acceptance_report`] | Assembles final acceptance report from stability data |
//! | [`seed_activated_phase_with_composition`] | Seeds a newly activated phase from its TPD minimizer |
//! | [`deactivate_phases_seed_only`] | Sets deactivated phase species to trace floor |
//! | [`compute_phase_totals`] | Sums moles by phase from log-mole data |
//! | [`initial_phase_activity`] | Computes initial activity for all phases |
//! | [`gas_solver`] | Convenience: one-step gas-phase equilibrium |
//! | [`gas_solver_for_T_range`] | Convenience: gas-phase equilibrium over T range |
//!
//! # Retained Compatibility Workflows
//!
//! The `gas_solver*` family preserves the historical gas-only setup:
//!
//! ```text
//!   SubsData ──> prepare_thermochemistry() ──> collect_gibbs_functions()
//!       │                                            │
//!       v                                            v
//!   required_element_composition()         EquilibriumLogMoles::new()
//!       │                                            │
//!       └──────────> set_problem() ──> solve() ──> result
//! ```
//!
//! # Non-obvious Details
//!
//! - The phase-control outer loop uses a **bounded iteration count** (`max_phase_iterations`)
//!   to guarantee termination even if phase transitions oscillate.
//! - `PHASE_CONTROL_TRACE_MOLE_FLOOR = 1e-300` is used when deactivating phases —
//!   this is intentionally extreme to avoid numerical underflow in the log-mole formulation.
//! - The `ActiveSetProjection` (from [`equilibrium_active_set`](super::equilibrium_active_set))
//!   constrains the solver to only consider species belonging to active phases.
//! - `finite_difference_jacobian` provides a fallback Jacobian when the symbolic
//!   engine is unavailable.
//!
//! # Related Modules
//!
//! - [`equilibrium_log_moles`](super::equilibrium_log_moles) — inner fixed-set solver
//! - [`equilibrium_active_set`](super::equilibrium_active_set) — active set projection
//! - [`equilibrium_validation`](super::equilibrium_validation) — candidate acceptance
//! - [`equilibrium_solver_policy`](super::equilibrium_solver_policy) — backend selection
//! - [`phase_equilibrium_problem`](super::phase_equilibrium_problem) — typed multiphase bridge
//!
//! # Compatibility Example
//!
//! ```rust, ignore
//! use KiThe::Thermodynamics::ChemEquilibrium::equilibrium_workflows::*;
//!
//! // Historical gas-only workflow; new code should use
//! // `ChemEquilibrium::prelude::PhaseEquilibriumPipelineRequest`.
//! let mut solver = gas_solver(
//!     vec!["CO2".to_string(), "CO".to_string(), "O2".to_string()],
//!     1500.0, 101325.0, Solvers::LM, Some("info"), true
//! );
//! solver.solve().unwrap();
//! ```
//!
use crate::Thermodynamics::ChemEquilibrium::equilibrium_active_set::ActiveSetProjection;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_activity::{
    PhaseActivityModel, phase_activity_models,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_constant_cross_validation::EquilibriumConstantCrossValidationStatus;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_constant_validation::EquilibriumConstantValidationMode;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_ids::PhaseIndex;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_log_moles::{
    EquilibriumLogMoles, EquilibriumSolveCandidate, GibbsFn, Phase, R, Solvers,
    compute_element_totals, reaction_phase_stoichiometry, species_to_phase_map,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_nonlinear::ReactionExtentError;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_phase_stability::{
    CanonicalPhaseState, IdealTpdProblem, PhaseStabilityGeometryCache,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_prepared_runner::PreparedEquilibriumRunner;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_problem::{
    EquilibriumConditions, EquilibriumProblem, LogMolesInitialGuess, PreparedEquilibriumProblem,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_solver_policy::{
    EquilibriumSolveReport, SolverPolicy,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_validation::EquilibriumCandidateReport;
use crate::Thermodynamics::User_substances::{LibraryPriority, Phases, SubsData};
use RustedSciThe::symbolic::symbolic_engine::Expr;
use log::info;
use nalgebra::{DMatrix, DVector};
use std::collections::{HashMap, HashSet};
use std::default::Default;
use std::fmt;
use std::rc::Rc;
use std::time::{Duration, Instant};

/// Lower bound used when phase-control helpers reconstruct log-mole values.
///
/// This is intentionally a named contract rather than a hidden numeric literal:
/// phase control is allowed to keep a trace amount alive, but it must do so
/// explicitly so future changes can audit the chosen floor in one place.
pub(crate) const PHASE_CONTROL_TRACE_MOLE_FLOOR: f64 = 1e-300;

/// Completes the shared thermochemistry preparation required by gas workflows.
///
/// The workflow boundary propagates database errors unchanged instead of
/// continuing with a partially populated `SubsData` instance.
fn prepare_thermochemistry(
    user_subs: &mut SubsData,
    temperature: f64,
) -> Result<(), ReactionExtentError> {
    user_subs.search_substances()?;
    user_subs.parse_all_thermal_coeffs()?;
    user_subs.extract_all_thermal_coeffs(temperature)?;
    user_subs.calculate_elem_composition_and_molar_mass(None)?;
    Ok(())
}

/// Extracts one standard-Gibbs closure per species in deterministic solver order.
fn collect_gibbs_functions(user_subs: &mut SubsData) -> Result<Vec<GibbsFn>, ReactionExtentError> {
    let mut functions = user_subs.calculate_dG0_fun_one_phase()?;
    let mut gibbs = Vec::with_capacity(user_subs.substances.len());

    for substance in &user_subs.substances {
        let function =
            functions
                .remove(substance)
                .ok_or_else(|| ReactionExtentError::InvalidProblem {
                    field: "gibbs_functions",
                    message: format!("missing standard Gibbs function for {substance}"),
                })?;
        gibbs.push(Rc::new(move |temperature: f64| function(temperature)) as GibbsFn);
    }

    Ok(gibbs)
}

/// Returns the species-by-element matrix required by the equilibrium formulation.
fn required_element_composition(user_subs: &SubsData) -> Result<DMatrix<f64>, ReactionExtentError> {
    user_subs
        .elem_composition_matrix
        .clone()
        .ok_or_else(|| ReactionExtentError::InvalidProblem {
            field: "element_composition",
            message: "thermochemistry preparation did not produce an elemental composition matrix"
                .to_string(),
        })
}

/// Total-mole policy for a TPD-derived activated phase.
///
/// The total is distributed by the accepted TPD minimizer, making the
/// physical amount and the phase composition explicit and inseparable.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum PhaseTotalSeedPolicy {
    /// Assign this positive total mole amount to the new phase.
    AbsoluteMoles {
        /// Total mole number shared among all components of the phase.
        moles: f64,
    },
    /// Use a fraction of the current total system inventory, bounded below by
    /// `minimum` total phase moles.
    RelativeToSystemTotal {
        /// Fraction of the current total system inventory.
        fraction: f64,
        /// Minimum total phase amount.
        minimum: f64,
    },
}

impl PhaseTotalSeedPolicy {
    fn total_moles(self, log_moles: &[f64]) -> Result<f64, ReactionExtentError> {
        match self {
            Self::AbsoluteMoles { moles } => {
                if !moles.is_finite() || moles <= 0.0 {
                    return Err(ReactionExtentError::InvalidProblem {
                        field: "phase_total_seed",
                        message: "absolute total phase seed must be finite and strictly positive"
                            .to_string(),
                    });
                }
                Ok(moles)
            }
            Self::RelativeToSystemTotal { fraction, minimum } => {
                if !fraction.is_finite() || fraction <= 0.0 || fraction > 1.0 {
                    return Err(ReactionExtentError::InvalidProblem {
                        field: "phase_total_seed_fraction",
                        message: "total phase seed fraction must lie in the interval (0, 1]"
                            .to_string(),
                    });
                }
                if !minimum.is_finite() || minimum <= 0.0 {
                    return Err(ReactionExtentError::InvalidProblem {
                        field: "phase_total_seed_minimum",
                        message: "total phase seed minimum must be finite and strictly positive"
                            .to_string(),
                    });
                }
                let total =
                    log_moles
                        .iter()
                        .map(|value| value.exp())
                        .try_fold(0.0, |sum, value| {
                            if !value.is_finite() {
                                Err(ReactionExtentError::InvalidProblem {
                                    field: "phase_total_seed",
                                    message: "system inventory contains non-finite mole numbers"
                                        .to_string(),
                                })
                            } else {
                                Ok(sum + value)
                            }
                        })?;
                if !total.is_finite() || total <= 0.0 {
                    return Ok(minimum);
                }
                Ok((total * fraction).max(minimum))
            }
        }
    }
}
pub use crate::Thermodynamics::ChemEquilibrium::equilibrium_phase_stability::{
    ElementPotentialReport, ElementalFeasibilityReport, PhaseStabilityConditions,
    PhaseStabilityLayout, PhaseStabilityReport, PhaseStabilityStatus, TpdMinimizerReport,
};

/// Policy used to construct the first active set of a phase-controlled solve.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum InitialPhaseSet {
    /// Derive activity from the caller's initial mole inventory.
    FromInitialMoles,
    /// Consider every declared phase as a candidate, while initially solving
    /// only phases with positive physical inventory. Inventory-free phases
    /// remain eligible for stability-driven activation after the first solve.
    ///
    /// The distinction matters for the log-moles formulation: an active phase
    /// must have a strictly positive phase total, whereas an absent phase is a
    /// valid boundary state represented by the outer active-set loop.
    AllCandidatePhases,
    /// Start with an explicit active set and permanently omit excluded phases
    /// from stability checks during this solve.
    Explicit {
        active: Vec<PhaseIndex>,
        excluded: Vec<PhaseIndex>,
    },
}

impl Default for InitialPhaseSet {
    fn default() -> Self {
        Self::FromInitialMoles
    }
}

/// Lifecycle state of one declared phase in the bounded outer loop.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum PhaseStatus {
    Active,
    Inactive,
    Excluded,
    Appeared,
    Disappeared,
}

impl PhaseStatus {
    pub(crate) fn is_active(self) -> bool {
        matches!(self, Self::Active | Self::Appeared)
    }

    fn is_candidate(self) -> bool {
        !matches!(self, Self::Excluded)
    }

    fn settled(self) -> Self {
        match self {
            Self::Appeared => Self::Active,
            Self::Disappeared => Self::Inactive,
            other => other,
        }
    }
}

/// Ordered typed state for every declared phase.
///
/// The nonlinear layer receives only [`Self::active_mask`]; all lifecycle
/// semantics remain outside residual and Jacobian callbacks.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct PhaseSet {
    statuses: Vec<PhaseStatus>,
}

impl PhaseSet {
    /// Builds the initial phase-status vector from a typed policy and the
    /// physical activity mask derived from the initial composition.
    ///
    /// The policy determines which phases start active, inactive, or excluded.
    /// At least one phase must be active after construction, or the problem
    /// has no feasible initial state.
    pub fn from_policy(
        policy: &InitialPhaseSet,
        derived_activity: &[bool],
    ) -> Result<Self, ReactionExtentError> {
        let phase_count = derived_activity.len();
        let statuses = match policy {
            InitialPhaseSet::FromInitialMoles => derived_activity
                .iter()
                .map(|&active| {
                    if active {
                        PhaseStatus::Active
                    } else {
                        PhaseStatus::Inactive
                    }
                })
                .collect(),
            InitialPhaseSet::AllCandidatePhases => vec![PhaseStatus::Active; phase_count],
            InitialPhaseSet::Explicit { active, excluded } => {
                let mut statuses = vec![PhaseStatus::Inactive; phase_count];
                for &phase in excluded {
                    if phase.index() >= phase_count {
                        return Err(ReactionExtentError::DimensionMismatch(format!(
                            "excluded phase {} is out of bounds for {phase_count} phases",
                            phase.index()
                        )));
                    }
                    statuses[phase.index()] = PhaseStatus::Excluded;
                }
                for &phase in active {
                    let index = phase.index();
                    if index >= phase_count {
                        return Err(ReactionExtentError::DimensionMismatch(format!(
                            "active phase {index} is out of bounds for {phase_count} phases"
                        )));
                    }
                    if statuses[index] == PhaseStatus::Excluded {
                        return Err(ReactionExtentError::InvalidProblem {
                            field: "initial_phase_set",
                            message: format!("phase {index} cannot be both active and excluded"),
                        });
                    }
                    statuses[index] = PhaseStatus::Active;
                }
                statuses
            }
        };
        if !statuses.iter().any(|status| status.is_active()) {
            return Err(ReactionExtentError::InvalidProblem {
                field: "initial_phase_set",
                message: "at least one phase must be active".to_string(),
            });
        }
        Ok(Self { statuses })
    }

    /// Returns the current status of one phase by typed index.
    pub fn status(&self, phase: PhaseIndex) -> PhaseStatus {
        self.statuses[phase.index()]
    }

    /// Returns a boolean mask where `true` means the phase is currently active.
    ///
    /// Active phases participate in the nonlinear solve; inactive and excluded
    /// phases are held at trace levels. The mask follows the declared phase
    /// ordering and is used by stability analysis, transition planning, and
    /// active-set projection.
    pub fn active_mask(&self) -> Vec<bool> {
        self.statuses
            .iter()
            .map(|status| status.is_active())
            .collect()
    }

    /// Returns typed indices of all currently active phases.
    ///
    /// Fails if the active-set is empty, which should never happen after
    /// [`from_policy`](Self::from_policy) has validated at least one active
    /// phase.
    pub fn active_phases(&self) -> Result<Vec<PhaseIndex>, ReactionExtentError> {
        self.statuses
            .iter()
            .enumerate()
            .filter_map(|(phase, status)| status.is_active().then_some(phase))
            .map(|phase| PhaseIndex::new(phase, self.statuses.len()))
            .collect()
    }

    pub(crate) fn is_candidate(&self, phase: usize) -> bool {
        self.statuses[phase].is_candidate()
    }

    /// Removes inventory-free phases before a positive log-moles solve.
    ///
    /// This is a policy normalization, not a thermodynamic decision. The
    /// removed phases remain candidates and can be activated later only after
    /// their accepted minimum TPD crosses the creation threshold. Returning
    /// the normalized phase indices keeps this boundary auditable for callers
    /// and tests without pretending that a solver restart already occurred.
    pub(crate) fn normalize_for_positive_solver(
        &mut self,
        physical_activity: &[bool],
    ) -> Result<Vec<PhaseIndex>, ReactionExtentError> {
        if physical_activity.len() != self.statuses.len() {
            return Err(ReactionExtentError::DimensionMismatch(format!(
                "initial phase normalization has {} activity flags for {} phases",
                physical_activity.len(),
                self.statuses.len()
            )));
        }

        let phase_count = self.statuses.len();
        let mut normalized = Vec::new();
        for (index, status) in self.statuses.iter_mut().enumerate() {
            if status.is_active() && !physical_activity[index] {
                *status = PhaseStatus::Inactive;
                normalized.push(PhaseIndex::new(index, phase_count)?);
            }
        }
        if !self.statuses.iter().any(|status| status.is_active()) {
            return Err(ReactionExtentError::InvalidProblem {
                field: "initial_phase_set",
                message:
                    "positive log-moles solve has no phase with positive physical initial inventory"
                        .to_string(),
            });
        }
        Ok(normalized)
    }

    /// Converts transient `Appeared`/`Disappeared` statuses into their settled
    /// equivalents (`Active`/`Inactive`) after a phase-control iteration.
    ///
    /// This must be called before the next transition decision so that the
    /// hysteresis policy sees a clean slate.
    pub(crate) fn settle_transitions(&mut self) {
        for status in &mut self.statuses {
            *status = status.settled();
        }
    }

    /// Marks one phase as having appeared in the current iteration.
    ///
    /// The transition is recorded as [`PhaseStatus::Appeared`]; call
    /// [`settle_transitions`](Self::settle_transitions) before the next
    /// classification pass.
    pub(crate) fn activate(&mut self, phase: PhaseIndex) {
        self.settle_transitions();
        self.statuses[phase.index()] = PhaseStatus::Appeared;
    }

    /// Marks one phase as having disappeared in the current iteration.
    ///
    /// The transition is recorded as [`PhaseStatus::Disappeared`]; call
    /// [`settle_transitions`](Self::settle_transitions) before the next
    /// classification pass.
    pub(crate) fn deactivate(&mut self, phase: PhaseIndex) {
        self.settle_transitions();
        self.statuses[phase.index()] = PhaseStatus::Disappeared;
    }
}

/// Computes typed phase-stability evidence from the canonical TPD foundation.
///
/// A numeric result is emitted only after constrained TPD minimization. The
/// outer workflow may subsequently apply hysteresis to `minimum_tpd`, but it
/// must never reinterpret a skipped evaluation as a zero-valued criterion.
#[allow(clippy::too_many_arguments)]
pub fn compute_phase_stability_reports(
    log_moles: &[f64],
    gibbs: &[GibbsFn],
    phases: &[Phase],
    species_phase: &[usize],
    element_composition: &DMatrix<f64>,
    temperature: f64,
    pressure: f64,
    p0: f64,
    phase_set: &PhaseSet,
) -> Result<Vec<PhaseStabilityReport>, ReactionExtentError> {
    let mut geometry_cache = PhaseStabilityGeometryCache::default();
    compute_phase_stability_reports_with_geometry_cache(
        log_moles,
        gibbs,
        phases,
        species_phase,
        element_composition,
        temperature,
        pressure,
        p0,
        phase_set,
        &mut geometry_cache,
    )
}

/// Cache-aware internal form of [`compute_phase_stability_reports`].
///
/// The cache is deliberately an orchestration concern: it retains only
/// immutable element-matrix geometry keyed by reference species. All current
/// chemical potentials, activities, and TPD values are reconstructed for the
/// accepted candidate state passed to this call.
#[allow(clippy::too_many_arguments)]
pub(crate) fn compute_phase_stability_reports_with_geometry_cache(
    log_moles: &[f64],
    gibbs: &[GibbsFn],
    phases: &[Phase],
    species_phase: &[usize],
    element_composition: &DMatrix<f64>,
    temperature: f64,
    pressure: f64,
    p0: f64,
    phase_set: &PhaseSet,
    geometry_cache: &mut PhaseStabilityGeometryCache,
) -> Result<Vec<PhaseStabilityReport>, ReactionExtentError> {
    let active = phase_set.active_mask();
    let species_count = log_moles.len();
    if gibbs.len() != species_count
        || species_phase.len() != species_count
        || element_composition.nrows() != species_count
        || active.len() != phases.len()
    {
        return Err(ReactionExtentError::DimensionMismatch(format!(
            "phase stability has {species_count} log-moles, {} Gibbs functions, {} phase labels, {} element rows, {} phases, and {} active flags",
            gibbs.len(),
            species_phase.len(),
            element_composition.nrows(),
            phases.len(),
            active.len()
        )));
    }
    for (parameter, value) in [
        ("temperature", temperature),
        ("pressure", pressure),
        ("reference_pressure", p0),
    ] {
        if !value.is_finite() || value <= 0.0 {
            return Err(ReactionExtentError::InvalidConditions { parameter, value });
        }
    }
    let candidate_mask = (0..phases.len())
        .map(|phase| phase_set.is_candidate(phase))
        .collect::<Vec<_>>();
    let state = CanonicalPhaseState::from_log_moles(
        log_moles,
        gibbs,
        phases,
        species_phase,
        element_composition,
        temperature,
        pressure,
        p0,
        &active,
        &candidate_mask,
    )?;
    let active = state.active();
    // Multiple declared gas phases still describe one fixed assemblage while
    // any gas is active. With a condensed-only reference, however, the gas
    // must become a real TPD candidate or complete condensation could never
    // be validated (and reverse evaporation could never be detected).
    let has_active_gas_assemblage = phases.iter().enumerate().any(|(phase_index, phase)| {
        active[phase_index] && matches!(phase.kind, PhaseActivityModel::IdealGas)
    });
    let conditions = PhaseStabilityConditions {
        temperature,
        pressure,
        reference_pressure: p0,
    };

    let mut reports = Vec::with_capacity(phases.len());
    for (phase_index, phase) in phases.iter().enumerate() {
        let phase_id = PhaseIndex::new(phase_index, phases.len())?;
        let layout = PhaseStabilityLayout {
            system_species_count: species_count,
            system_phase_count: phases.len(),
            element_count: state.element_composition().ncols(),
            phase_component_indices: phase.species.clone(),
        };
        if !state.candidates()[phase_index] {
            reports.push(PhaseStabilityReport {
                phase: phase_id,
                active: false,
                status: PhaseStabilityStatus::ExcludedByPolicy,
                conditions,
                layout,
                minimum_tpd: None,
                incipient_composition: None,
                element_potentials: None,
                elemental_feasibility: None,
                minimizer: None,
            });
            continue;
        }
        if matches!(phase.kind, PhaseActivityModel::IdealGas) && has_active_gas_assemblage {
            reports.push(PhaseStabilityReport {
                phase: phase_id,
                active: active[phase_index],
                status: PhaseStabilityStatus::FixedGasAssemblage,
                conditions,
                layout,
                minimum_tpd: None,
                incipient_composition: None,
                element_potentials: None,
                elemental_feasibility: None,
                minimizer: None,
            });
            continue;
        }
        let reference_species: Vec<usize> = (0..species_count)
            .filter(|&species| {
                let reference_phase = state.species_phase()[species];
                active[reference_phase]
                    && reference_phase != phase_index
                    && state.moles()[species] > PHASE_CONTROL_TRACE_MOLE_FLOOR
            })
            .collect();
        if reference_species.is_empty() {
            if !active[phase_index] {
                return Err(ReactionExtentError::ValidationNotApplicable {
                    path: "phase_stability",
                    message: format!(
                        "inactive candidate phase {phase_index} has no active reference assemblage for a TPD stability test"
                    ),
                });
            }
            reports.push(PhaseStabilityReport {
                phase: phase_id,
                active: active[phase_index],
                status: PhaseStabilityStatus::ActiveWithoutReferenceAssemblage,
                conditions,
                layout,
                minimum_tpd: None,
                incipient_composition: None,
                element_potentials: None,
                elemental_feasibility: None,
                minimizer: None,
            });
            continue;
        }

        let element_count = state.element_composition().ncols();
        let candidate_elements =
            DMatrix::from_fn(phase.species.len(), element_count, |row, element| {
                state.element_composition()[(phase.species[row], element)]
            });
        let reference_mu = DVector::from_iterator(
            reference_species.len(),
            reference_species
                .iter()
                .map(|&species| state.chemical_potentials()[species]),
        );
        let geometry = geometry_cache.geometry_for(&state, &reference_species)?;
        let fit = geometry.fit_potentials(&reference_mu, &reference_species)?;
        let tpd = IdealTpdProblem::new(
            phase.kind,
            phase
                .species
                .iter()
                .map(|&species| state.standard_gibbs()[species])
                .collect(),
            candidate_elements,
            fit.potentials().to_vec(),
            temperature,
            pressure,
            p0,
        )?
        .constrained_minimum_with_geometry(geometry)?;
        reports.push(PhaseStabilityReport {
            phase: phase_id,
            active: active[phase_index],
            status: PhaseStabilityStatus::Evaluated,
            conditions,
            layout,
            minimum_tpd: Some(tpd.minimum().minimum_tpd()),
            incipient_composition: Some(tpd.minimum().incipient_composition().to_vec()),
            element_potentials: Some(fit.into()),
            elemental_feasibility: Some(ElementalFeasibilityReport::from(tpd.feasibility())),
            minimizer: Some(TpdMinimizerReport::from(&tpd)),
        });
    }

    Ok(reports)
}

/// Assembles the final multiphase acceptance report from phase-control
/// history, stability evidence, and hysteresis thresholds.
///
/// The report bundles the phase-control trace, per-phase stability reports,
/// and a complementarity summary into one auditable result. Callers use it
/// to decide whether the outer phase-control loop has converged.
pub fn build_multiphase_acceptance_report(
    phase_control: PhaseControlledSolveReport,
    phase_stability: Vec<PhaseStabilityReport>,
    dg_create: Option<f64>,
    dg_keep: Option<f64>,
) -> Result<MultiphaseAcceptanceReport, ReactionExtentError> {
    if phase_control.final_phase_set.active_mask().len() != phase_stability.len() {
        return Err(ReactionExtentError::DimensionMismatch(format!(
            "acceptance bundle has {} phase-set entries and {} stability reports",
            phase_control.final_phase_set.active_mask().len(),
            phase_stability.len()
        )));
    }
    let complementarity = MultiphaseComplementarityReport::from_stability_reports(
        &phase_stability,
        dg_create,
        dg_keep,
    );
    Ok(MultiphaseAcceptanceReport {
        final_validation: phase_control.final_validation.clone(),
        phase_control,
        phase_stability,
        complementarity,
    })
}

/// Seeds a newly activated phase from its TPD minimizer composition.
///
/// This canonical path interprets the policy result as one total phase amount.
/// It distributes that total according to
/// `incipient_composition`, applies the trace floor only to produce valid
/// log-mole coordinates, and renormalizes the positive values to preserve the
/// effective requested phase total exactly.
pub fn seed_activated_phase_with_composition(
    log_moles: &mut [f64],
    phase_id: PhaseIndex,
    species_phase: &[usize],
    incipient_composition: &[f64],
    seed_policy: PhaseTotalSeedPolicy,
) -> Result<(), ReactionExtentError> {
    if log_moles.len() != species_phase.len() {
        return Err(ReactionExtentError::DimensionMismatch(format!(
            "TPD phase seed has {} log-moles and {} phase labels",
            log_moles.len(),
            species_phase.len()
        )));
    }
    let species = species_phase
        .iter()
        .enumerate()
        .filter_map(|(index, &phase)| (phase == phase_id.index()).then_some(index))
        .collect::<Vec<_>>();
    if species.is_empty() || species.len() != incipient_composition.len() {
        return Err(ReactionExtentError::DimensionMismatch(format!(
            "TPD phase seed has {} phase species and {} composition entries",
            species.len(),
            incipient_composition.len(),
        )));
    }
    if incipient_composition
        .iter()
        .any(|value| !value.is_finite() || *value < 0.0)
    {
        return Err(ReactionExtentError::InvalidCandidate {
            field: "phase_seed_composition",
            message: "TPD incipient composition must be finite and non-negative".to_string(),
        });
    }
    let composition_sum = incipient_composition.iter().sum::<f64>();
    if (composition_sum - 1.0).abs() > 1e-10 {
        return Err(ReactionExtentError::InvalidCandidate {
            field: "phase_seed_composition",
            message: format!("TPD incipient composition must sum to one, got {composition_sum}"),
        });
    }

    let requested_total = seed_policy.total_moles(log_moles)?;
    let effective_total =
        requested_total.max(species.len() as f64 * PHASE_CONTROL_TRACE_MOLE_FLOOR);
    let floor_fraction = PHASE_CONTROL_TRACE_MOLE_FLOOR / effective_total;
    let mut weights = incipient_composition
        .iter()
        .map(|&value| value.max(floor_fraction))
        .collect::<Vec<_>>();
    let weight_sum = weights.iter().sum::<f64>();
    if !weight_sum.is_finite() || weight_sum <= 0.0 {
        return Err(ReactionExtentError::InvalidCandidate {
            field: "phase_seed_composition",
            message: "TPD seed composition normalization is invalid".to_string(),
        });
    }
    for weight in &mut weights {
        *weight /= weight_sum;
    }
    for (species_index, weight) in species.into_iter().zip(weights) {
        let moles = effective_total * weight;
        if !moles.is_finite() || moles < PHASE_CONTROL_TRACE_MOLE_FLOOR {
            return Err(ReactionExtentError::InvalidCandidate {
                field: "phase_seed_moles",
                message: "TPD phase seed produced an invalid positive mole amount".to_string(),
            });
        }
        log_moles[species_index] = moles.ln();
    }
    Ok(())
}

/// Sets log-mole coordinates for deactivated phases to a trace floor.
///
/// This is a seed-only operation: it does not modify the phase-set status.
/// The caller is responsible for updating the phase-set separately. The
/// trace floor must be finite and strictly positive.
pub fn deactivate_phases_seed_only(
    y: &mut [f64], // log-moles (modified in place)
    phases_to_deactivate: &[PhaseIndex],
    species_phase: &[usize],
    trace_floor: f64,
) -> Result<(), ReactionExtentError> {
    if !trace_floor.is_finite() || trace_floor <= 0.0 {
        return Err(ReactionExtentError::InvalidProblem {
            field: "trace_floor",
            message: "phase-control trace floor must be finite and strictly positive".to_string(),
        });
    }
    let floor = trace_floor.ln();
    for &phase_id in phases_to_deactivate {
        for (index, &phase) in species_phase.iter().enumerate() {
            if phase == phase_id.index() {
                y[index] = floor;
            }
        }
    }
    Ok(())
}

/// Why a phase-control pass ended without a transition.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PhaseStabilityReason {
    /// No active phase fell below the destruction threshold and no inactive
    /// phase qualified for creation.
    NoPhaseTransitionNeeded,
}

/// Typed result of one phase-control transition decision.
///
/// The phase-control loop is intentionally bounded to one explicit action per
/// restart. That keeps the restart logic auditable and prevents a single outer
/// iteration from trying to reshuffle several phases at once.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PhaseTransitionPlan {
    /// Deactivate one phase that has drifted below the destruction threshold.
    Deactivate { phase: PhaseIndex },
    /// Activate one missing phase whose creation score is sufficiently low.
    Activate { phase: PhaseIndex },
    /// Keep one marginally small active phase because the hysteresis band says
    /// it is still thermodynamically justified.
    Hold { phase: PhaseIndex },
    /// No phase transition is needed for this restart.
    NoTransition { reason: PhaseStabilityReason },
}

/// Physical reason recorded for one active-set transition.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum PhaseTransitionReason {
    UnstableInactivePhase {
        minimum_tpd: f64,
    },
    VanishingUnstableActivePhase {
        phase_moles: f64,
        minimum_tpd: f64,
    },
    /// The positive interior problem failed at a phase boundary, and a
    /// validated solve without that phase proved it thermodynamically absent.
    BoundaryUnstableActivePhase {
        initial_phase_moles: f64,
        minimum_tpd: f64,
    },
}

/// One auditable active-set change between completed nonlinear solves.
#[derive(Debug, Clone, PartialEq)]
pub struct PhaseTransitionRecord {
    pub iteration: usize,
    /// Wall time spent in the phase-control pass that produced this change.
    /// This includes the fixed-set candidate and the stability/transition
    /// decision, and is zero only for synthetic records.
    pub transition_duration: Duration,
    pub activated: Vec<PhaseIndex>,
    pub deactivated: Vec<PhaseIndex>,
    pub phase_totals: Vec<f64>,
    /// Per-phase constrained TPD minima in phase-index order. `None` means
    /// the phase was not a TPD candidate under the declared workflow policy.
    pub minimum_tpds: Vec<Option<f64>>,
    /// TPD minimizer used to seed an activated phase. Deactivation records do
    /// not carry one. Retaining this independently of `restart_seed` makes
    /// the physical transition auditable without reverse-engineering floored
    /// log coordinates.
    pub incipient_composition: Option<Vec<f64>>,
    pub reason: PhaseTransitionReason,
    pub previous_phase_set: PhaseSet,
    pub new_phase_set: PhaseSet,
    /// Log-mole seed passed to the next fixed-active-set solve.
    pub restart_seed: Vec<f64>,
    /// Backend cascade that produced the accepted fixed-set candidate.
    pub nonlinear_report: EquilibriumSolveReport,
    /// Acceptance evidence for the fixed-set solution that triggered this
    /// transition. Seed vectors themselves are never treated as solutions.
    pub candidate_validation: EquilibriumCandidateReport,
}

/// Evidence produced by one bounded phase-control solve.
#[derive(Debug, Clone, PartialEq)]
pub struct PhaseControlledSolveReport {
    pub iterations: usize,
    pub initial_active_phases: Vec<PhaseIndex>,
    pub final_active_phases: Vec<PhaseIndex>,
    pub initial_phase_set: PhaseSet,
    pub final_phase_set: PhaseSet,
    pub transitions: Vec<PhaseTransitionRecord>,
    pub final_validation: EquilibriumCandidateReport,
    pub nonlinear_reports: Vec<EquilibriumSolveReport>,
}

/// Numerical complementarity evidence for the final multiphase state.
#[derive(Debug, Clone, PartialEq)]
pub struct MultiphaseComplementarityReport {
    pub supported_phase_count: usize,
    pub active_supported_phase_count: usize,
    pub inactive_supported_phase_count: usize,
    pub max_active_minimum_tpd: Option<f64>,
    pub min_inactive_minimum_tpd: Option<f64>,
    pub max_inactive_minimum_tpd: Option<f64>,
    pub dg_create: Option<f64>,
    pub dg_keep: Option<f64>,
    pub satisfied: bool,
}

/// Final acceptance bundle for one multiphase equilibrium solve.
#[derive(Debug, Clone, PartialEq)]
pub struct MultiphaseAcceptanceReport {
    pub phase_control: PhaseControlledSolveReport,
    pub phase_stability: Vec<PhaseStabilityReport>,
    pub complementarity: MultiphaseComplementarityReport,
    pub final_validation: EquilibriumCandidateReport,
}

/// One stable summary row for the final multiphase acceptance bundle.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MultiphaseAcceptanceRow {
    /// Logical section name, e.g. `acceptance`, `validation`, or `complementarity`.
    pub section: &'static str,
    /// Stable row label.
    pub label: String,
    /// Human-readable value.
    pub value: String,
}

impl fmt::Display for MultiphaseAcceptanceRow {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "[{}] {} = {}", self.section, self.label, self.value)
    }
}

/// One stable summary row for a phase-control solve report.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PhaseControlledSolveRow {
    /// Logical section name, e.g. `phase_control`, `phase_set`, or `validation`.
    pub section: &'static str,
    /// Stable row label.
    pub label: String,
    /// Human-readable value.
    pub value: String,
}

impl fmt::Display for PhaseControlledSolveRow {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "[{}] {} = {}", self.section, self.label, self.value)
    }
}

impl PhaseControlledSolveReport {
    /// Returns stable summary rows for CLI output and snapshot tests.
    pub fn summary_rows(&self) -> Vec<PhaseControlledSolveRow> {
        let mut rows = vec![
            PhaseControlledSolveRow {
                section: "phase_control",
                label: "iterations".to_string(),
                value: self.iterations.to_string(),
            },
            PhaseControlledSolveRow {
                section: "phase_control",
                label: "transitions".to_string(),
                value: self.transitions.len().to_string(),
            },
            PhaseControlledSolveRow {
                section: "phase_control",
                label: "initial_active_phases".to_string(),
                value: self
                    .initial_active_phases
                    .iter()
                    .map(|phase| phase.index().to_string())
                    .collect::<Vec<_>>()
                    .join(", "),
            },
            PhaseControlledSolveRow {
                section: "phase_control",
                label: "final_active_phases".to_string(),
                value: self
                    .final_active_phases
                    .iter()
                    .map(|phase| phase.index().to_string())
                    .collect::<Vec<_>>()
                    .join(", "),
            },
            PhaseControlledSolveRow {
                section: "validation",
                label: "residual_l2_norm".to_string(),
                value: format!("{:.6e}", self.final_validation.residual_l2_norm),
            },
            PhaseControlledSolveRow {
                section: "validation",
                label: "max_abs_element_balance_error".to_string(),
                value: format!(
                    "{:.6e}",
                    self.final_validation.max_abs_element_balance_error
                ),
            },
            PhaseControlledSolveRow {
                section: "validation",
                label: "min_moles".to_string(),
                value: format!("{:.6e}", self.final_validation.min_moles),
            },
            PhaseControlledSolveRow {
                section: "nonlinear",
                label: "accepted_reports".to_string(),
                value: self.nonlinear_reports.len().to_string(),
            },
        ];

        if let Some(last_transition) = self.transitions.last() {
            rows.push(PhaseControlledSolveRow {
                section: "phase_control",
                label: "last_transition_reason".to_string(),
                value: format!("{:?}", last_transition.reason),
            });
        }

        rows
    }
}

impl fmt::Display for PhaseControlledSolveReport {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for row in self.summary_rows() {
            writeln!(f, "{row}")?;
        }
        Ok(())
    }
}

impl MultiphaseComplementarityReport {
    fn from_stability_reports(
        stability: &[PhaseStabilityReport],
        dg_create: Option<f64>,
        dg_keep: Option<f64>,
    ) -> Self {
        let mut supported_phase_count = 0;
        let mut active_supported_phase_count = 0;
        let mut inactive_supported_phase_count = 0;
        let mut max_active_minimum_tpd = None;
        let mut min_inactive_minimum_tpd = None;
        let mut max_inactive_minimum_tpd = None;
        let mut satisfied = true;

        for report in stability {
            if let Some(force) = report.minimum_tpd {
                supported_phase_count += 1;
                if report.active {
                    active_supported_phase_count += 1;
                    max_active_minimum_tpd =
                        Some(max_active_minimum_tpd.map_or(force, |value: f64| value.max(force)));
                    if let Some(dg_keep) = dg_keep {
                        if force > dg_keep {
                            satisfied = false;
                        }
                    }
                } else {
                    inactive_supported_phase_count += 1;
                    min_inactive_minimum_tpd =
                        Some(min_inactive_minimum_tpd.map_or(force, |value: f64| value.min(force)));
                    max_inactive_minimum_tpd =
                        Some(max_inactive_minimum_tpd.map_or(force, |value: f64| value.max(force)));
                    if let Some(dg_create) = dg_create {
                        if force < dg_create {
                            satisfied = false;
                        }
                    }
                }
            }
        }

        Self {
            supported_phase_count,
            active_supported_phase_count,
            inactive_supported_phase_count,
            max_active_minimum_tpd,
            min_inactive_minimum_tpd,
            max_inactive_minimum_tpd,
            dg_create,
            dg_keep,
            satisfied,
        }
    }
}

impl MultiphaseAcceptanceReport {
    /// Returns stable summary rows for CLI output and snapshot tests.
    pub fn summary_rows(&self) -> Vec<MultiphaseAcceptanceRow> {
        let mut rows = vec![
            MultiphaseAcceptanceRow {
                section: "acceptance",
                label: "phase_iterations".to_string(),
                value: self.phase_control.iterations.to_string(),
            },
            MultiphaseAcceptanceRow {
                section: "acceptance",
                label: "phase_transitions".to_string(),
                value: self.phase_control.transitions.len().to_string(),
            },
            MultiphaseAcceptanceRow {
                section: "acceptance",
                label: "supported_phases".to_string(),
                value: self.complementarity.supported_phase_count.to_string(),
            },
            MultiphaseAcceptanceRow {
                section: "acceptance",
                label: "active_supported_phases".to_string(),
                value: self
                    .complementarity
                    .active_supported_phase_count
                    .to_string(),
            },
            MultiphaseAcceptanceRow {
                section: "acceptance",
                label: "inactive_supported_phases".to_string(),
                value: self
                    .complementarity
                    .inactive_supported_phase_count
                    .to_string(),
            },
            MultiphaseAcceptanceRow {
                section: "validation",
                label: "residual_l2_norm".to_string(),
                value: format!("{:.6e}", self.final_validation.residual_l2_norm),
            },
            MultiphaseAcceptanceRow {
                section: "validation",
                label: "max_abs_element_balance_error".to_string(),
                value: format!(
                    "{:.6e}",
                    self.final_validation.max_abs_element_balance_error
                ),
            },
            MultiphaseAcceptanceRow {
                section: "validation",
                label: "min_moles".to_string(),
                value: format!("{:.6e}", self.final_validation.min_moles),
            },
            MultiphaseAcceptanceRow {
                section: "complementarity",
                label: "satisfied".to_string(),
                value: self.complementarity.satisfied.to_string(),
            },
        ];

        if let Some(value) = self.complementarity.max_active_minimum_tpd {
            rows.push(MultiphaseAcceptanceRow {
                section: "complementarity",
                label: "max_active_minimum_tpd_j_per_mol".to_string(),
                value: format!("{:.6e}", value),
            });
        }
        if let Some(value) = self.complementarity.min_inactive_minimum_tpd {
            rows.push(MultiphaseAcceptanceRow {
                section: "complementarity",
                label: "min_inactive_minimum_tpd_j_per_mol".to_string(),
                value: format!("{:.6e}", value),
            });
        }
        if let Some(value) = self.complementarity.max_inactive_minimum_tpd {
            rows.push(MultiphaseAcceptanceRow {
                section: "complementarity",
                label: "max_inactive_minimum_tpd_j_per_mol".to_string(),
                value: format!("{:.6e}", value),
            });
        }

        // Keep per-phase TPD evidence separate from aggregate
        // complementarity. A small aggregate can otherwise hide a rejected
        // minimizer or a poor elemental/KKT residual in one candidate phase.
        for stability in &self.phase_stability {
            let phase = stability.phase.index();
            rows.push(MultiphaseAcceptanceRow {
                section: "phase_stability",
                label: format!("phase_{phase}_status"),
                value: format!("{:?}", stability.status),
            });
            if let Some(value) = stability.minimum_tpd {
                rows.push(MultiphaseAcceptanceRow {
                    section: "phase_stability",
                    label: format!("phase_{phase}_minimum_tpd_j_per_mol"),
                    value: format!("{value:.6e}"),
                });
            }
            if let Some(feasibility) = &stability.elemental_feasibility {
                rows.push(MultiphaseAcceptanceRow {
                    section: "phase_stability",
                    label: format!("phase_{phase}_elemental_feasibility_residual"),
                    value: format!(
                        "{:.6e} <= {:.6e}",
                        feasibility.max_abs_residual, feasibility.residual_tolerance
                    ),
                });
            }
            if let Some(minimizer) = &stability.minimizer {
                rows.push(MultiphaseAcceptanceRow {
                    section: "phase_stability",
                    label: format!("phase_{phase}_kkt_residual"),
                    value: format!("{:.6e}", minimizer.max_abs_kkt_residual),
                });
            }
        }
        rows
    }
}

impl fmt::Display for MultiphaseAcceptanceReport {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for row in self.summary_rows() {
            writeln!(f, "{row}")?;
        }
        Ok(())
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum PhaseHysteresisPolicy {
    TemperatureScaled {
        create_rt_factor: f64,
        keep_rt_factor: f64,
    },
    Explicit {
        dg_create: f64,
        dg_keep: f64,
    },
}

impl Default for PhaseHysteresisPolicy {
    fn default() -> Self {
        Self::TemperatureScaled {
            create_rt_factor: -1e-6,
            keep_rt_factor: 1e-8,
        }
    }
}

/// Manages phase activation/deactivation decisions during multiphase equilibrium.
///
/// The phase manager owns hysteresis thresholds, the destruction epsilon, and
/// the policy for constructing the initial active set. It is deliberately
/// separate from the nonlinear solver state.
#[derive(Debug, Clone)]
pub struct PhaseManager {
    /// Destruction threshold: phases with total moles below this value are candidates for deactivation.
    pub phase_eps: f64,
    /// Hysteresis policy controlling the gap between phase creation and retention thresholds.
    pub phase_hysteresis: PhaseHysteresisPolicy,
    /// Maximum number of outer phase-control iterations before giving up.
    pub max_phase_iterations: usize,
    /// Policy for constructing the initial active set from the declared phases.
    pub initial_phase_set: InitialPhaseSet,
}

impl PhaseManager {
    pub fn new(phase_eps: f64, dg_create: f64, dg_keep: f64) -> Self {
        Self {
            phase_eps,
            phase_hysteresis: PhaseHysteresisPolicy::Explicit { dg_create, dg_keep },
            max_phase_iterations: 16,
            initial_phase_set: InitialPhaseSet::default(),
        }
    }

    pub fn with_temperature_scaled_hysteresis(
        phase_eps: f64,
        create_rt_factor: f64,
        keep_rt_factor: f64,
    ) -> Self {
        Self {
            phase_eps,
            phase_hysteresis: PhaseHysteresisPolicy::TemperatureScaled {
                create_rt_factor,
                keep_rt_factor,
            },
            max_phase_iterations: 16,
            initial_phase_set: InitialPhaseSet::default(),
        }
    }

    /// Validates policy fields that do not depend on the number of declared
    /// phases. Phase-index bounds are checked later when the resolved layout
    /// is available; malformed scalar policy values must fail at the public
    /// facade boundary instead of inside an outer-loop iteration.
    pub fn validate(&self) -> Result<(), ReactionExtentError> {
        if !self.phase_eps.is_finite() || self.phase_eps <= 0.0 {
            return Err(ReactionExtentError::InvalidProblem {
                field: "phase_eps",
                message: "phase_eps must be finite and strictly positive".to_string(),
            });
        }
        if self.max_phase_iterations == 0 {
            return Err(ReactionExtentError::InvalidProblem {
                field: "max_phase_iterations",
                message: "max_phase_iterations must be greater than zero".to_string(),
            });
        }
        self.thresholds_at(298.15)?;

        if let InitialPhaseSet::Explicit { active, excluded } = &self.initial_phase_set {
            if active.iter().any(|phase| excluded.contains(phase)) {
                return Err(ReactionExtentError::InvalidProblem {
                    field: "initial_phase_set",
                    message: "a phase cannot be both active and excluded".to_string(),
                });
            }
            for (position, phase) in active.iter().enumerate() {
                if active.iter().skip(position + 1).any(|other| other == phase) {
                    return Err(ReactionExtentError::InvalidProblem {
                        field: "initial_phase_set.active",
                        message: format!("phase {} is listed more than once", phase.index()),
                    });
                }
            }
            for (position, phase) in excluded.iter().enumerate() {
                if excluded
                    .iter()
                    .skip(position + 1)
                    .any(|other| other == phase)
                {
                    return Err(ReactionExtentError::InvalidProblem {
                        field: "initial_phase_set.excluded",
                        message: format!("phase {} is listed more than once", phase.index()),
                    });
                }
            }
        }
        Ok(())
    }

    /// Validates this policy against one resolved phase count.
    ///
    /// Scalar policy validation cannot check typed phase indices because the
    /// resolved layout is not known when `PhaseControlPolicy` is constructed.
    /// The workflow calls this immediately before active-set construction so
    /// an out-of-range policy fails at the boundary, never during a later
    /// transition.
    pub fn validate_for_phase_count(&self, phase_count: usize) -> Result<(), ReactionExtentError> {
        self.validate()?;
        if phase_count == 0 {
            return Err(ReactionExtentError::InvalidProblem {
                field: "phase_count",
                message: "at least one declared phase is required".to_string(),
            });
        }
        if let InitialPhaseSet::Explicit { active, excluded } = &self.initial_phase_set {
            for (field, phases) in [
                ("initial_phase_set.active", active),
                ("initial_phase_set.excluded", excluded),
            ] {
                if let Some(phase) = phases.iter().find(|phase| phase.index() >= phase_count) {
                    return Err(ReactionExtentError::DimensionMismatch(format!(
                        "{field} contains phase {} outside resolved phase count {phase_count}",
                        phase.index()
                    )));
                }
            }
        }
        Ok(())
    }

    pub fn set_explicit_hysteresis(&mut self, dg_create: f64, dg_keep: f64) {
        self.phase_hysteresis = PhaseHysteresisPolicy::Explicit { dg_create, dg_keep };
    }

    pub fn set_temperature_scaled_hysteresis(
        &mut self,
        create_rt_factor: f64,
        keep_rt_factor: f64,
    ) {
        self.phase_hysteresis = PhaseHysteresisPolicy::TemperatureScaled {
            create_rt_factor,
            keep_rt_factor,
        };
    }

    pub fn thresholds_at(&self, temperature: f64) -> Result<(f64, f64), ReactionExtentError> {
        match self.phase_hysteresis {
            PhaseHysteresisPolicy::TemperatureScaled {
                create_rt_factor,
                keep_rt_factor,
            } => {
                if !temperature.is_finite() || temperature <= 0.0 {
                    return Err(ReactionExtentError::InvalidProblem {
                        field: "temperature",
                        message: format!(
                            "phase hysteresis requires a positive finite temperature, got {temperature}"
                        ),
                    });
                }
                let dg_create = create_rt_factor * R * temperature;
                let dg_keep = keep_rt_factor * R * temperature;
                if !dg_create.is_finite() || !dg_keep.is_finite() {
                    return Err(ReactionExtentError::InvalidProblem {
                        field: "phase_hysteresis",
                        message: "phase hysteresis thresholds must be finite".to_string(),
                    });
                }
                if dg_create >= dg_keep {
                    return Err(ReactionExtentError::InvalidProblem {
                        field: "phase_hysteresis",
                        message: format!(
                            "phase creation threshold {dg_create} must be strictly below keep threshold {dg_keep}"
                        ),
                    });
                }
                Ok((dg_create, dg_keep))
            }
            PhaseHysteresisPolicy::Explicit { dg_create, dg_keep } => {
                if !dg_create.is_finite() || !dg_keep.is_finite() {
                    return Err(ReactionExtentError::InvalidProblem {
                        field: "phase_hysteresis",
                        message: "phase hysteresis thresholds must be finite".to_string(),
                    });
                }
                if dg_create >= dg_keep {
                    return Err(ReactionExtentError::InvalidProblem {
                        field: "phase_hysteresis",
                        message: format!(
                            "phase creation threshold {dg_create} must be strictly below keep threshold {dg_keep}"
                        ),
                    });
                }
                Ok((dg_create, dg_keep))
            }
        }
    }
    /// Detect candidate phases whose accepted minimum TPD requests creation.
    /// Detect phases that must be removed (n_phase < eps)
    pub fn detect_phase_destruction(&self, n_phase: &[f64], active: &[bool]) -> Vec<usize> {
        n_phase
            .iter()
            .enumerate()
            .filter(|&(p, &n)| active[p] && n < self.phase_eps)
            .map(|(p, _)| p)
            .collect()
    }

    /// Classifies phase transitions for one outer restart.
    ///
    /// The current heuristic keeps the loop bounded by selecting at most one
    /// phase transition per pass. Deactivation has priority when a phase is
    /// already below the destruction threshold.
    pub fn classify_phases(
        &self,
        phase_totals: &[f64],
        stability: &[PhaseStabilityReport],
        phase_set: &PhaseSet,
    ) -> Result<PhaseTransitionPlan, ReactionExtentError> {
        match self.phase_hysteresis {
            PhaseHysteresisPolicy::Explicit { dg_create, dg_keep } => self
                .classify_phases_with_thresholds(phase_totals, stability, phase_set, dg_create, dg_keep),
            PhaseHysteresisPolicy::TemperatureScaled { .. } => Err(
                ReactionExtentError::InvalidProblem {
                    field: "phase_hysteresis",
                    message: "temperature-scaled hysteresis must be resolved with classify_phases_at_temperature".to_string(),
                },
            ),
        }
    }

    fn classify_phases_with_thresholds(
        &self,
        phase_totals: &[f64],
        stability: &[PhaseStabilityReport],
        phase_set: &PhaseSet,
        dg_create: f64,
        dg_keep: f64,
    ) -> Result<PhaseTransitionPlan, ReactionExtentError> {
        let active = phase_set.active_mask();
        if phase_totals.len() != active.len() || stability.len() != active.len() {
            return Err(ReactionExtentError::DimensionMismatch(format!(
                "phase classification has {} totals, {} stability reports, and {} active flags",
                phase_totals.len(),
                stability.len(),
                active.len()
            )));
        }
        for (index, report) in stability.iter().enumerate() {
            if report.phase.index() != index || report.active != active[index] {
                return Err(ReactionExtentError::InvalidProblem {
                    field: "phase_stability",
                    message: format!(
                        "stability report {index} does not match the active-set ordering"
                    ),
                });
            }
        }

        let active_count = active.iter().filter(|&&is_active| is_active).count();
        if active_count > 1 {
            if let Some((phase, _)) = phase_totals
                .iter()
                .enumerate()
                .filter(|&(phase, &total)| {
                    active[phase]
                        && total < self.phase_eps
                        && stability[phase]
                            .minimum_tpd
                            .is_some_and(|force| force > dg_keep)
                })
                .max_by(|(phase_a, _), (phase_b, _)| {
                    stability[*phase_a]
                        .minimum_tpd
                        .partial_cmp(&stability[*phase_b].minimum_tpd)
                        .unwrap_or(std::cmp::Ordering::Equal)
                })
            {
                return Ok(PhaseTransitionPlan::Deactivate {
                    phase: stability[phase].phase,
                });
            }
        }

        if let Some((phase, _)) = stability
            .iter()
            .enumerate()
            .filter(|&(phase, report)| {
                !active[phase]
                    && phase_set.is_candidate(phase)
                    && report.minimum_tpd.is_some_and(|force| force < dg_create)
            })
            .min_by(|(_, a), (_, b)| {
                a.minimum_tpd
                    .partial_cmp(&b.minimum_tpd)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
        {
            return Ok(PhaseTransitionPlan::Activate {
                phase: stability[phase].phase,
            });
        }

        if let Some((phase, _)) = phase_totals
            .iter()
            .enumerate()
            .find(|&(phase, &total)| active[phase] && total < self.phase_eps)
        {
            return Ok(PhaseTransitionPlan::Hold {
                phase: stability[phase].phase,
            });
        }

        Ok(PhaseTransitionPlan::NoTransition {
            reason: PhaseStabilityReason::NoPhaseTransitionNeeded,
        })
    }

    /// Classifies phase transitions using temperature-scaled hysteresis thresholds.
    ///
    /// This is a convenience wrapper that resolves the temperature-dependent
    /// thresholds via [`thresholds_at`](Self::thresholds_at) and delegates to
    /// the private [`classify_phases_with_thresholds`] method.
    pub fn classify_phases_at_temperature(
        &self,
        temperature: f64,
        phase_totals: &[f64],
        stability: &[PhaseStabilityReport],
        phase_set: &PhaseSet,
    ) -> Result<PhaseTransitionPlan, ReactionExtentError> {
        let (dg_create, dg_keep) = self.thresholds_at(temperature)?;
        self.classify_phases_with_thresholds(phase_totals, stability, phase_set, dg_create, dg_keep)
    }
}
impl Default for PhaseManager {
    fn default() -> Self {
        Self {
            phase_eps: 1e-30,
            phase_hysteresis: PhaseHysteresisPolicy::default(),
            max_phase_iterations: 16,
            initial_phase_set: InitialPhaseSet::default(),
        }
    }
}
/// Computes total physical moles in each phase from log-mole coordinates.
///
/// The result has one entry per phase, ordered by phase index. This is a
/// pure function that does not inspect solver state or modify any cached
/// problem data.
pub fn compute_phase_totals(
    y: &[f64],               // log-moles
    species_phase: &[usize], // m
) -> Vec<f64> {
    let m = y.len();
    let n_phases = species_phase.iter().copied().max().unwrap_or(0) + 1;

    let mut n_phase = vec![0.0; n_phases];

    for i in 0..m {
        let p = species_phase[i];
        n_phase[p] += y[i].exp();
    }

    n_phase
}

/// Derives a phase activity mask from physical component mole numbers.
///
/// This is the physical boundary for the first outer-loop iteration. Zero
/// inventory means that a candidate phase starts absent; it must not become
/// active merely because the nonlinear solver later represents the same zero
/// with a positive trace seed.
pub fn initial_phase_activity_from_moles(
    moles: &[f64],
    species_phase: &[usize],
    phase_count: usize,
    phase_eps: f64,
) -> Result<Vec<bool>, ReactionExtentError> {
    if moles.len() != species_phase.len() {
        return Err(ReactionExtentError::DimensionMismatch(format!(
            "phase activity has {} mole values and {} phase labels",
            moles.len(),
            species_phase.len()
        )));
    }
    if !phase_eps.is_finite() || phase_eps < 0.0 {
        return Err(ReactionExtentError::InvalidProblem {
            field: "phase_eps",
            message: "phase activity threshold must be finite and non-negative".to_string(),
        });
    }

    let mut totals = vec![0.0; phase_count];
    for (species, (&moles, &phase)) in moles.iter().zip(species_phase).enumerate() {
        if !moles.is_finite() || moles < 0.0 {
            return Err(ReactionExtentError::InvalidCandidate {
                field: "initial_phase_activity",
                message: format!("species {species} has invalid physical initial moles {moles:e}"),
            });
        }
        let total = totals.get_mut(phase).ok_or_else(|| {
            ReactionExtentError::DimensionMismatch(format!(
                "species {species} refers to missing phase {phase}"
            ))
        })?;
        *total += moles;
    }
    Ok(totals.into_iter().map(|total| total > phase_eps).collect())
}

/// Derives a phase activity mask from a numerical log-mole seed.
///
/// This helper remains useful for explicit numerical warm starts. First-pass
/// phase-control initialization instead uses
/// [`initial_phase_activity_from_moles`], preserving the physical distinction
/// between an absent phase and a trace-coordinate representation of zero.
pub fn initial_phase_activity(
    log_moles: &[f64],
    species_phase: &[usize],
    phase_count: usize,
    phase_eps: f64,
) -> Result<Vec<bool>, ReactionExtentError> {
    let mut physical_moles = Vec::with_capacity(log_moles.len());
    for (species, &log_moles) in log_moles.iter().enumerate() {
        let moles = log_moles.exp();
        if !moles.is_finite() {
            return Err(ReactionExtentError::InvalidCandidate {
                field: "initial_phase_activity",
                message: format!("species {species} reconstructs non-finite moles"),
            });
        }
        physical_moles.push(moles);
    }
    initial_phase_activity_from_moles(&physical_moles, species_phase, phase_count, phase_eps)
}

/// Detects phase-control cycles by comparing the current phase set against
/// all previously visited sets.
///
/// A cycle is declared when the settled phase set has been seen before. This
/// prevents infinite loops in the outer phase-control iteration. The
/// iteration number is preserved in the error for diagnostics.
pub(crate) fn reject_repeated_phase_set(
    visited: &mut HashSet<PhaseSet>,
    phase_set: &PhaseSet,
    iteration: usize,
) -> Result<(), ReactionExtentError> {
    let mut settled = phase_set.clone();
    settled.settle_transitions();
    if visited.insert(settled.clone()) {
        return Ok(());
    }
    Err(ReactionExtentError::PhaseControlCycleDetected {
        iteration,
        active_phases: settled
            .active_phases()?
            .into_iter()
            .map(PhaseIndex::index)
            .collect(),
    })
}

/// Validates that a phase-set candidate has consistent dimensions and that
/// every active phase has positive physical inventory.
///
/// This is called after each nonlinear solve to ensure the candidate phase
/// configuration is physically meaningful before accepting a transition.
pub(crate) fn validate_phase_set_candidate(
    phase_totals: &[f64],
    active: &[bool],
    phase_eps: f64,
) -> Result<(), ReactionExtentError> {
    if phase_totals.len() != active.len() {
        return Err(ReactionExtentError::DimensionMismatch(format!(
            "phase candidate has {} totals and {} active flags",
            phase_totals.len(),
            active.len()
        )));
    }
    for (phase, (&total, &is_active)) in phase_totals.iter().zip(active).enumerate() {
        if !total.is_finite() || total < 0.0 {
            return Err(ReactionExtentError::InvalidCandidate {
                field: "phase_totals",
                message: format!("phase {phase} has invalid total {total}"),
            });
        }
        if !is_active && total > phase_eps {
            return Err(ReactionExtentError::InvalidCandidate {
                field: "phase_active_set",
                message: format!(
                    "inactive phase {phase} contains {total} mol, above the {phase_eps} mol threshold"
                ),
            });
        }
    }
    Ok(())
}

/*

///
///  loop:
///     solve equilibrium (NR / trust region)
///
///     compute n_phase
///
///     if inactive phases detected:
///         deactivate phases
///         restart solver
///
///     if ΔG < 0 for missing phase:
///         activate phase
///         restart solver
///
///     break
///
///



*/
impl EquilibriumLogMoles {
    /// Solves exactly one fixed active set and expands it back to the declared
    /// species ordering without publishing intermediate state.
    fn solve_fixed_active_set_candidate(
        &mut self,
        active: &[bool],
        full_seed: &[f64],
    ) -> Result<EquilibriumSolveCandidate, ReactionExtentError> {
        if active.len() != self.phases.len() || full_seed.len() != self.n0.len() {
            return Err(ReactionExtentError::DimensionMismatch(format!(
                "fixed active set has {} flags and {} seed values for {} phases and {} species",
                active.len(),
                full_seed.len(),
                self.phases.len(),
                self.n0.len()
            )));
        }
        let projection = ActiveSetProjection::build(
            &self.phases,
            &self.species_phase,
            &self.elem_composition,
            active,
            self.solver_settings.solver_params.tol,
        )?;
        // These totals are the closed-system physical invariant. They are
        // deliberately recomputed from the full composition, never inferred
        // from the reduced log-space seed or any trace-floor values.
        let full_element_totals = compute_element_totals(&self.elem_composition, &self.n0)?
            .iter()
            .copied()
            .collect::<Vec<_>>();
        projection.validate_element_totals_representable(
            &full_element_totals,
            self.solver_settings.solver_params.tol,
        )?;
        let reduced_seed = projection.project_log_moles(full_seed)?;

        let reduced_components = projection
            .active_species
            .iter()
            .map(|species| self.subs_data.substances[species.index()].clone())
            .collect();
        let reduced_initial_moles = reduced_seed
            .iter()
            .map(|value| value.exp())
            .collect::<Vec<_>>();
        let reduced_gibbs = projection
            .active_species
            .iter()
            .map(|species| self.gibbs[species.index()].clone())
            .collect();
        let reduced_gibbs_sym = if self.gibbs_sym.len() == self.n0.len() {
            projection
                .active_species
                .iter()
                .map(|species| self.gibbs_sym[species.index()].clone())
                .collect()
        } else {
            Vec::new()
        };
        let conditions = EquilibriumConditions::new(self.T, self.P, self.p0)?;
        let reduced_problem = EquilibriumProblem::new(
            reduced_components,
            reduced_initial_moles,
            LogMolesInitialGuess::new(reduced_seed)?,
            projection.element_composition.clone(),
            reduced_gibbs,
            projection.phases.clone(),
            conditions,
        )?;
        let reduced_totals = projection.reduced_element_totals(&full_element_totals)?;
        let prepared = PreparedEquilibriumProblem::new_with_element_totals(
            reduced_problem,
            Some(reduced_totals),
        )?;
        let mut runner = PreparedEquilibriumRunner::new(prepared, reduced_gibbs_sym)?;
        runner.configure().clone_from(&self.solver_settings);
        // Independent K_eq validation is applied to canonical fixed-species
        // problems. It must not silently reinterpret an active-set projection.
        runner.configure().keq_validation_mode = EquilibriumConstantValidationMode::Off;
        let outcome = runner.solve()?;
        let reduced = outcome.solution;
        let floor_log = PHASE_CONTROL_TRACE_MOLE_FLOOR.ln();
        let log_moles = projection.scatter_log_moles(reduced.log_moles(), floor_log)?;
        let moles = log_moles
            .iter()
            .map(|value| value.exp())
            .collect::<Vec<_>>();
        let mole_table = self
            .subs_data
            .substances
            .iter()
            .cloned()
            .zip(moles.iter().copied())
            .map(|(species, moles)| (species, vec![moles]))
            .collect();

        Ok(EquilibriumSolveCandidate {
            log_moles,
            moles,
            mole_table,
            validation_report: reduced.validation().clone(),
            solve_report: outcome.solve_report,
            keq_validation_status: None,
        })
    }

    #[deprecated(
        note = "use phase_equilibrium_workflow::solve_resolved_pt; this mutable method is compatibility-only"
    )]
    pub fn solve_with_phase_control(&mut self) -> Result<(), ReactionExtentError> {
        let keq_applicable =
            self.phases.len() == 1 && matches!(self.phases[0].kind, PhaseActivityModel::IdealGas);
        if self.solver_settings.keq_validation_mode == EquilibriumConstantValidationMode::Required
            && !keq_applicable
        {
            return Err(ReactionExtentError::ValidationNotApplicable {
                path: "equilibrium_constant_cross_validation",
                message: "Required K_eq validation supports only one ideal-gas phase; bounded phase control with condensed phases cannot satisfy this contract".to_string(),
            });
        }
        let species_phase = self.species_phase.clone();
        let (mut seed, _) = self.resolved_initial_guess()?;
        let derived_activity = if self.phase_active_mask.len() == self.phases.len() {
            self.phase_active_mask.clone()
        } else {
            initial_phase_activity_from_moles(
                &self.n0,
                &species_phase,
                self.phases.len(),
                self.phase_manager.phase_eps,
            )?
        };
        let max_phase_iterations = self.phase_manager.max_phase_iterations;
        if max_phase_iterations == 0 {
            return Err(ReactionExtentError::PhaseControlDidNotConverge { iterations: 0 });
        }
        let mut phase_set =
            PhaseSet::from_policy(&self.phase_manager.initial_phase_set, &derived_activity)?;
        // Normalize the user policy at the physical boundary before entering
        // the strictly-positive log-moles solver. This is not phase removal:
        // inventory-free candidates remain eligible for stability activation.
        phase_set.normalize_for_positive_solver(&derived_activity)?;
        let mut visited = HashSet::new();
        let mut settled_initial = phase_set.clone();
        settled_initial.settle_transitions();
        let initial_phase_set = settled_initial.clone();
        visited.insert(settled_initial);
        let initial_active_phases = phase_set.active_phases()?;
        let mut transitions = Vec::new();
        let mut nonlinear_reports = Vec::new();

        for iteration in 0..max_phase_iterations {
            let transition_started = Instant::now();
            phase_set.settle_transitions();
            let phase_active = phase_set.active_mask();
            let mut candidate = self.solve_fixed_active_set_candidate(&phase_active, &seed)?;
            if self.solver_settings.keq_validation_mode
                == EquilibriumConstantValidationMode::WhenApplicable
                && !keq_applicable
            {
                candidate.keq_validation_status = Some(
                    EquilibriumConstantCrossValidationStatus::ValidatorNotApplicable {
                        message: "bounded phase control includes more than one phase or a non-ideal phase; independent K_eq validation is limited to one ideal-gas phase".to_string(),
                    },
                );
            }
            nonlinear_reports.push(candidate.solve_report.clone());
            let phase_manager = &self.phase_manager;
            let mut y = candidate.log_moles.clone();
            // 2. Compute phase totals
            let n_phase = compute_phase_totals(&y, &species_phase);
            let stability = compute_phase_stability_reports(
                &y,
                &self.gibbs,
                &self.phases,
                &species_phase,
                &self.elem_composition,
                self.T,
                self.P,
                self.p0,
                &phase_set,
            )?;
            let transition_plan = phase_manager
                .classify_phases_at_temperature(self.T, &n_phase, &stability, &phase_set)?;
            let minimum_tpds = stability
                .iter()
                .map(|report| report.minimum_tpd)
                .collect::<Vec<_>>();

            match transition_plan {
                PhaseTransitionPlan::Deactivate { phase } => {
                    let phase_index = phase.index();
                    let previous_phase_set = phase_set.clone();
                    let minimum_tpd = stability[phase_index].minimum_tpd.ok_or_else(|| {
                        ReactionExtentError::InvalidCandidate {
                            field: "phase_stability",
                            message: format!(
                                "phase {phase_index} was selected for deactivation without a TPD minimum"
                            ),
                        }
                    })?;
                    let phase_moles = n_phase[phase_index];
                    info!("phase {phase_index} destroyed");
                    deactivate_phases_seed_only(
                        y.as_mut_slice(),
                        &[phase],
                        &species_phase,
                        PHASE_CONTROL_TRACE_MOLE_FLOOR,
                    )?;
                    phase_set.deactivate(phase);
                    let new_phase_set = phase_set.clone();
                    transitions.push(PhaseTransitionRecord {
                        iteration,
                        transition_duration: transition_started.elapsed(),
                        activated: Vec::new(),
                        deactivated: vec![phase],
                        phase_totals: n_phase,
                        minimum_tpds,
                        incipient_composition: None,
                        reason: PhaseTransitionReason::VanishingUnstableActivePhase {
                            phase_moles,
                            minimum_tpd,
                        },
                        previous_phase_set,
                        new_phase_set,
                        restart_seed: y.clone(),
                        nonlinear_report: candidate.solve_report.clone(),
                        candidate_validation: candidate.validation_report.clone(),
                    });
                    reject_repeated_phase_set(&mut visited, &phase_set, iteration + 1)?;
                    seed = y;
                    continue; // restart solver
                }
                PhaseTransitionPlan::Activate { phase } => {
                    let phase_index = phase.index();
                    let previous_phase_set = phase_set.clone();
                    let minimum_tpd = stability[phase_index].minimum_tpd.ok_or_else(|| {
                        ReactionExtentError::InvalidCandidate {
                            field: "phase_stability",
                            message: format!(
                                "phase {phase_index} was selected for activation without a TPD minimum"
                            ),
                        }
                    })?;
                    let composition = stability[phase_index]
                        .incipient_composition
                        .as_deref()
                        .ok_or_else(|| ReactionExtentError::InvalidCandidate {
                            field: "phase_stability",
                            message: format!(
                                "phase {phase_index} was selected for activation without a TPD minimizer composition"
                            ),
                        })?;
                    info!("phase {phase_index} created");
                    seed_activated_phase_with_composition(
                        y.as_mut_slice(),
                        phase,
                        &species_phase,
                        composition,
                        PhaseTotalSeedPolicy::RelativeToSystemTotal {
                            fraction: 1e-8,
                            minimum: PHASE_CONTROL_TRACE_MOLE_FLOOR,
                        },
                    )?;
                    phase_set.activate(phase);
                    let new_phase_set = phase_set.clone();
                    transitions.push(PhaseTransitionRecord {
                        iteration,
                        transition_duration: transition_started.elapsed(),
                        activated: vec![phase],
                        deactivated: Vec::new(),
                        phase_totals: n_phase,
                        minimum_tpds,
                        incipient_composition: Some(composition.to_vec()),
                        reason: PhaseTransitionReason::UnstableInactivePhase { minimum_tpd },
                        previous_phase_set,
                        new_phase_set,
                        restart_seed: y.clone(),
                        nonlinear_report: candidate.solve_report.clone(),
                        candidate_validation: candidate.validation_report.clone(),
                    });
                    reject_repeated_phase_set(&mut visited, &phase_set, iteration + 1)?;
                    seed = y;
                    continue; // restart solver
                }
                PhaseTransitionPlan::Hold { phase } => {
                    info!("phase {} held by hysteresis", phase.index());
                    validate_phase_set_candidate(&n_phase, &phase_active, phase_manager.phase_eps)?;
                    let report = PhaseControlledSolveReport {
                        iterations: iteration + 1,
                        initial_active_phases,
                        final_active_phases: phase_set.active_phases()?,
                        initial_phase_set,
                        final_phase_set: phase_set.clone(),
                        transitions,
                        final_validation: candidate.validation_report.clone(),
                        nonlinear_reports,
                    };
                    self.phase_active_mask = phase_active;
                    self.initial_guess = Some(y);
                    self.publish_solve_candidate(candidate);
                    self.last_phase_control_report = Some(report);
                    return Ok(());
                }
                PhaseTransitionPlan::NoTransition { reason } => {
                    info!("phase set stable: {reason:?}");
                    validate_phase_set_candidate(&n_phase, &phase_active, phase_manager.phase_eps)?;
                    let report = PhaseControlledSolveReport {
                        iterations: iteration + 1,
                        initial_active_phases,
                        final_active_phases: phase_set.active_phases()?,
                        initial_phase_set,
                        final_phase_set: phase_set.clone(),
                        transitions,
                        final_validation: candidate.validation_report.clone(),
                        nonlinear_reports,
                    };
                    // Nothing changed → converged equilibrium.
                    self.phase_active_mask = phase_active;
                    self.initial_guess = Some(y);
                    self.publish_solve_candidate(candidate);
                    self.last_phase_control_report = Some(report);
                    return Ok(());
                }
            }
        }

        Err(ReactionExtentError::PhaseControlDidNotConverge {
            iterations: max_phase_iterations,
        })
    }
}
///////////////////////////SYMBOLIC//////////////////////////////////////////////////////////////////////
/// Generates symbolic residual expressions for equilibrium equations
///
/// Creates symbolic expressions for residuals that can be used for
/// analytical differentiation or symbolic manipulation.
pub fn multiphase_equilibrium_residual_generator_sym(
    reactions: DMatrix<f64>,  // m × r
    elements: DMatrix<f64>,   // m × E
    element_totals: Vec<f64>, // E
    gibbs_sym: Vec<Expr>,     // m
    phases: Vec<Phase>,

    pressure: f64,
    p0: f64,
) -> Result<Vec<Expr>, ReactionExtentError> {
    multiphase_equilibrium_residual_generator_sym_at_temperature(
        reactions,
        elements,
        element_totals,
        gibbs_sym,
        phases,
        pressure,
        p0,
        Expr::Var("T".to_string()),
    )
}

/// Builds the same symbolic fixed-phase residual at an explicit temperature
/// expression.
///
/// Fixed-`P,T` uses the `T` equation parameter through the public wrapper
/// above. Coupled `P,H` instead supplies its bounded temperature expression,
/// so RST can differentiate the full system with respect to the additional
/// temperature coordinate rather than receiving a fixed-`P,T` Jacobian.
pub(crate) fn multiphase_equilibrium_residual_generator_sym_at_temperature(
    reactions: DMatrix<f64>,
    elements: DMatrix<f64>,
    element_totals: Vec<f64>,
    gibbs_sym: Vec<Expr>,
    phases: Vec<Phase>,
    pressure: f64,
    p0: f64,
    temperature: Expr,
) -> Result<Vec<Expr>, ReactionExtentError> {
    let m = reactions.nrows();
    let r = reactions.ncols();
    let e = elements.ncols();

    if gibbs_sym.len() != m || element_totals.len() != e {
        return Err(ReactionExtentError::DimensionMismatch(format!(
            "symbolic residual generator has {m} species, {} Gibbs expressions, and {} element totals for {e} elements",
            gibbs_sym.len(),
            element_totals.len(),
        )));
    }

    let species_phase = species_to_phase_map(&phases, m)?;
    let delta_n = reaction_phase_stoichiometry(&reactions, &phases);
    let rt = Expr::Const(R) * temperature;

    // log-mole symbolic variables (y_i = ln n_i)
    let y = Expr::IndexedVars(m, "y").0;
    if y.len() != m {
        return Err(ReactionExtentError::ResidualEvaluation(
            "y length mismatch".to_string(),
        ));
    }

    // --- species mole numbers as n_i = exp(y_i) ---
    let mut n: Vec<Expr> = Vec::with_capacity(m);
    for i in 0..m {
        n.push(y[i].clone().exp());
    }

    // --- phase mole totals ---
    let mut n_phase: Vec<Expr> = vec![Expr::Const(0.0); phases.len()];
    for i in 0..m {
        let j = species_phase[i];
        n_phase[j] += n[i].clone();
    }

    // --- phase offsets from the same activity contract as numeric residuals ---
    let phi = phase_activity_models(&phases)
        .into_iter()
        .map(|model| model.log_phase_offset(pressure, p0).map(Expr::Const))
        .collect::<Result<Vec<_>, _>>()?;

    // Build residuals: r reaction eqs followed by e element balance eqs
    let mut f: Vec<Expr> = vec![Expr::Const(0.0); r + e];

    for k in 0..r {
        let mut sum_ln_n = Expr::Const(0.0);
        let mut sum_ln_n_phase = Expr::Const(0.0);
        let mut sum_phi = Expr::Const(0.0);
        let mut dg0 = Expr::Const(0.0);

        // species contributions: ν_{ik} * y_i (since y_i = ln n_i) and ν_{ik} * g0_i(T)
        for i in 0..m {
            let nu = reactions[(i, k)];
            if nu == 0.0 {
                continue;
            }
            let nu_expr = Expr::Const(nu);
            sum_ln_n += nu_expr.clone() * y[i].clone();
            dg0 += nu_expr.clone() * gibbs_sym[i].clone();
        }

        // phase contributions: Δn_{kj} * ln(N_j) and Δn_{kj} * φ_j
        for j in 0..phases.len() {
            let dnk = delta_n[k][j];
            if dnk == 0.0 {
                continue;
            }
            let dnk_expr = Expr::Const(dnk);
            let nj = n_phase[j].clone();
            sum_ln_n_phase += dnk_expr.clone() * nj.ln();
            sum_phi += dnk_expr.clone() * phi[j].clone();
        }

        let ln_k = -dg0 / rt.clone();
        f[k] = (sum_ln_n - sum_ln_n_phase + sum_phi - ln_k).simplify();
    }

    // element balance equations: sum_i a_{i,el} * n_i - element_totals[el]
    for el in 0..e {
        let mut sum_el = Expr::Const(0.0);
        for i in 0..m {
            sum_el += Expr::Const(elements[(i, el)]) * n[i].clone();
        }
        f[r + el] = (sum_el - Expr::Const(element_totals[el])).simplify();
    }

    Ok(f)
}

////////////////////////CONVENIENCE FUNCTIONS////////////////////////////////
/// Convenience function for gas-phase equilibrium calculations
///
/// Sets up complete equilibrium solver for gas-phase systems with
/// automatic database searching and Gibbs function generation.
///
/// # Arguments
/// * `subs` - List of substance names
/// * `T` - Temperature in K
/// * `P` - Pressure in Pa
/// * `solver` - Choice of numerical solver
/// * `loglevel` - Logging level (None to disable)
/// * `scaling` - Whether to enable equation scaling
/// the function don't run calculation by itself only prepare data
/// for calculation
#[deprecated(
    note = "use ChemEquilibrium::prelude::solve_resolved_pt; retained for compatibility characterization"
)]
pub fn gas_solver(
    subs: Vec<String>,
    T: f64,
    P: f64,
    solver: Solvers,
    loglevel: Option<&str>,
    scaling: bool,
) -> Result<EquilibriumLogMoles, ReactionExtentError> {
    let mut user_subs = SubsData::new();

    user_subs.substances = subs;
    let map_of_phases = user_subs
        .substances
        .iter()
        .map(|s| (s.clone(), Some(Phases::Gas)))
        .collect();
    user_subs.map_of_phases = map_of_phases;
    user_subs.set_library_priority("NASA_gas".to_string(), LibraryPriority::Priority);
    user_subs.set_library_priority("NUIG_thermo".to_string(), LibraryPriority::Permitted);

    prepare_thermochemistry(&mut user_subs, T)?;
    let g_vec = collect_gibbs_functions(&mut user_subs)?;
    let element_composition = required_element_composition(&user_subs)?;
    let mut instance = EquilibriumLogMoles::new();

    instance.elem_composition = element_composition;
    instance.subs_data = user_subs;
    instance.with_loglevel(loglevel);
    instance.solver_settings.solver = solver;
    // This compatibility helper accepts a legacy `Solvers` selector. Make
    // that contract explicit instead of letting symbolic context silently
    // promote the call into the RST production cascade.
    instance.solver_settings.solver_policy = Some(SolverPolicy::legacy_default(solver));
    instance.solver_settings.scaling_flag = scaling;
    instance.gibbs = g_vec;

    instance.P = P;
    instance.T = T;
    Ok(instance)
}

#[deprecated(
    note = "use ChemEquilibrium::prelude::solve_resolved_pt; retained for compatibility characterization"
)]
pub fn gas_solver_from_elements(
    elements: Vec<String>,
    map_of_nonzero_moles: HashMap<String, f64>,
    T: f64,
    P: f64,
    solver: Solvers,
    loglevel: Option<&str>,
    scaling: bool,
) -> Result<EquilibriumLogMoles, ReactionExtentError> {
    let mut user_subs = SubsData::new();

    user_subs.search_by_elements_only(elements)?;
    let subs = user_subs.substances.clone();
    user_subs.parse_all_thermal_coeffs()?;
    user_subs.extract_all_thermal_coeffs(T)?;
    user_subs.calculate_elem_composition_and_molar_mass(None)?;
    let g_vec = collect_gibbs_functions(&mut user_subs)?;
    let element_composition = required_element_composition(&user_subs)?;
    let mut instance = EquilibriumLogMoles::new();
    instance.with_loglevel(loglevel);
    instance.elem_composition = element_composition;
    instance.subs_data = user_subs;
    info!("calculating n0");
    instance.set_n0_from_non_zero_map(map_of_nonzero_moles)?;
    let n0 = instance.n0.clone();

    instance.solver_settings.solver = solver;
    instance.solver_settings.solver_policy = Some(SolverPolicy::legacy_default(solver));
    instance.solver_settings.scaling_flag = scaling;
    instance.gibbs = g_vec;
    let species: Vec<usize> = (0..subs.len()).map(|x| x as usize).collect();

    instance.set_problem(
        n0,
        vec![Phase {
            kind: PhaseActivityModel::IdealGas,
            species: species,
        }],
        P,
    )?;
    info!(" phase data{:?}", instance.phases);
    instance.P = P;
    instance.T = T;

    Ok(instance)
}

#[deprecated(
    note = "use the typed temperature-postprocessing workflow after solve_resolved_pt; retained for compatibility"
)]
pub fn gas_solver_for_T_range(
    subs: Vec<String>,
    n0: Vec<f64>,
    P: f64,
    T_start: f64,
    T_end: f64,
    T_step: f64,
    solver: Solvers,
    loglevel: Option<&str>,
    parallel: bool,
) -> Result<EquilibriumLogMoles, ReactionExtentError> {
    if subs.len() != n0.len() {
        return Err(ReactionExtentError::DimensionMismatch(format!(
            "gas temperature-range workflow received {} species and {} initial-mole values",
            subs.len(),
            n0.len()
        )));
    }
    let mut user_subs = SubsData::new();

    user_subs.substances = subs.clone();
    let map_of_phases = user_subs
        .substances
        .iter()
        .map(|s| (s.clone(), Some(Phases::Gas)))
        .collect();
    user_subs.map_of_phases = map_of_phases;

    user_subs.set_library_priority("NASA_gas".to_string(), LibraryPriority::Priority);
    user_subs.set_library_priority("NUIG_thermo".to_string(), LibraryPriority::Permitted);

    prepare_thermochemistry(&mut user_subs, T_start)?;
    let g_vec = collect_gibbs_functions(&mut user_subs)?;
    let element_composition = required_element_composition(&user_subs)?;
    let species = user_subs.substances.clone();
    let initial_log_moles = LogMolesInitialGuess::from_initial_moles(&n0)?;
    let conditions = EquilibriumConditions::new(T_start, P, P)?;
    let problem = EquilibriumProblem::new(
        species,
        n0,
        initial_log_moles,
        element_composition,
        g_vec,
        vec![Phase {
            kind: PhaseActivityModel::IdealGas,
            species: (0..user_subs.substances.len()).collect(),
        }],
        conditions,
    )?;
    let mut instance = EquilibriumLogMoles::from_problem(problem)?;
    instance.subs_data = user_subs;
    instance.solver_settings.scaling_flag = false;
    instance.with_loglevel(loglevel);
    instance.solver_settings.solver = solver;
    // Preserve the historical meaning of this compatibility selector: the
    // requested legacy method runs first, followed only by legacy fallbacks.
    instance.solver_settings.solver_policy = Some(SolverPolicy::legacy_default(solver));
    if parallel {
        instance.solve_for_T_range_par2(T_start, T_end, T_step)?
    } else {
        instance.solve_for_T_range(T_start, T_end, T_step)?
    };
    Ok(instance)
}

#[deprecated(
    note = "use the typed temperature-postprocessing workflow after solve_resolved_pt; retained for compatibility"
)]
pub fn gas_solver_for_T_range_for_elements(
    elements: Vec<String>,
    map_of_nonzero_moles: HashMap<String, f64>,
    P: f64,
    T_start: f64,
    T_end: f64,
    T_step: f64,
    solver: Solvers,
    loglevel: Option<&str>,
) -> Result<EquilibriumLogMoles, ReactionExtentError> {
    let mut user_subs = SubsData::new();
    user_subs.search_by_elements_only(elements)?;
    user_subs.parse_all_thermal_coeffs()?;
    user_subs.calculate_elem_composition_and_molar_mass(None)?;
    let subs = user_subs.substances.clone();
    let map_of_phases = subs
        .clone()
        .iter()
        .map(|s| (s.clone(), Some(Phases::Gas)))
        .collect();
    user_subs.map_of_phases = map_of_phases;
    prepare_thermochemistry(&mut user_subs, T_start)?;
    let g_vec = collect_gibbs_functions(&mut user_subs)?;

    let n0 = {
        let subs = user_subs.substances.clone();
        let mut n0 = vec![1e-8; subs.len()];
        for (i, subs_i) in subs.iter().enumerate() {
            if let Some(moles_i) = map_of_nonzero_moles.get(subs_i) {
                n0[i] = *moles_i;
            }
        }
        n0
    };
    let element_composition = required_element_composition(&user_subs)?;
    let initial_log_moles = LogMolesInitialGuess::from_initial_moles(&n0)?;
    let conditions = EquilibriumConditions::new(T_start, P, P)?;
    let problem = EquilibriumProblem::new(
        user_subs.substances.clone(),
        n0,
        initial_log_moles,
        element_composition,
        g_vec,
        vec![Phase {
            kind: PhaseActivityModel::IdealGas,
            species: (0..user_subs.substances.len()).collect(),
        }],
        conditions,
    )?;
    let mut instance = EquilibriumLogMoles::from_problem(problem)?;
    instance.subs_data = user_subs;
    instance.solver_settings.scaling_flag = false;
    instance.with_loglevel(loglevel);
    instance.solver_settings.solver = solver;
    instance.solver_settings.solver_policy = Some(SolverPolicy::legacy_default(solver));
    instance.solve_for_T_range(T_start, T_end, T_step)?;
    Ok(instance)
}

////////////////////////misc///////////////////////////////
/// Computes finite difference approximation of Jacobian matrix
///
/// Used for testing or when analytical Jacobian is not available.
/// Less accurate and slower than analytical Jacobian.
pub fn finite_difference_jacobian<F>(
    f: &F,
    x: &[f64],
    eps: f64,
) -> Result<DMatrix<f64>, ReactionExtentError>
where
    F: Fn(&[f64]) -> Result<Vec<f64>, ReactionExtentError>,
{
    let f0 = f(x)?;
    let n = x.len();
    let m = f0.len();

    let mut j = DMatrix::<f64>::zeros(m, n);

    for k in 0..n {
        let mut x_pert = x.to_vec();
        x_pert[k] += eps;

        let f1 = f(&x_pert)?;

        for i in 0..m {
            j[(i, k)] = (f1[i] - f0[i]) / eps;
        }
    }

    Ok(j)
}
