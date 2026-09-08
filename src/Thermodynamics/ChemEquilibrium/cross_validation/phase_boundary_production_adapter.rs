//! Production evidence adapters for pure-phase boundary cross-validation.
//!
//! The independent `K_eq` boundary validator deliberately knows nothing about
//! the production outer loop. This module is the narrow, one-way adapter that
//! reads an immutable accepted solution and turns its published lifecycle
//! evidence into [`CanonicalPurePhaseEvidence`]. It never reruns a solver,
//! reconstructs a phase mask, or treats numerical trace moles as physical
//! inventory.
//!
//! Two sources of boundary evidence are intentionally distinct:
//!
//! - a stable inactive candidate uses its final evaluated TPD report;
//! - an appearing candidate uses the TPD stored on its activation transition.
//! - a disappearing candidate uses the reduced-boundary TPD stored on its
//!   `BoundaryUnstableActivePhase` transition.
//!
//! A final phase report is not a substitute for the gas-only boundary TPD that
//! caused an activation or deactivation decision.

use crate::Thermodynamics::ChemEquilibrium::equilibrium_activity::PhaseActivityModel;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_nonlinear::ReactionExtentError;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_workflows::{
    PhaseStatus, PhaseTransitionReason, PhaseTransitionRecord,
};
use crate::Thermodynamics::ChemEquilibrium::phase_boundary_validation::{
    CanonicalPurePhaseEvidence, canonical_pure_phase_evidence_from_stability_report,
};
use crate::Thermodynamics::ChemEquilibrium::phase_equilibrium_solution::MultiphaseEquilibriumSolution;
use crate::Thermodynamics::User_PhaseOrSolution::PhaseModel;
use crate::Thermodynamics::phase_layout::PhaseId;
use std::ops::Range;

/// Identifies the gas reference assemblage and one pure condensed candidate
/// within an accepted production solution.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PurePhaseProductionEvidenceRequest {
    pub gas_phase: PhaseId,
    pub candidate_phase: PhaseId,
}

impl PurePhaseProductionEvidenceRequest {
    /// Creates an explicit request. The referenced phases are validated only
    /// against a concrete accepted solution because their component topology
    /// is solution-layout dependent.
    pub fn new(gas_phase: PhaseId, candidate_phase: PhaseId) -> Self {
        Self {
            gas_phase,
            candidate_phase,
        }
    }
}

/// Extracts canonical evidence for a candidate that remained inactive after
/// bounded phase control.
///
/// The final evaluated stability report is a valid gas-only boundary witness
/// only because the candidate is inactive in the accepted solution.
pub fn stable_inactive_evidence_from_solution(
    solution: &MultiphaseEquilibriumSolution,
    request: &PurePhaseProductionEvidenceRequest,
) -> Result<CanonicalPurePhaseEvidence, ReactionExtentError> {
    let layout = pure_phase_layout(solution, request)?;
    if layout.candidate_active {
        return Err(not_applicable(
            "stable-inactive evidence requested for an active candidate phase",
        ));
    }

    let acceptance = solution.acceptance_report().ok_or_else(|| {
        not_applicable("production pure-phase evidence requires bounded phase-control acceptance")
    })?;
    let stability = acceptance
        .phase_stability
        .iter()
        .find(|report| report.phase == layout.candidate_phase_index)
        .ok_or_else(|| {
            not_applicable("final acceptance report lacks the candidate stability report")
        })?;

    canonical_pure_phase_evidence_from_stability_report(
        layout.candidate_phase_index.index(),
        layout.gas_species,
        layout.candidate_name,
        false,
        layout.gas_moles,
        layout.candidate_moles,
        stability,
    )
}

/// Extracts canonical evidence for a candidate activated by bounded phase
/// control.
///
/// The boundary TPD comes from the recorded activation transition, while the
/// final topology and composition come from the immutable accepted solution.
/// This preserves the physical ordering of the lifecycle argument.
pub fn activation_evidence_from_solution(
    solution: &MultiphaseEquilibriumSolution,
    request: &PurePhaseProductionEvidenceRequest,
) -> Result<CanonicalPurePhaseEvidence, ReactionExtentError> {
    let layout = pure_phase_layout(solution, request)?;
    if !layout.candidate_active {
        return Err(not_applicable(
            "activation evidence requested for a candidate that is not active in the accepted solution",
        ));
    }

    let phase_control = solution.phase_control_report().ok_or_else(|| {
        not_applicable("activation evidence requires a bounded phase-control report")
    })?;
    let transition = phase_control
        .transitions
        .iter()
        .find(|record| record.activated.contains(&layout.candidate_phase_index))
        .ok_or_else(|| {
            not_applicable("accepted candidate has no recorded activation transition")
        })?;

    let (previously_active, now_active) = transition_phase_direction(
        transition,
        layout.candidate_phase_index.index(),
        "activation",
    )?;
    if previously_active || !now_active {
        return Err(not_applicable(
            "activation transition does not change the candidate phase from inactive to active",
        ));
    }

    let boundary_tpd = transition
        .minimum_tpds
        .get(layout.candidate_phase_index.index())
        .and_then(|value| *value)
        .filter(|value| value.is_finite())
        .ok_or_else(|| not_applicable("activation transition lacks a finite candidate TPD"))?;
    let boundary_gas_moles = boundary_gas_moles_from_transition(
        transition,
        layout.gas_component_range.clone(),
        "activation",
    )?;

    Ok(CanonicalPurePhaseEvidence {
        gas_species: layout.gas_species,
        candidate_name: layout.candidate_name,
        candidate_active: true,
        gas_moles: layout.gas_moles,
        boundary_gas_moles,
        candidate_moles: layout.candidate_moles,
        boundary_minimum_tpd: Some(boundary_tpd),
    })
}

/// Extracts canonical evidence for a pure candidate removed after validated
/// reduced-boundary recovery.
///
/// Disappearance is identified by the immutable lifecycle transition, not by
/// the candidate's final numerical trace coordinate. The transition must move
/// the requested candidate from active to inactive for the physical reason
/// [`PhaseTransitionReason::BoundaryUnstableActivePhase`]. Its reduced-set TPD
/// and restart seed are the boundary witness; final topology and physical
/// composition still come from the accepted solution.
pub fn disappearance_evidence_from_solution(
    solution: &MultiphaseEquilibriumSolution,
    request: &PurePhaseProductionEvidenceRequest,
) -> Result<CanonicalPurePhaseEvidence, ReactionExtentError> {
    let layout = pure_phase_layout(solution, request)?;
    if layout.candidate_active {
        return Err(not_applicable(
            "disappearance evidence requested for a candidate that is still active in the accepted solution",
        ));
    }

    let phase_control = solution.phase_control_report().ok_or_else(|| {
        not_applicable("disappearance evidence requires a bounded phase-control report")
    })?;
    let transition = phase_control
        .transitions
        .iter()
        .find(|record| record.deactivated.contains(&layout.candidate_phase_index))
        .ok_or_else(|| {
            not_applicable("inactive candidate has no recorded deactivation transition")
        })?;

    if !transition.activated.is_empty()
        || transition.deactivated.as_slice() != [layout.candidate_phase_index]
    {
        return Err(not_applicable(
            "disappearance transition must contain exactly one deactivated phase and no activation",
        ));
    }
    let (previously_active, now_active) = transition_phase_direction(
        transition,
        layout.candidate_phase_index.index(),
        "disappearance",
    )?;
    if !previously_active || now_active {
        return Err(not_applicable(
            "disappearance transition does not change the candidate phase from active to inactive",
        ));
    }

    let reason_tpd = match transition.reason {
        PhaseTransitionReason::BoundaryUnstableActivePhase { minimum_tpd, .. }
            if minimum_tpd.is_finite() =>
        {
            minimum_tpd
        }
        PhaseTransitionReason::BoundaryUnstableActivePhase { .. } => {
            return Err(not_applicable(
                "boundary-recovery disappearance reason contains a non-finite TPD",
            ));
        }
        _ => {
            return Err(not_applicable(
                "candidate deactivation was not caused by validated reduced-boundary recovery",
            ));
        }
    };
    let recorded_tpd = transition
        .minimum_tpds
        .get(layout.candidate_phase_index.index())
        .and_then(|value| *value)
        .filter(|value| value.is_finite())
        .ok_or_else(|| {
            not_applicable("disappearance transition lacks a finite candidate boundary TPD")
        })?;
    if recorded_tpd != reason_tpd {
        return Err(not_applicable(
            "disappearance transition reason and per-phase TPD evidence disagree",
        ));
    }

    let boundary_gas_moles = boundary_gas_moles_from_transition(
        transition,
        layout.gas_component_range.clone(),
        "disappearance",
    )?;

    Ok(CanonicalPurePhaseEvidence {
        gas_species: layout.gas_species,
        candidate_name: layout.candidate_name,
        candidate_active: false,
        gas_moles: layout.gas_moles,
        boundary_gas_moles,
        candidate_moles: layout.candidate_moles,
        boundary_minimum_tpd: Some(recorded_tpd),
    })
}

fn transition_phase_direction(
    transition: &PhaseTransitionRecord,
    phase_position: usize,
    route: &str,
) -> Result<(bool, bool), ReactionExtentError> {
    let previous_mask = transition.previous_phase_set.active_mask();
    let new_mask = transition.new_phase_set.active_mask();
    let previously_active = previous_mask.get(phase_position).copied().ok_or_else(|| {
        not_applicable(format!(
            "{route} transition previous phase set does not match the accepted layout"
        ))
    })?;
    let now_active = new_mask.get(phase_position).copied().ok_or_else(|| {
        not_applicable(format!(
            "{route} transition new phase set does not match the accepted layout"
        ))
    })?;
    Ok((previously_active, now_active))
}

/// Reconstructs the boundary gas composition from a transition restart seed.
///
/// The transition's log-mole restart seed carries the gas-only boundary state
/// that justified the activation or deactivation decision. This helper selects
/// the gas component range, exponentiates the log-moles into physical moles,
/// and rejects any non-finite or non-positive result.
fn boundary_gas_moles_from_transition(
    transition: &PhaseTransitionRecord,
    gas_component_range: Range<usize>,
    route: &str,
) -> Result<Vec<f64>, ReactionExtentError> {
    let boundary_log_moles = transition
        .restart_seed
        .get(gas_component_range)
        .ok_or_else(|| {
            not_applicable(format!(
                "{route} restart seed does not match the gas layout"
            ))
        })?;
    let boundary_gas_moles = boundary_log_moles
        .iter()
        .map(|log_moles| log_moles.exp())
        .collect::<Vec<_>>();
    if boundary_gas_moles
        .iter()
        .any(|moles| !moles.is_finite() || *moles <= 0.0)
    {
        return Err(not_applicable(format!(
            "{route} restart seed contains non-finite or non-positive gas moles"
        )));
    }
    Ok(boundary_gas_moles)
}

struct PurePhaseProductionLayout {
    gas_species: Vec<String>,
    gas_moles: Vec<f64>,
    gas_component_range: Range<usize>,
    candidate_name: String,
    candidate_moles: f64,
    candidate_active: bool,
    candidate_phase_index: crate::Thermodynamics::ChemEquilibrium::equilibrium_ids::PhaseIndex,
}

/// Extracts the gas/candidate layout from an accepted production solution.
///
/// Locates the requested ideal-gas reference phase and the one-component pure
/// condensed candidate phase in the accepted layout, projects their physical
/// mole numbers and component ranges, and derives whether the candidate is
/// active from its published lifecycle status. Any mismatch between the
/// request and the accepted solution (missing phase, wrong activity model,
/// infeasible moles) is rejected as validation-not-applicable.
fn pure_phase_layout(
    solution: &MultiphaseEquilibriumSolution,
    request: &PurePhaseProductionEvidenceRequest,
) -> Result<PurePhaseProductionLayout, ReactionExtentError> {
    let gas = solution
        .phases()
        .iter()
        .find(|phase| phase.id() == &request.gas_phase)
        .ok_or_else(|| {
            not_applicable("requested gas phase is absent from the accepted solution")
        })?;
    if gas.activity_model() != PhaseActivityModel::IdealGas {
        return Err(not_applicable(
            "pure-phase boundary validation requires an ideal-gas reference phase",
        ));
    }

    let candidate = solution
        .phases()
        .iter()
        .find(|phase| phase.id() == &request.candidate_phase)
        .ok_or_else(|| {
            not_applicable("requested candidate phase is absent from the accepted solution")
        })?;
    if candidate.phase_model() != PhaseModel::PureCondensed
        || candidate.activity_model() != PhaseActivityModel::IdealSolution
        || candidate.component_range().len() != 1
    {
        return Err(not_applicable(
            "pure-phase boundary validation requires a one-component pure condensed candidate",
        ));
    }

    let gas_range = gas.component_range();
    let candidate_index = candidate.component_range().start;
    let components = solution.metadata().components();
    let physical_moles = solution.component_moles();
    if physical_moles.len() != components.len() || candidate_index >= physical_moles.len() {
        return Err(not_applicable(
            "accepted component moles do not match the production phase layout",
        ));
    }

    let gas_species = components[gas_range.clone()]
        .iter()
        .map(|component| component.substance().to_string())
        .collect::<Vec<_>>();
    let gas_moles = physical_moles[gas_range].to_vec();
    let candidate_name = components[candidate_index].substance().to_string();
    let candidate_moles = physical_moles[candidate_index];
    if gas_species.is_empty()
        || gas_moles
            .iter()
            .any(|moles| !moles.is_finite() || *moles < 0.0)
        || !candidate_moles.is_finite()
        || candidate_moles < 0.0
    {
        return Err(not_applicable(
            "accepted solution has non-finite or negative pure-phase evidence moles",
        ));
    }

    let candidate_status = solution
        .phase_status(&request.candidate_phase)
        .ok_or_else(|| not_applicable("accepted solution lacks the candidate phase status"))?;
    let candidate_active = matches!(
        candidate_status,
        PhaseStatus::Active | PhaseStatus::Appeared
    );

    Ok(PurePhaseProductionLayout {
        gas_species,
        gas_moles,
        gas_component_range: gas.component_range(),
        candidate_name,
        candidate_moles,
        candidate_active,
        candidate_phase_index: candidate.index(),
    })
}

/// Constructs a `ValidationNotApplicable` error scoped to this adapter.
///
/// Centralizes the path label so every rejection emitted by the production
/// adapter points at the same diagnostic namespace and carries the same
/// user-facing message.
fn not_applicable(message: impl Into<String>) -> ReactionExtentError {
    ReactionExtentError::ValidationNotApplicable {
        path: "pure_phase_boundary_production_adapter",
        message: message.into(),
    }
}
