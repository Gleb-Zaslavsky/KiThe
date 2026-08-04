//! Immutable active-set orchestration for bounded fixed-`P,T` equilibrium.
//!
//! The numerical inner problem is prepared once per settled active phase set,
//! while the outer loop owns only seeds, phase transitions, and reports. Each
//! active-set entry also retains its RST symbolic problem, so a temperature
//! sweep updates the shared `T` parameter instead of rebuilding equations and
//! lambdified closures. The runner never constructs the historical
//! `EquilibriumLogMoles` god object. That object remains available for
//! compatibility callers; this runner is the canonical bridge path.

use std::collections::{HashMap, HashSet};
use std::time::{Duration, Instant};

use crate::Thermodynamics::ChemEquilibrium::equilibrium_active_set::ActiveSetProjection;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_constant_cross_validation::EquilibriumConstantCrossValidationStatus;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_constant_validation::EquilibriumConstantValidationMode;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_ids::PhaseIndex;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_log_moles::{
    EquilibriumSolverSettings, GibbsFn,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_nonlinear::ReactionExtentError;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_prepared_runner::PreparedEquilibriumRunner;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_problem::{
    EquilibriumProblem, LogMolesInitialGuess, PreparedEquilibriumProblem,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_rst_backend::{
    prepare_rst_symbolic_problem_from_prepared, RstPreparedProblem,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_solver_policy::{
    EquilibriumSolveReport, SolverBackend,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_validation::EquilibriumCandidateReport;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_workflows::{
    build_multiphase_acceptance_report, compute_phase_stability_reports, compute_phase_totals,
    initial_phase_activity_from_moles, reject_repeated_phase_set, seed_activated_phase,
    validate_phase_set_candidate, MultiphaseAcceptanceReport, PhaseControlledSolveReport,
    PhaseManager, PhaseSeedPolicy, PhaseSet, PhaseStatus, PhaseTransitionPlan,
    PhaseTransitionReason, PhaseTransitionRecord, PHASE_CONTROL_TRACE_MOLE_FLOOR,
};
use RustedSciThe::symbolic::symbolic_engine::Expr;

/// Immutable-runner output consumed by the phase-aware public result facade.
#[derive(Debug)]
pub(crate) struct PreparedPhaseControlOutcome {
    pub(crate) solution:
        crate::Thermodynamics::ChemEquilibrium::equilibrium_problem::EquilibriumSolution,
    pub(crate) solve_report: EquilibriumSolveReport,
    pub(crate) keq_validation_status: Option<EquilibriumConstantCrossValidationStatus>,
    pub(crate) phase_control_report: PhaseControlledSolveReport,
    pub(crate) acceptance_report: MultiphaseAcceptanceReport,
    pub(crate) phase_statuses: Vec<PhaseStatus>,
    pub(crate) projection_build: Duration,
    /// Time spent constructing missing reduced formulations or rebuilding
    /// their symbolic backend for the current active-set snapshot.
    pub(crate) formulation_build: Duration,
    pub(crate) validation_duration: Duration,
    pub(crate) rst_symbolic_reused: bool,
}

pub(crate) struct PreparedActiveSetCandidate {
    pub(crate) log_moles: Vec<f64>,
    /// Active mask actually used by the fixed-set formulation. Normally it
    /// equals the lifecycle mask; a monolithic P,H probe may solve a wider
    /// trace-seeded mask to discover a phase whose enthalpy branch is needed.
    pub(crate) solved_active_mask: Vec<bool>,
    pub(crate) validation_report: EquilibriumCandidateReport,
    pub(crate) solve_report: EquilibriumSolveReport,
    pub(crate) keq_validation_status: Option<EquilibriumConstantCrossValidationStatus>,
    /// Full-layout Gibbs capabilities and conditions at this accepted
    /// candidate. Phase stability must use these values rather than the
    /// runner's construction temperature; P,H candidates solve temperature
    /// and composition together.
    pub(crate) stability_gibbs: Vec<GibbsFn>,
    pub(crate) conditions:
        crate::Thermodynamics::ChemEquilibrium::equilibrium_problem::EquilibriumConditions,
    pub(crate) projection_build: Duration,
    pub(crate) formulation_build: Duration,
    pub(crate) validation_duration: Duration,
    pub(crate) rst_symbolic_reused: bool,
}

struct PreparedBoundaryRecovery {
    phase: PhaseIndex,
    candidate: PreparedActiveSetCandidate,
    phase_set: PhaseSet,
    stability:
        Vec<crate::Thermodynamics::ChemEquilibrium::equilibrium_workflows::PhaseStabilityReport>,
}

/// Numeric and symbolic preparation retained for one active phase mask.
///
/// Canonical fixed-P,T RST graphs receive standard Gibbs values as equation
/// parameters. A coefficient-interval transition therefore reuses this entry;
/// only a changed active mask needs another graph.
struct PreparedActiveSetCacheEntry {
    prepared: PreparedEquilibriumProblem,
    symbolic_standard_gibbs: Vec<Expr>,
    rst_problem: Option<RstPreparedProblem>,
    /// Cumulative time spent building this active-set entry, including later
    /// symbolic rebuilds when a native coefficient interval changes.
    formulation_build: Duration,
}

/// Runs phase transitions around immutable prepared fixed-set solves.
pub(crate) struct PreparedPhaseControlRunner {
    /// Prepared topology and conserved inventory shared by every temperature
    /// point. Its retarget operation replaces only conditions, seed, and Gibbs
    /// closures while retaining the reaction basis.
    prepared: PreparedEquilibriumProblem,
    symbolic_standard_gibbs: Vec<Expr>,
    solver_settings: EquilibriumSolverSettings,
    phase_manager: PhaseManager,
    timing_enabled: bool,
    continuation_seed: Option<Vec<f64>>,
    continuation_phase_set: Option<PhaseSet>,
    projection_cache: HashMap<Vec<bool>, ActiveSetProjection>,
    prepared_active_set_cache: HashMap<Vec<bool>, PreparedActiveSetCacheEntry>,
}

impl PreparedPhaseControlRunner {
    /// Creates the runner from one validated fixed-`P,T` problem.
    pub(crate) fn new(
        problem: EquilibriumProblem,
        symbolic_standard_gibbs: Vec<Expr>,
        timing_enabled: bool,
    ) -> Result<Self, ReactionExtentError> {
        problem.validate()?;
        if !symbolic_standard_gibbs.is_empty()
            && symbolic_standard_gibbs.len() != problem.species().len()
        {
            return Err(ReactionExtentError::DimensionMismatch(format!(
                "symbolic Gibbs snapshot has {} entries for {} species",
                symbolic_standard_gibbs.len(),
                problem.species().len()
            )));
        }
        Ok(Self {
            prepared: PreparedEquilibriumProblem::new(problem)?,
            symbolic_standard_gibbs,
            solver_settings: EquilibriumSolverSettings::default(),
            phase_manager: PhaseManager::default(),
            timing_enabled,
            continuation_seed: None,
            continuation_phase_set: None,
            projection_cache: HashMap::new(),
            prepared_active_set_cache: HashMap::new(),
        })
    }

    /// Retargets numeric conditions and thermochemistry while retaining the
    /// symbolic capability snapshot captured during initial phase resolution.
    ///
    /// This is the temperature-range path: fixed-P,T RST graphs parameterize
    /// `G0_i`, so recalculating symbolic polynomials at every point is both
    /// redundant and a source of avoidable rebuilds.
    pub(crate) fn retarget_numeric(
        &mut self,
        conditions: crate::Thermodynamics::ChemEquilibrium::equilibrium_problem::EquilibriumConditions,
        seed: LogMolesInitialGuess,
        gibbs: Vec<crate::Thermodynamics::ChemEquilibrium::equilibrium_log_moles::GibbsFn>,
    ) -> Result<(), ReactionExtentError> {
        self.prepared = self
            .prepared
            .retarget_with_gibbs(conditions, seed.clone(), gibbs)?;
        self.continuation_seed = Some(seed.as_slice().to_vec());
        Ok(())
    }

    /// Carries the accepted phase state into the next temperature point.
    pub(crate) fn set_continuation_phase_set(
        &mut self,
        phase_set: PhaseSet,
    ) -> Result<(), ReactionExtentError> {
        phase_set.active_phases()?;
        if phase_set.active_mask().len() != self.prepared.problem().phases().len() {
            return Err(ReactionExtentError::DimensionMismatch(format!(
                "continuation phase set has {} phases for {} prepared phases",
                phase_set.active_mask().len(),
                self.prepared.problem().phases().len()
            )));
        }
        self.continuation_phase_set = Some(phase_set);
        Ok(())
    }

    /// Number of distinct active-set projections prepared so far.
    pub(crate) fn projection_cache_size(&self) -> usize {
        self.projection_cache.len()
    }

    /// Number of reduced formulations retained for active masks encountered
    /// during the sweep.
    pub(crate) fn prepared_active_set_cache_size(&self) -> usize {
        self.prepared_active_set_cache.len()
    }

    /// Number of active-set entries with a retained RST symbolic problem.
    pub(crate) fn rst_prepared_cache_size(&self) -> usize {
        self.prepared_active_set_cache
            .values()
            .filter(|entry| entry.rst_problem.is_some())
            .count()
    }

    /// Returns deterministic per-entry build timing for release reports.
    ///
    /// `HashMap` iteration order is deliberately hidden behind active-mask
    /// sorting so repeated runs produce comparable evidence. The mask uses
    /// the canonical phase order retained by the prepared problem.
    pub(crate) fn prepared_active_set_cache_timings(&self) -> Vec<(Vec<bool>, Duration)> {
        let mut entries = self
            .prepared_active_set_cache
            .iter()
            .map(|(active_mask, entry)| (active_mask.clone(), entry.formulation_build))
            .collect::<Vec<_>>();
        entries.sort_by(|left, right| left.0.cmp(&right.0));
        entries
    }

    /// Checks whether the configured solver policy selects a RustedSciThe
    /// backend for the supplied symbolic expression vector.
    ///
    /// Returns `true` when at least one backend in the resolved policy is
    /// a RustedSciThe variant and the symbolic expressions are non-empty.
    /// This guards the RST problem preparation path: when no symbolic backend
    /// is present, the runner skips the expensive symbolic construction and
    /// uses only the analytical formulation.
    fn uses_rst_backend(&self, symbols: &[Expr]) -> bool {
        self.solver_settings
            .solver_policy
            .as_ref()
            .map(|policy| {
                policy
                    .ordered_backends()
                    .iter()
                    .any(|backend| matches!(backend, SolverBackend::RustedSciThe(_)))
            })
            .unwrap_or(!symbols.is_empty())
    }

    /// Mutates numerical policy only; the equilibrium problem remains fixed.
    pub(crate) fn configure_solver(&mut self) -> &mut EquilibriumSolverSettings {
        &mut self.solver_settings
    }

    /// Mutates phase-control policy only; problem data remains fixed.
    pub(crate) fn configure_phase_control(&mut self) -> &mut PhaseManager {
        &mut self.phase_manager
    }

    /// Borrows the immutable full-layout problem used as the active-set
    /// projection source. The P,H adapter uses only its topology and conserved
    /// inventory; its nonlinear equations are built separately.
    pub(crate) fn prepared_problem(&self) -> &PreparedEquilibriumProblem {
        &self.prepared
    }

    /// Clones the validated backend settings for a fixed-set adapter.
    pub(crate) fn solver_settings(&self) -> EquilibriumSolverSettings {
        self.solver_settings.clone()
    }

    /// Executes the bounded outer loop and publishes only the final accepted
    /// immutable solution plus complete transition evidence.
    pub(crate) fn solve(&mut self) -> Result<PreparedPhaseControlOutcome, ReactionExtentError> {
        self.solve_with_fixed_active_solver(|runner, active, seed, species_phase, totals| {
            runner.solve_active_set(active, seed, species_phase, totals)
        })
    }

    /// Executes the common active-set lifecycle with an injected fixed-set
    /// candidate solver.
    ///
    /// The callback owns only one immutable fixed active set. This runner
    /// remains the sole owner of hysteresis, boundary recovery, cycle checks,
    /// transition accounting, and final transactional publication. P,T uses
    /// [`Self::solve_active_set`]; coupled P,H will use the same lifecycle
    /// after providing a monolithic candidate adapter.
    pub(crate) fn solve_with_fixed_active_solver<F>(
        &mut self,
        mut solve_active_set: F,
    ) -> Result<PreparedPhaseControlOutcome, ReactionExtentError>
    where
        F: FnMut(
            &mut Self,
            &[bool],
            &[f64],
            &[usize],
            &[f64],
        ) -> Result<PreparedActiveSetCandidate, ReactionExtentError>,
    {
        self.solver_settings.validate()?;
        self.phase_manager
            .validate_for_phase_count(self.prepared.problem().phases().len())?;

        let keq_applicable = self.prepared.problem().phases().len() == 1
            && matches!(
                self.prepared.problem().phases()[0].kind,
                crate::Thermodynamics::ChemEquilibrium::equilibrium_log_moles::PhaseKind::IdealGas
            );
        if self.solver_settings.keq_validation_mode == EquilibriumConstantValidationMode::Required
            && !keq_applicable
        {
            return Err(ReactionExtentError::ValidationNotApplicable {
                path: "equilibrium_constant_cross_validation",
                message: "Required K_eq validation supports only one ideal-gas phase; bounded phase control with condensed phases cannot satisfy this contract".to_string(),
            });
        }

        let species_phase =
            crate::Thermodynamics::ChemEquilibrium::equilibrium_log_moles::species_to_phase_map(
                self.prepared.problem().phases(),
                self.prepared.problem().species().len(),
            )?;
        let mut seed = self.continuation_seed.take().unwrap_or_else(|| {
            self.prepared
                .problem()
                .initial_log_moles()
                .as_slice()
                .to_vec()
        });
        let mut phase_set = if let Some(phase_set) = self.continuation_phase_set.take() {
            phase_set
        } else {
            let derived_activity = initial_phase_activity_from_moles(
                self.prepared.problem().initial_moles(),
                &species_phase,
                self.prepared.problem().phases().len(),
                self.phase_manager.phase_eps,
            )?;
            let mut phase_set =
                PhaseSet::from_policy(&self.phase_manager.initial_phase_set, &derived_activity)?;
            phase_set.normalize_for_positive_solver(&derived_activity)?;
            phase_set
        };
        let max_phase_iterations = self.phase_manager.max_phase_iterations;
        phase_set.settle_transitions();
        let initial_phase_set = phase_set.clone();
        let initial_active_phases = phase_set.active_phases()?;
        let mut visited = HashSet::new();
        visited.insert(phase_set.clone());
        let mut transitions = Vec::new();
        let mut nonlinear_reports = Vec::new();
        let mut projection_build = Duration::ZERO;
        let mut formulation_build = Duration::ZERO;
        let mut validation_duration = Duration::ZERO;
        let mut rst_symbolic_reused = false;
        let full_element_totals = self.prepared.element_totals().to_vec();

        for iteration in 0..max_phase_iterations {
            let transition_started = Instant::now();
            phase_set.settle_transitions();
            let phase_active = phase_set.active_mask();
            let mut candidate = match solve_active_set(
                self,
                &phase_active,
                &seed,
                &species_phase,
                &full_element_totals,
            ) {
                Ok(candidate) => {
                    projection_build += candidate.projection_build;
                    formulation_build += candidate.formulation_build;
                    validation_duration += candidate.validation_duration;
                    rst_symbolic_reused |= candidate.rst_symbolic_reused;
                    candidate
                }
                Err(primary_error) => {
                    let Some(recovery) = self.recover_active_boundary(
                        &phase_set,
                        &seed,
                        &species_phase,
                        &full_element_totals,
                        &mut solve_active_set,
                    )?
                    else {
                        return Err(primary_error);
                    };

                    let PreparedBoundaryRecovery {
                        phase,
                        mut candidate,
                        phase_set: recovered_phase_set,
                        stability,
                    } = recovery;
                    projection_build += candidate.projection_build;
                    formulation_build += candidate.formulation_build;
                    validation_duration += candidate.validation_duration;
                    rst_symbolic_reused |= candidate.rst_symbolic_reused;
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
                    let phase_totals = compute_phase_totals(&candidate.log_moles, &species_phase);
                    let driving_forces = stability
                        .iter()
                        .map(|report| report.driving_force)
                        .collect::<Vec<_>>();
                    let driving_force =
                        stability[phase.index()].driving_force.ok_or_else(|| {
                            ReactionExtentError::InvalidCandidate {
                                field: "phase_stability",
                                message: format!(
                                    "boundary recovery selected phase {} without a driving force",
                                    phase.index()
                                ),
                            }
                        })?;
                    let initial_phase_moles = self
                        .prepared
                        .problem()
                        .initial_moles()
                        .iter()
                        .zip(species_phase.iter())
                        .filter(|(_, phase_index)| **phase_index == phase.index())
                        .map(|(moles, _)| *moles)
                        .sum::<f64>();
                    let previous_phase_set = phase_set.clone();
                    phase_set = recovered_phase_set;
                    transitions.push(PhaseTransitionRecord {
                        iteration,
                        transition_duration: transition_started.elapsed(),
                        activated: Vec::new(),
                        deactivated: vec![phase],
                        phase_totals,
                        driving_forces,
                        reason: PhaseTransitionReason::BoundaryUnstableActivePhase {
                            initial_phase_moles,
                            driving_force,
                        },
                        previous_phase_set,
                        new_phase_set: phase_set.clone(),
                        restart_seed: candidate.log_moles.clone(),
                        nonlinear_report: candidate.solve_report.clone(),
                        candidate_validation: candidate.validation_report.clone(),
                    });
                    reject_repeated_phase_set(&mut visited, &phase_set, iteration + 1)?;
                    seed = candidate.log_moles.clone();
                    continue;
                }
            };
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
            let mut y = candidate.log_moles.clone();
            let phase_totals = compute_phase_totals(&y, &species_phase);
            let mut phase_active = phase_set.active_mask();
            let probe_stability = compute_phase_stability_reports(
                &y,
                &candidate.stability_gibbs,
                self.prepared.problem().phases(),
                &species_phase,
                self.prepared.problem().element_composition(),
                candidate.conditions.temperature(),
                candidate.conditions.pressure(),
                candidate.conditions.reference_pressure(),
                &phase_set,
            )?;
            if candidate.solved_active_mask.len() != phase_active.len() {
                return Err(ReactionExtentError::DimensionMismatch(format!(
                    "candidate solved-active mask has {} phases but lifecycle has {}",
                    candidate.solved_active_mask.len(),
                    phase_active.len()
                )));
            }
            if phase_active
                .iter()
                .zip(&candidate.solved_active_mask)
                .any(|(&lifecycle_active, &solved_active)| lifecycle_active && !solved_active)
            {
                return Err(ReactionExtentError::InvalidCandidate {
                    field: "phase_active_set",
                    message:
                        "a fixed-set candidate cannot silently remove a lifecycle-active phase"
                            .to_string(),
                });
            }
            let mut probe_expanded = false;
            for phase_index in 0..phase_active.len() {
                if phase_active[phase_index] || !candidate.solved_active_mask[phase_index] {
                    continue;
                }
                let phase = PhaseIndex::new(phase_index, phase_active.len())?;
                let driving_force =
                    probe_stability[phase_index].driving_force.ok_or_else(|| {
                        ReactionExtentError::InvalidCandidate {
                            field: "phase_stability",
                            message: format!(
                            "monolithic active-set probe expanded phase {} without a driving force",
                            phase.index()
                        ),
                        }
                    })?;
                let previous_phase_set = phase_set.clone();
                phase_set.activate(phase);
                phase_active = phase_set.active_mask();
                probe_expanded = true;
                transitions.push(PhaseTransitionRecord {
                    iteration,
                    transition_duration: transition_started.elapsed(),
                    activated: vec![phase],
                    deactivated: Vec::new(),
                    phase_totals: phase_totals.clone(),
                    driving_forces: probe_stability
                        .iter()
                        .map(|report| report.driving_force)
                        .collect(),
                    reason: PhaseTransitionReason::UnstableInactivePhase { driving_force },
                    previous_phase_set,
                    new_phase_set: phase_set.clone(),
                    restart_seed: y.clone(),
                    nonlinear_report: candidate.solve_report.clone(),
                    candidate_validation: candidate.validation_report.clone(),
                });
                reject_repeated_phase_set(&mut visited, &phase_set, iteration + 1)?;
            }
            let stability = if probe_expanded {
                compute_phase_stability_reports(
                    &y,
                    &candidate.stability_gibbs,
                    self.prepared.problem().phases(),
                    &species_phase,
                    self.prepared.problem().element_composition(),
                    candidate.conditions.temperature(),
                    candidate.conditions.pressure(),
                    candidate.conditions.reference_pressure(),
                    &phase_set,
                )?
            } else {
                probe_stability
            };
            let transition_plan = self.phase_manager.classify_phases_at_temperature(
                candidate.conditions.temperature(),
                &phase_totals,
                &stability,
                &phase_set,
            )?;
            let driving_forces = stability
                .iter()
                .map(|report| report.driving_force)
                .collect::<Vec<_>>();

            match transition_plan {
                PhaseTransitionPlan::Deactivate { phase } => {
                    let phase_index = phase.index();
                    let previous_phase_set = phase_set.clone();
                    let driving_force = stability[phase_index].driving_force.ok_or_else(|| {
                        ReactionExtentError::InvalidCandidate {
                            field: "phase_stability",
                            message: format!(
                                "phase {phase_index} was selected for deactivation without a driving force"
                            ),
                        }
                    })?;
                    let phase_moles = phase_totals[phase_index];
                    crate::Thermodynamics::ChemEquilibrium::equilibrium_workflows::deactivate_phases_seed_only(
                        &mut y,
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
                        phase_totals,
                        driving_forces,
                        reason: PhaseTransitionReason::VanishingUnstableActivePhase {
                            phase_moles,
                            driving_force,
                        },
                        previous_phase_set,
                        new_phase_set,
                        restart_seed: y.clone(),
                        nonlinear_report: candidate.solve_report.clone(),
                        candidate_validation: candidate.validation_report.clone(),
                    });
                    reject_repeated_phase_set(&mut visited, &phase_set, iteration + 1)?;
                    seed = y;
                }
                PhaseTransitionPlan::Activate { phase } => {
                    let previous_phase_set = phase_set.clone();
                    let driving_force =
                        stability[phase.index()].driving_force.ok_or_else(|| {
                            ReactionExtentError::InvalidCandidate {
                                field: "phase_stability",
                                message: format!(
                                    "phase {} was selected for activation without a driving force",
                                    phase.index()
                                ),
                            }
                        })?;
                    seed_activated_phase(
                        &mut y,
                        phase,
                        &species_phase,
                        PhaseSeedPolicy::RelativeToSystemTotal {
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
                        phase_totals,
                        driving_forces,
                        reason: PhaseTransitionReason::UnstableInactivePhase { driving_force },
                        previous_phase_set,
                        new_phase_set,
                        restart_seed: y.clone(),
                        nonlinear_report: candidate.solve_report.clone(),
                        candidate_validation: candidate.validation_report.clone(),
                    });
                    reject_repeated_phase_set(&mut visited, &phase_set, iteration + 1)?;
                    seed = y;
                }
                PhaseTransitionPlan::Hold { phase: _ } => {
                    validate_phase_set_candidate(
                        &phase_totals,
                        &phase_active,
                        self.phase_manager.phase_eps,
                    )?;
                    return self.finish(
                        candidate,
                        phase_set,
                        initial_phase_set,
                        initial_active_phases,
                        transitions,
                        nonlinear_reports,
                        iteration + 1,
                        stability,
                        projection_build,
                        formulation_build,
                        validation_duration,
                        rst_symbolic_reused,
                    );
                }
                PhaseTransitionPlan::NoTransition { .. } => {
                    validate_phase_set_candidate(
                        &phase_totals,
                        &phase_active,
                        self.phase_manager.phase_eps,
                    )?;
                    return self.finish(
                        candidate,
                        phase_set,
                        initial_phase_set,
                        initial_active_phases,
                        transitions,
                        nonlinear_reports,
                        iteration + 1,
                        stability,
                        projection_build,
                        formulation_build,
                        validation_duration,
                        rst_symbolic_reused,
                    );
                }
            }
        }

        Err(ReactionExtentError::PhaseControlDidNotConverge {
            iterations: max_phase_iterations,
        })
    }

    /// Recovers a phase that has no positive interior equilibrium.
    ///
    /// A fixed active set may fail even though the physical equilibrium is
    /// valid on its boundary (for example, all liquid evaporates). We only
    /// retry by removing an active phase when it had positive input inventory;
    /// the alternative active set solves successfully; and its stability
    /// report places the removed phase above the keep threshold. This narrow
    /// contract prevents the recovery path from hiding unrelated numerical or
    /// data errors behind a convenient phase deletion.
    ///
    /// Iterates over active phases that had positive input inventory, removes
    /// each one in turn, and re-solves the reduced active set. Recovery is
    /// accepted only when the reduced solve succeeds and the removed phase's
    /// stability driving force exceeds the keep threshold. Returns `None` when
    /// no phase satisfies all recovery preconditions.
    fn recover_active_boundary<F>(
        &mut self,
        phase_set: &PhaseSet,
        seed: &[f64],
        species_phase: &[usize],
        full_element_totals: &[f64],
        solve_active_set: &mut F,
    ) -> Result<Option<PreparedBoundaryRecovery>, ReactionExtentError>
    where
        F: FnMut(
            &mut Self,
            &[bool],
            &[f64],
            &[usize],
            &[f64],
        ) -> Result<PreparedActiveSetCandidate, ReactionExtentError>,
    {
        let active = phase_set.active_mask();
        if active.iter().filter(|&&is_active| is_active).count() <= 1 {
            return Ok(None);
        }
        for phase_index in 0..active.len() {
            if !active[phase_index] {
                continue;
            }
            let initial_phase_moles = self
                .prepared
                .problem()
                .initial_moles()
                .iter()
                .zip(species_phase.iter())
                .filter(|(_, phase)| **phase == phase_index)
                .map(|(moles, _)| *moles)
                .sum::<f64>();
            if initial_phase_moles <= self.phase_manager.phase_eps {
                continue;
            }

            let phase = PhaseIndex::new(phase_index, active.len())?;
            let mut reduced_phase_set = phase_set.clone();
            reduced_phase_set.deactivate(phase);
            reduced_phase_set.settle_transitions();
            let reduced_active = reduced_phase_set.active_mask();
            let Ok(candidate) = solve_active_set(
                self,
                &reduced_active,
                seed,
                species_phase,
                full_element_totals,
            ) else {
                continue;
            };
            let Ok(stability) = compute_phase_stability_reports(
                &candidate.log_moles,
                &candidate.stability_gibbs,
                self.prepared.problem().phases(),
                species_phase,
                self.prepared.problem().element_composition(),
                candidate.conditions.temperature(),
                candidate.conditions.pressure(),
                candidate.conditions.reference_pressure(),
                &reduced_phase_set,
            ) else {
                continue;
            };
            let phase_total =
                compute_phase_totals(&candidate.log_moles, species_phase)[phase_index];
            let Some(driving_force) = stability[phase_index].driving_force else {
                continue;
            };
            let (_, dg_keep) = self
                .phase_manager
                .thresholds_at(candidate.conditions.temperature())?;
            if phase_total < self.phase_manager.phase_eps && driving_force > dg_keep {
                return Ok(Some(PreparedBoundaryRecovery {
                    phase,
                    candidate,
                    phase_set: reduced_phase_set,
                    stability,
                }));
            }
        }
        Ok(None)
    }

    /// Solves one fixed active set and returns a validated candidate.
    ///
    /// Builds or retrieves the cached active-set projection, prepares the
    /// reduced numerical problem (and optionally the RST symbolic problem),
    /// then dispatches through the solver cascade. The result includes the
    /// solved log-moles, validation report, backend cascade evidence, and
    /// timing for the projection and formulation builds. The active mask is
    /// retained so the outer phase-control loop can associate each candidate
    /// with its phase set.
    fn solve_active_set(
        &mut self,
        active: &[bool],
        full_seed: &[f64],
        species_phase: &[usize],
        full_element_totals: &[f64],
    ) -> Result<PreparedActiveSetCandidate, ReactionExtentError> {
        let cache_key = active.to_vec();
        let projection_cached = self.projection_cache.contains_key(&cache_key);
        let projection_started = (!projection_cached && self.timing_enabled).then(Instant::now);
        if !projection_cached {
            let projection = ActiveSetProjection::build(
                self.prepared.problem().phases(),
                species_phase,
                self.prepared.problem().element_composition(),
                active,
                self.solver_settings.solver_params.tol,
            )?;
            self.projection_cache.insert(cache_key.clone(), projection);
        }
        let projection = self
            .projection_cache
            .get(&cache_key)
            .expect("projection cache entry inserted above");
        let projection_build = projection_started
            .map(|started| started.elapsed())
            .unwrap_or(Duration::ZERO);
        projection.validate_element_totals_representable(
            full_element_totals,
            self.solver_settings.solver_params.tol,
        )?;
        let reduced_seed = projection.project_log_moles(full_seed)?;
        let reduced_components = projection
            .active_species
            .iter()
            .map(|species| self.prepared.problem().components()[species.index()].clone())
            .collect();
        let reduced_initial_moles = reduced_seed.iter().map(|value| value.exp()).collect();
        let reduced_gibbs: Vec<_> = projection
            .active_species
            .iter()
            .map(|species| self.prepared.problem().gibbs()[species.index()].clone())
            .collect();
        let reduced_symbols = if self.symbolic_standard_gibbs.is_empty() {
            Vec::new()
        } else {
            projection
                .active_species
                .iter()
                .map(|species| self.symbolic_standard_gibbs[species.index()].clone())
                .collect()
        };
        let reduced_totals = projection.reduced_element_totals(full_element_totals)?;
        let reduced_seed = LogMolesInitialGuess::new(reduced_seed.clone())?;
        let current_conditions = self.prepared.problem().conditions();
        let solver_settings = self.solver_settings.clone();
        let use_rst = self.uses_rst_backend(&reduced_symbols);
        let mut rst_symbolic_reused = false;
        let mut formulation_build = Duration::ZERO;

        if !self.prepared_active_set_cache.contains_key(&cache_key) {
            let formulation_started = self.timing_enabled.then(Instant::now);
            let reduced_problem = EquilibriumProblem::new(
                reduced_components,
                reduced_initial_moles,
                reduced_seed.clone(),
                projection.element_composition.clone(),
                reduced_gibbs.clone(),
                projection.phases.clone(),
                current_conditions,
            )?;
            let prepared = PreparedEquilibriumProblem::new_with_element_totals(
                reduced_problem,
                Some(reduced_totals),
            )?;
            let entry_build = formulation_started
                .map(|started| started.elapsed())
                .unwrap_or(Duration::ZERO);
            self.prepared_active_set_cache.insert(
                cache_key.clone(),
                PreparedActiveSetCacheEntry {
                    prepared,
                    symbolic_standard_gibbs: Vec::new(),
                    rst_problem: None,
                    formulation_build: entry_build,
                },
            );
            formulation_build += entry_build;
        }

        // The cache entry owns both equation preparations. Retargeting makes
        // a cheap immutable numeric copy for this point, while the RST object
        // itself remains borrowed and receives the new `T` plus `G0_i` values.
        let outcome = {
            let entry = self
                .prepared_active_set_cache
                .get_mut(&cache_key)
                .expect("active-set formulation inserted above");
            let prepared = entry.prepared.retarget_with_gibbs(
                current_conditions,
                reduced_seed.clone(),
                reduced_gibbs,
            )?;
            if use_rst {
                if entry.rst_problem.is_none() {
                    let formulation_started = self.timing_enabled.then(Instant::now);
                    entry.rst_problem = Some(prepare_rst_symbolic_problem_from_prepared(
                        &prepared,
                        &reduced_symbols,
                    )?);
                    entry.symbolic_standard_gibbs = reduced_symbols.clone();
                    let symbolic_build = formulation_started
                        .map(|started| started.elapsed())
                        .unwrap_or(Duration::ZERO);
                    entry.formulation_build += symbolic_build;
                    formulation_build += symbolic_build;
                } else {
                    rst_symbolic_reused = true;
                }
                entry
                    .rst_problem
                    .as_mut()
                    .expect("RST problem was created or retained above")
                    .set_fixed_pt_thermochemistry_from_prepared(&prepared)?;
            } else {
                // An explicitly legacy-only policy must not retain an unused
                // symbolic backend object or report it as reusable.
                entry.rst_problem = None;
                entry.symbolic_standard_gibbs.clear();
            }

            let mut runner = PreparedEquilibriumRunner::new(prepared, reduced_symbols.clone())?;
            runner.configure().clone_from(&solver_settings);
            runner.configure().keq_validation_mode = EquilibriumConstantValidationMode::Off;

            if use_rst {
                let rst_problem = entry
                    .rst_problem
                    .as_ref()
                    .expect("RST problem prepared above");
                runner.solve_from_seed_with_rst(reduced_seed, rst_problem)?
            } else {
                runner.solve_from_seed(reduced_seed)?
            }
        };
        let log_moles = projection.scatter_log_moles(
            outcome.solution.log_moles(),
            PHASE_CONTROL_TRACE_MOLE_FLOOR.ln(),
        )?;
        Ok(PreparedActiveSetCandidate {
            log_moles,
            solved_active_mask: active.to_vec(),
            validation_report: outcome.solution.validation().clone(),
            solve_report: outcome.solve_report,
            keq_validation_status: outcome.keq_validation_status,
            stability_gibbs: self.prepared.problem().gibbs().to_vec(),
            conditions: self.prepared.problem().conditions(),
            projection_build,
            formulation_build,
            validation_duration: outcome.validation_duration,
            rst_symbolic_reused,
        })
    }

    /// Assembles the final [`PreparedPhaseControlOutcome`] from the accepted
    /// active-set candidate and the outer-loop transition history.
    ///
    /// Computes phase stability reports for the accepted phase set, builds the
    /// multiphase acceptance report, and packages everything into a single
    /// immutable outcome value. The initial and final phase sets are retained
    /// so the caller can inspect the complete transition sequence without
    /// re-running the outer loop.
    fn finish(
        &self,
        candidate: PreparedActiveSetCandidate,
        phase_set: PhaseSet,
        initial_phase_set: PhaseSet,
        initial_active_phases: Vec<PhaseIndex>,
        transitions: Vec<PhaseTransitionRecord>,
        nonlinear_reports: Vec<EquilibriumSolveReport>,
        iterations: usize,
        phase_stability: Vec<
            crate::Thermodynamics::ChemEquilibrium::equilibrium_workflows::PhaseStabilityReport,
        >,
        projection_build: Duration,
        formulation_build: Duration,
        validation_duration: Duration,
        rst_symbolic_reused: bool,
    ) -> Result<PreparedPhaseControlOutcome, ReactionExtentError> {
        let phase_control_report = PhaseControlledSolveReport {
            iterations,
            initial_active_phases,
            final_active_phases: phase_set.active_phases()?,
            initial_phase_set,
            final_phase_set: phase_set.clone(),
            transitions,
            final_validation: candidate.validation_report.clone(),
            nonlinear_reports,
        };
        let (dg_create, dg_keep) = self
            .phase_manager
            .thresholds_at(candidate.conditions.temperature())?;
        let acceptance_report = build_multiphase_acceptance_report(
            phase_control_report.clone(),
            phase_stability,
            Some(dg_create),
            Some(dg_keep),
        )?;
        let final_mask = phase_set.active_mask();
        let phase_statuses = (0..final_mask.len())
            .map(|index| PhaseIndex::new(index, final_mask.len()))
            .collect::<Result<Vec<_>, _>>()?
            .into_iter()
            .map(|phase| phase_set.status(phase))
            .collect::<Vec<_>>();
        let solve_report = candidate.solve_report.clone();
        let log_moles = candidate.log_moles;
        let moles = log_moles.iter().map(|value| value.exp()).collect();
        let solution =
            crate::Thermodynamics::ChemEquilibrium::equilibrium_problem::EquilibriumSolution::new(
                log_moles,
                moles,
                candidate.conditions,
                candidate.validation_report,
            )?;
        Ok(PreparedPhaseControlOutcome {
            solution,
            solve_report,
            keq_validation_status: candidate.keq_validation_status,
            phase_control_report,
            acceptance_report,
            phase_statuses,
            projection_build,
            formulation_build,
            validation_duration,
            rst_symbolic_reused,
        })
    }
}
