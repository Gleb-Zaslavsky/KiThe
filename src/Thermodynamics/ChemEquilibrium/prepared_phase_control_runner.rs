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
use crate::Thermodynamics::ChemEquilibrium::equilibrium_log_moles::EquilibriumSolverSettings;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_nonlinear::ReactionExtentError;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_prepared_runner::PreparedEquilibriumRunner;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_problem::{
    EquilibriumProblem, LogMolesInitialGuess, PreparedEquilibriumProblem,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_rst_backend::{
    RstPreparedProblem, prepare_rst_symbolic_problem_from_prepared,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_solver_policy::{
    EquilibriumSolveReport, SolverBackend,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_validation::EquilibriumCandidateReport;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_workflows::{
    MultiphaseAcceptanceReport, PHASE_CONTROL_TRACE_MOLE_FLOOR, PhaseControlledSolveReport,
    PhaseManager, PhaseSeedPolicy, PhaseSet, PhaseStatus, PhaseTransitionPlan,
    PhaseTransitionReason, PhaseTransitionRecord, build_multiphase_acceptance_report,
    compute_phase_stability_reports, compute_phase_totals, initial_phase_activity_from_moles,
    reject_repeated_phase_set, seed_activated_phase, validate_phase_set_candidate,
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
    pub(crate) validation_duration: Duration,
    pub(crate) rst_symbolic_reused: bool,
}

#[derive(Debug)]
struct PreparedActiveSetCandidate {
    log_moles: Vec<f64>,
    validation_report: EquilibriumCandidateReport,
    solve_report: EquilibriumSolveReport,
    keq_validation_status: Option<EquilibriumConstantCrossValidationStatus>,
    projection_build: Duration,
    validation_duration: Duration,
    rst_symbolic_reused: bool,
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
/// The symbolic snapshot is compared before reuse. This matters because a
/// temperature interval may cross a thermochemical coefficient boundary: the
/// active mask can remain unchanged while the symbolic equation itself must
/// still be rebuilt.
struct PreparedActiveSetCacheEntry {
    prepared: PreparedEquilibriumProblem,
    symbolic_standard_gibbs: Vec<Expr>,
    rst_problem: Option<RstPreparedProblem>,
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

    /// Retargets the immutable prepared formulation for the next temperature.
    ///
    /// The elemental inventory and reaction basis are intentionally retained.
    /// Only temperature-dependent conditions, Gibbs closures, and the seed are
    /// replaced. A changed active mask will populate one additional projection
    /// cache entry when it is first solved.
    pub(crate) fn retarget(
        &mut self,
        conditions: crate::Thermodynamics::ChemEquilibrium::equilibrium_problem::EquilibriumConditions,
        seed: LogMolesInitialGuess,
        gibbs: Vec<crate::Thermodynamics::ChemEquilibrium::equilibrium_log_moles::GibbsFn>,
        symbolic_standard_gibbs: Vec<Expr>,
    ) -> Result<(), ReactionExtentError> {
        self.prepared = self
            .prepared
            .retarget_with_gibbs(conditions, seed.clone(), gibbs)?;
        if !symbolic_standard_gibbs.is_empty()
            && symbolic_standard_gibbs.len() != self.prepared.problem().species().len()
        {
            return Err(ReactionExtentError::DimensionMismatch(format!(
                "symbolic Gibbs snapshot has {} entries for {} species",
                symbolic_standard_gibbs.len(),
                self.prepared.problem().species().len()
            )));
        }
        self.symbolic_standard_gibbs = symbolic_standard_gibbs;
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

    /// Executes the bounded outer loop and publishes only the final accepted
    /// immutable solution plus complete transition evidence.
    pub(crate) fn solve(&mut self) -> Result<PreparedPhaseControlOutcome, ReactionExtentError> {
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
        let mut validation_duration = Duration::ZERO;
        let mut rst_symbolic_reused = false;
        let full_element_totals = self.prepared.element_totals().to_vec();

        for iteration in 0..max_phase_iterations {
            phase_set.settle_transitions();
            let phase_active = phase_set.active_mask();
            let mut candidate = match self.solve_active_set(
                &phase_active,
                &seed,
                &species_phase,
                &full_element_totals,
            ) {
                Ok(candidate) => {
                    projection_build += candidate.projection_build;
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
            let stability = compute_phase_stability_reports(
                &y,
                self.prepared.problem().gibbs(),
                self.prepared.problem().phases(),
                &species_phase,
                self.prepared.problem().element_composition(),
                self.prepared.problem().conditions().temperature(),
                self.prepared.problem().conditions().pressure(),
                self.prepared.problem().conditions().reference_pressure(),
                &phase_set,
            )?;
            let transition_plan = self.phase_manager.classify_phases_at_temperature(
                self.prepared.problem().conditions().temperature(),
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
    fn recover_active_boundary(
        &mut self,
        phase_set: &PhaseSet,
        seed: &[f64],
        species_phase: &[usize],
        full_element_totals: &[f64],
    ) -> Result<Option<PreparedBoundaryRecovery>, ReactionExtentError> {
        let active = phase_set.active_mask();
        if active.iter().filter(|&&is_active| is_active).count() <= 1 {
            return Ok(None);
        }
        let (_, dg_keep) = self
            .phase_manager
            .thresholds_at(self.prepared.problem().conditions().temperature())?;

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
            let Ok(candidate) =
                self.solve_active_set(&reduced_active, seed, species_phase, full_element_totals)
            else {
                continue;
            };
            let Ok(stability) = compute_phase_stability_reports(
                &candidate.log_moles,
                self.prepared.problem().gibbs(),
                self.prepared.problem().phases(),
                species_phase,
                self.prepared.problem().element_composition(),
                self.prepared.problem().conditions().temperature(),
                self.prepared.problem().conditions().pressure(),
                self.prepared.problem().conditions().reference_pressure(),
                &reduced_phase_set,
            ) else {
                continue;
            };
            let phase_total =
                compute_phase_totals(&candidate.log_moles, species_phase)[phase_index];
            let Some(driving_force) = stability[phase_index].driving_force else {
                continue;
            };
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
        let current_temperature = current_conditions.temperature();
        let solver_settings = self.solver_settings.clone();
        let use_rst = self.uses_rst_backend(&reduced_symbols);
        let mut rst_symbolic_reused = false;

        if !self.prepared_active_set_cache.contains_key(&cache_key) {
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
            self.prepared_active_set_cache.insert(
                cache_key.clone(),
                PreparedActiveSetCacheEntry {
                    prepared,
                    symbolic_standard_gibbs: Vec::new(),
                    rst_problem: None,
                },
            );
        }

        // The cache entry owns both equation preparations. Retargeting makes
        // a cheap immutable numeric copy for this point, while the RST object
        // itself remains borrowed for the solve and receives only a new `T`.
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
                let symbolic_changed =
                    entry.rst_problem.is_none() || entry.symbolic_standard_gibbs != reduced_symbols;
                rst_symbolic_reused = !symbolic_changed;
                if symbolic_changed {
                    entry.rst_problem = Some(prepare_rst_symbolic_problem_from_prepared(
                        &prepared,
                        &reduced_symbols,
                    )?);
                    entry.symbolic_standard_gibbs = reduced_symbols.clone();
                } else {
                    entry
                        .rst_problem
                        .as_mut()
                        .expect("cached RST problem exists when symbolic snapshot is unchanged")
                        .set_temperature(current_temperature)?;
                }
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
            validation_report: outcome.solution.validation().clone(),
            solve_report: outcome.solve_report,
            keq_validation_status: outcome.keq_validation_status,
            projection_build,
            validation_duration: outcome.validation_duration,
            rst_symbolic_reused,
        })
    }

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
            .thresholds_at(self.prepared.problem().conditions().temperature())?;
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
                self.prepared.problem().conditions(),
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
            validation_duration,
            rst_symbolic_reused,
        })
    }
}
