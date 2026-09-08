//! Immutable active-set orchestration for bounded fixed-`P,T` equilibrium.
//!
//! The numerical inner problem is prepared once per settled active phase set,
//! while the outer loop owns only seeds, phase transitions, and reports. Each
//! active-set entry also retains its RST symbolic problem, so a temperature
//! sweep updates the shared `T` parameter instead of rebuilding equations and
//! lambdified closures. The runner never constructs the historical
//! `EquilibriumLogMoles` god object. That object remains available for
//! compatibility callers; this runner is the canonical bridge path.

#[cfg(test)]
use std::cell::Cell;
use std::collections::{HashMap, HashSet};
use std::time::{Duration, Instant};

use crate::Thermodynamics::ChemEquilibrium::equilibrium_active_set::ActiveSetProjection;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_constant_cross_validation::EquilibriumConstantCrossValidationStatus;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_constant_validation::EquilibriumConstantValidationMode;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_diagnostics::{
    EquilibriumDiagnosticEvent, EquilibriumDiagnosticsCollector, EquilibriumDiagnosticsMode,
    EquilibriumDiagnosticsOptions, EquilibriumDiagnosticsReport, PhaseStabilityDiagnostic,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_ids::PhaseIndex;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_log_moles::{
    EquilibriumSolverSettings, GibbsFn,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_nonlinear::ReactionExtentError;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_phase_stability::{
    PhaseStabilityGeometryCache, PhaseStabilityGeometryCacheStats, PhaseStabilityReport,
};
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
    PhaseManager, PhaseSet, PhaseStatus, PhaseTotalSeedPolicy, PhaseTransitionPlan,
    PhaseTransitionReason, PhaseTransitionRecord, build_multiphase_acceptance_report,
    compute_phase_stability_reports_with_geometry_cache, compute_phase_totals,
    initial_phase_activity_from_moles, reject_repeated_phase_set,
    seed_activated_phase_with_composition, validate_phase_set_candidate,
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
    pub(crate) diagnostics: EquilibriumDiagnosticsReport,
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

/// State consumed by one phase-control attempt that must survive an error.
///
/// Prepared projections and symbolic graphs are internal structural
/// memoization: they never publish a solution and are always retargeted before
/// use. Continuation, however, is a physical seed selected from an earlier
/// accepted point. Consuming it on an unsuccessful attempt would propagate a
/// failed lifecycle transition into the next request.
#[derive(Debug, Clone)]
struct ContinuationCheckpoint {
    seed: Option<Vec<f64>>,
    phase_set: Option<PhaseSet>,
}

/// A failure injected only by unit tests after a local phase transition has
/// been assembled, but before its restart seed can be committed.
///
/// This is deliberately not part of the production policy.  It verifies the
/// transactional boundary where a phase-control attempt may have rich
/// diagnostic evidence while its mutable continuation state must still be
/// discarded on failure.
#[cfg(test)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum TestFailPoint {
    AfterPhaseTransitionBeforeCommit,
}

// Test-only injection state for end-to-end prepared-workflow tests. A
// thread-local counter keeps parallel test threads isolated and avoids adding
// a fault-injection knob to the production solve options.
#[cfg(test)]
thread_local! {
    static TEST_FAILPOINT_TRANSITIONS_REMAINING: Cell<Option<usize>> = const { Cell::new(None) };
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
    /// Reuses only active elemental geometry for TPD; all state-dependent
    /// chemical potentials and minimizers remain recomputed per candidate.
    phase_stability_geometry_cache: PhaseStabilityGeometryCache,
    diagnostics: EquilibriumDiagnosticsCollector,
    #[cfg(test)]
    test_failpoint: Option<TestFailPoint>,
}

impl PreparedPhaseControlRunner {
    /// Test-only convenience constructor with diagnostics disabled.
    #[cfg(test)]
    pub(crate) fn new(
        problem: EquilibriumProblem,
        symbolic_standard_gibbs: Vec<Expr>,
        timing_enabled: bool,
    ) -> Result<Self, ReactionExtentError> {
        Self::new_with_diagnostics(
            problem,
            symbolic_standard_gibbs,
            timing_enabled,
            EquilibriumDiagnosticsOptions::disabled(),
        )
    }

    /// Creates the runner with one explicitly configured observational trace.
    ///
    /// This keeps diagnostics at the orchestration boundary; numerical
    /// residual and TPD primitives remain pure functions.
    pub(crate) fn new_with_diagnostics(
        problem: EquilibriumProblem,
        symbolic_standard_gibbs: Vec<Expr>,
        timing_enabled: bool,
        diagnostics: EquilibriumDiagnosticsOptions,
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
            phase_stability_geometry_cache: PhaseStabilityGeometryCache::default(),
            diagnostics: EquilibriumDiagnosticsCollector::new(diagnostics),
            #[cfg(test)]
            test_failpoint: None,
        })
    }

    /// Arms a one-shot failure at the transition/restart commit boundary.
    ///
    /// Tests use this to exercise rollback after a meaningful transition has
    /// already been diagnosed.  The hook is compiled out of production
    /// builds, so it cannot alter normal solver behavior or API semantics.
    #[cfg(test)]
    pub(crate) fn arm_test_failpoint_once(&mut self, failpoint: TestFailPoint) {
        self.test_failpoint = Some(failpoint);
    }

    /// Arms the same boundary failure for the next prepared runner created on
    /// this test thread after `completed_transitions` transitions. This is the
    /// narrow seam used by real-data range stories; it remains test-only and
    /// does not become part of the public solver configuration.
    #[cfg(test)]
    pub(crate) fn arm_test_failpoint_after_transitions(completed_transitions: usize) {
        TEST_FAILPOINT_TRANSITIONS_REMAINING.with(|remaining| {
            remaining.set(Some(completed_transitions));
        });
    }

    #[cfg(test)]
    fn trip_test_failpoint(&mut self) -> Result<(), ReactionExtentError> {
        if self.test_failpoint == Some(TestFailPoint::AfterPhaseTransitionBeforeCommit) {
            self.test_failpoint = None;
            return Err(ReactionExtentError::InvalidProblem {
                field: "test_failpoint",
                message: "test-only failure after phase transition before commit".to_string(),
            });
        }
        Ok(())
    }

    #[cfg(test)]
    fn trip_global_test_failpoint(&mut self) -> Result<(), ReactionExtentError> {
        let should_trip = TEST_FAILPOINT_TRANSITIONS_REMAINING.with(|remaining| {
            let Some(count) = remaining.get() else {
                return false;
            };
            if count == 0 {
                remaining.set(None);
                true
            } else {
                remaining.set(Some(count - 1));
                false
            }
        });
        if should_trip {
            let backend = crate::Thermodynamics::ChemEquilibrium::equilibrium_solver_policy::
                SolverBackend::Legacy(
                    crate::Thermodynamics::ChemEquilibrium::equilibrium_log_moles::Solvers::NR,
                );
            return Err(ReactionExtentError::AllBackendsFailed {
                attempts: vec![
                    crate::Thermodynamics::ChemEquilibrium::equilibrium_solver_policy::
                        SolverAttemptReport {
                        backend,
                        outcome: crate::Thermodynamics::ChemEquilibrium::equilibrium_solver_policy::
                            SolverAttemptOutcome::Failed {
                            kind: crate::Thermodynamics::ChemEquilibrium::equilibrium_solver_policy::
                                SolverAttemptFailureKind::Solver,
                            reason: "test-only failure after phase transition before commit"
                                .to_string(),
                        },
                        metrics: None,
                    },
                ],
            });
        }
        Ok(())
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

    /// Installs a numerical basin seed without carrying an accepted phase set.
    ///
    /// This is intentionally distinct from continuation. Phase activity is
    /// still derived from the current physical inventory and configured
    /// `InitialPhaseSet`, so activation/deactivation evidence is reproduced by
    /// the current transaction instead of inherited from the seed source.
    pub(crate) fn set_independent_basin_seed(
        &mut self,
        seed: LogMolesInitialGuess,
    ) -> Result<(), ReactionExtentError> {
        if seed.as_slice().len() != self.prepared.problem().species().len() {
            return Err(ReactionExtentError::DimensionMismatch(format!(
                "independent basin seed has {} entries for {} prepared species",
                seed.as_slice().len(),
                self.prepared.problem().species().len()
            )));
        }
        self.continuation_seed = Some(seed.into_inner());
        self.continuation_phase_set = None;
        Ok(())
    }

    /// Clears any carried physical continuation before an independent first
    /// range point. A numerical seed alone is not continuation: only a seed
    /// paired with a previously accepted phase set may influence lifecycle
    /// initialization.
    pub(crate) fn clear_continuation_state(&mut self) {
        self.continuation_seed = None;
        self.continuation_phase_set = None;
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

    /// Rebuild/reuse evidence for cached SVD/range/null-space geometries used
    /// by the canonical TPD phase-stability service.
    pub(crate) fn phase_stability_geometry_cache_statistics(
        &self,
    ) -> PhaseStabilityGeometryCacheStats {
        self.phase_stability_geometry_cache.statistics()
    }

    /// Replaces diagnostics for the next solve transaction. Range templates
    /// call this only after the preceding point has been published.
    pub(crate) fn set_diagnostics_options(&mut self, options: EquilibriumDiagnosticsOptions) {
        self.diagnostics.set_options(options);
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
        let checkpoint = ContinuationCheckpoint {
            seed: self.continuation_seed.clone(),
            phase_set: self.continuation_phase_set.clone(),
        };
        let result = self.solve_with_fixed_active_solver_inner(&mut solve_active_set);
        if let Err(error) = &result {
            // Only an accepted result may advance continuation. Structural
            // cache entries are not solution publication and remain valid
            // memoization for the same immutable layout.
            self.continuation_seed = checkpoint.seed;
            self.continuation_phase_set = checkpoint.phase_set;
            self.diagnostics.record(
                EquilibriumDiagnosticsMode::PhaseLifecycle,
                EquilibriumDiagnosticEvent::ContinuationRestored {
                    retained_phase_set: self.continuation_phase_set.clone(),
                    retained_seed: self.continuation_seed.is_some(),
                },
            );
            self.diagnostics.record(
                EquilibriumDiagnosticsMode::Summary,
                EquilibriumDiagnosticEvent::SolveFailed {
                    message: error.to_string(),
                    continuation_restored: true,
                },
            );
        }
        result
    }

    fn solve_with_fixed_active_solver_inner<F>(
        &mut self,
        mut solve_active_set: &mut F,
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
        self.diagnostics.record(
            EquilibriumDiagnosticsMode::Summary,
            EquilibriumDiagnosticEvent::SolveStarted {
                conditions: self.prepared.problem().conditions(),
                initial_phase_set: initial_phase_set.clone(),
            },
        );
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
            self.diagnostics.record(
                EquilibriumDiagnosticsMode::PhaseLifecycle,
                EquilibriumDiagnosticEvent::OuterIterationStarted {
                    iteration,
                    active_phase_set: phase_set.clone(),
                },
            );
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
                    self.diagnostics.record(
                        EquilibriumDiagnosticsMode::Detailed,
                        EquilibriumDiagnosticEvent::ActiveSetCandidateRejected {
                            iteration,
                            active_phase_set: phase_set.clone(),
                            message: primary_error.to_string(),
                        },
                    );
                    let Some(recovery) = self.recover_active_boundary(
                        iteration,
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
                    let minimum_tpds = stability
                        .iter()
                        .map(|report| report.minimum_tpd)
                        .collect::<Vec<_>>();
                    let minimum_tpd = stability[phase.index()].minimum_tpd.ok_or_else(|| {
                        ReactionExtentError::InvalidCandidate {
                            field: "phase_stability",
                            message: format!(
                                "boundary recovery selected phase {} without a TPD minimum",
                                phase.index()
                            ),
                        }
                    })?;
                    // Report the actual phase inventory immediately before
                    // the recovery probe. Construction-time inventory can be
                    // stale after accepted continuation across a T range.
                    let initial_phase_moles = seed
                        .iter()
                        .zip(species_phase.iter())
                        .filter(|(_, phase_index)| **phase_index == phase.index())
                        .map(|(log_moles, _)| log_moles.exp())
                        .sum::<f64>();
                    let previous_phase_set = phase_set.clone();
                    phase_set = recovered_phase_set;
                    let (dg_create, dg_keep) = self
                        .phase_manager
                        .thresholds_at(candidate.conditions.temperature())?;
                    self.record_stability_diagnostics(
                        iteration,
                        dg_create,
                        dg_keep,
                        &phase_totals,
                        &stability,
                    );
                    transitions.push(PhaseTransitionRecord {
                        iteration,
                        transition_duration: transition_started.elapsed(),
                        activated: Vec::new(),
                        deactivated: vec![phase],
                        phase_totals,
                        minimum_tpds,
                        incipient_composition: None,
                        reason: PhaseTransitionReason::BoundaryUnstableActivePhase {
                            initial_phase_moles,
                            minimum_tpd,
                        },
                        previous_phase_set: previous_phase_set.clone(),
                        new_phase_set: phase_set.clone(),
                        restart_seed: candidate.log_moles.clone(),
                        nonlinear_report: candidate.solve_report.clone(),
                        candidate_validation: candidate.validation_report.clone(),
                    });
                    self.diagnostics.record(
                        EquilibriumDiagnosticsMode::Detailed,
                        EquilibriumDiagnosticEvent::RecoveryProbeAccepted {
                            iteration,
                            phase_index: phase.index(),
                            previous_phase_set,
                            new_phase_set: phase_set.clone(),
                        },
                    );
                    #[cfg(test)]
                    self.trip_test_failpoint()?;
                    #[cfg(test)]
                    self.trip_global_test_failpoint()?;
                    self.reject_repeated_phase_set(&mut visited, &phase_set, iteration + 1)?;
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
            self.diagnostics.record(
                EquilibriumDiagnosticsMode::Detailed,
                EquilibriumDiagnosticEvent::ActiveSetCandidateAccepted {
                    iteration,
                    active_phase_set: phase_set.clone(),
                    validation: candidate.validation_report.clone(),
                    backend_summary: candidate.solve_report.summary(),
                },
            );
            let probe_stability = compute_phase_stability_reports_with_geometry_cache(
                &y,
                &candidate.stability_gibbs,
                self.prepared.problem().phases(),
                &species_phase,
                self.prepared.problem().element_composition(),
                candidate.conditions.temperature(),
                candidate.conditions.pressure(),
                candidate.conditions.reference_pressure(),
                &phase_set,
                &mut self.phase_stability_geometry_cache,
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
            // An all-active monolithic P,H result is a bounded numerical probe,
            // not an activation decision. Keep only a phase whose canonical
            // TPD is below the creation threshold, restore every other probed
            // phase to trace, and restart from the TPD minimizer composition.
            let mut probe_expanded = false;
            let (dg_create, _) = self
                .phase_manager
                .thresholds_at(candidate.conditions.temperature())?;
            let (_, dg_keep) = self
                .phase_manager
                .thresholds_at(candidate.conditions.temperature())?;
            self.record_stability_diagnostics(
                iteration,
                dg_create,
                dg_keep,
                &phase_totals,
                &probe_stability,
            );
            let initially_inactive = phase_active
                .iter()
                .enumerate()
                .filter(|(_, is_active)| !**is_active)
                .map(|(phase_index, _)| PhaseIndex::new(phase_index, phase_active.len()))
                .collect::<Result<Vec<_>, _>>()?;
            let mut tpd_activation = None;
            for phase_index in 0..phase_active.len() {
                if phase_active[phase_index] || !candidate.solved_active_mask[phase_index] {
                    continue;
                }
                let phase = PhaseIndex::new(phase_index, phase_active.len())?;
                let minimum_tpd = probe_stability[phase_index].minimum_tpd.ok_or_else(|| {
                    ReactionExtentError::InvalidCandidate {
                        field: "phase_stability",
                        message: format!(
                            "monolithic active-set probe expanded phase {} without a TPD minimum",
                            phase.index()
                        ),
                    }
                })?;
                if minimum_tpd >= dg_create {
                    continue;
                }
                let composition = probe_stability[phase_index]
                    .incipient_composition
                    .as_deref()
                    .ok_or_else(|| ReactionExtentError::InvalidCandidate {
                        field: "phase_stability",
                        message: format!(
                            "monolithic active-set probe expanded phase {} without a TPD minimizer composition",
                            phase.index()
                        ),
                    })?;
                // One lifecycle transition per outer pass keeps the phase-set
                // fingerprint, hysteresis, and rollback semantics deterministic.
                tpd_activation = Some((phase, minimum_tpd, composition.to_vec()));
                break;
            }
            if !initially_inactive.is_empty() {
                crate::Thermodynamics::ChemEquilibrium::equilibrium_workflows::deactivate_phases_seed_only(
                    &mut y,
                    &initially_inactive,
                    &species_phase,
                    PHASE_CONTROL_TRACE_MOLE_FLOOR,
                )?;
            }
            if let Some((phase, minimum_tpd, composition)) = tpd_activation {
                let previous_phase_set = phase_set.clone();
                seed_activated_phase_with_composition(
                    &mut y,
                    phase,
                    &species_phase,
                    &composition,
                    PhaseTotalSeedPolicy::RelativeToSystemTotal {
                        fraction: 1e-8,
                        minimum: PHASE_CONTROL_TRACE_MOLE_FLOOR,
                    },
                )?;
                phase_set.activate(phase);
                phase_active = phase_set.active_mask();
                probe_expanded = true;
                transitions.push(PhaseTransitionRecord {
                    iteration,
                    transition_duration: transition_started.elapsed(),
                    activated: vec![phase],
                    deactivated: Vec::new(),
                    phase_totals: compute_phase_totals(&y, &species_phase),
                    minimum_tpds: probe_stability
                        .iter()
                        .map(|report| report.minimum_tpd)
                        .collect(),
                    incipient_composition: Some(composition.clone()),
                    reason: PhaseTransitionReason::UnstableInactivePhase { minimum_tpd },
                    previous_phase_set: previous_phase_set.clone(),
                    new_phase_set: phase_set.clone(),
                    restart_seed: y.clone(),
                    nonlinear_report: candidate.solve_report.clone(),
                    candidate_validation: candidate.validation_report.clone(),
                });
                self.record_transition_diagnostics(
                    iteration,
                    &[phase],
                    &[],
                    PhaseTransitionReason::UnstableInactivePhase { minimum_tpd },
                    previous_phase_set,
                    phase_set.clone(),
                    Some(composition),
                );
                #[cfg(test)]
                self.trip_test_failpoint()?;
                #[cfg(test)]
                self.trip_global_test_failpoint()?;
                self.reject_repeated_phase_set(&mut visited, &phase_set, iteration + 1)?;
            } else if candidate.solved_active_mask != phase_active {
                // The probe located a numerical branch but supplied no TPD
                // evidence for publishing a wider physical phase set. It is
                // deliberately rejected: a caller may select the independent
                // nested P,H route, but monolithic phase control must never
                // turn neutral probe occupancy into an accepted phase.
                return Err(ReactionExtentError::PhMonolithicPhaseProbeRejected {
                    message: format!(
                        "the all-active recovery probe found no inactive phase below the TPD creation threshold {dg_create:e}; the probe cannot be published as an accepted phase set"
                    ),
                });
            }
            if probe_expanded {
                // The probe only selected the physical restart state. Solve
                // the newly activated fixed set before it can be accepted or
                // contribute final stability evidence.
                seed = y;
                continue;
            }
            let stability = probe_stability;
            let transition_plan = self.phase_manager.classify_phases_at_temperature(
                candidate.conditions.temperature(),
                &phase_totals,
                &stability,
                &phase_set,
            )?;
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
                        minimum_tpds,
                        incipient_composition: None,
                        reason: PhaseTransitionReason::VanishingUnstableActivePhase {
                            phase_moles,
                            minimum_tpd,
                        },
                        previous_phase_set: previous_phase_set.clone(),
                        new_phase_set: new_phase_set.clone(),
                        restart_seed: y.clone(),
                        nonlinear_report: candidate.solve_report.clone(),
                        candidate_validation: candidate.validation_report.clone(),
                    });
                    self.record_transition_diagnostics(
                        iteration,
                        &[],
                        &[phase],
                        PhaseTransitionReason::VanishingUnstableActivePhase {
                            phase_moles,
                            minimum_tpd,
                        },
                        previous_phase_set.clone(),
                        new_phase_set.clone(),
                        None,
                    );
                    #[cfg(test)]
                    self.trip_test_failpoint()?;
                    #[cfg(test)]
                    self.trip_global_test_failpoint()?;
                    self.reject_repeated_phase_set(&mut visited, &phase_set, iteration + 1)?;
                    seed = y;
                }
                PhaseTransitionPlan::Activate { phase } => {
                    let previous_phase_set = phase_set.clone();
                    let minimum_tpd = stability[phase.index()].minimum_tpd.ok_or_else(|| {
                        ReactionExtentError::InvalidCandidate {
                            field: "phase_stability",
                            message: format!(
                                "phase {} was selected for activation without a TPD minimum",
                                phase.index()
                            ),
                        }
                    })?;
                    let composition = stability[phase.index()]
                        .incipient_composition
                        .as_deref()
                        .ok_or_else(|| ReactionExtentError::InvalidCandidate {
                            field: "phase_stability",
                            message: format!(
                                "phase {} was selected for activation without a TPD minimizer composition",
                                phase.index()
                            ),
                        })?;
                    seed_activated_phase_with_composition(
                        &mut y,
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
                        // Activation records describe the actual restart
                        // state.  The TPD seed has already been installed in
                        // `y`, so retaining the pre-seed trace total here
                        // would make the report contradict `restart_seed`.
                        phase_totals: compute_phase_totals(&y, &species_phase),
                        minimum_tpds,
                        incipient_composition: Some(composition.to_vec()),
                        reason: PhaseTransitionReason::UnstableInactivePhase { minimum_tpd },
                        previous_phase_set: previous_phase_set.clone(),
                        new_phase_set: new_phase_set.clone(),
                        restart_seed: y.clone(),
                        nonlinear_report: candidate.solve_report.clone(),
                        candidate_validation: candidate.validation_report.clone(),
                    });
                    self.record_transition_diagnostics(
                        iteration,
                        &[phase],
                        &[],
                        PhaseTransitionReason::UnstableInactivePhase { minimum_tpd },
                        previous_phase_set,
                        new_phase_set,
                        Some(composition.to_vec()),
                    );
                    #[cfg(test)]
                    self.trip_test_failpoint()?;
                    #[cfg(test)]
                    self.trip_global_test_failpoint()?;
                    self.reject_repeated_phase_set(&mut visited, &phase_set, iteration + 1)?;
                    seed = y;
                }
                PhaseTransitionPlan::Hold { phase } => {
                    self.diagnostics.record(
                        EquilibriumDiagnosticsMode::PhaseLifecycle,
                        EquilibriumDiagnosticEvent::TransitionHeldByHysteresis {
                            iteration,
                            phase_index: phase.index(),
                        },
                    );
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

        self.diagnostics.record(
            EquilibriumDiagnosticsMode::PhaseLifecycle,
            EquilibriumDiagnosticEvent::PhaseControlBudgetExhausted {
                max_outer_iterations: max_phase_iterations,
            },
        );
        Err(ReactionExtentError::PhaseControlDidNotConverge {
            iterations: max_phase_iterations,
        })
    }

    /// Recovers a phase that has no positive interior equilibrium.
    ///
    /// A fixed active set may fail even though the physical equilibrium is
    /// valid on its boundary (for example, all liquid evaporates). We only
    /// retry by removing an active phase when its current continuation seed
    /// carries positive inventory;
    /// the alternative active set solves successfully; and its stability
    /// report places the removed phase above the keep threshold. This narrow
    /// contract prevents the recovery path from hiding unrelated numerical or
    /// data errors behind a convenient phase deletion.
    ///
    /// Iterates over active phases with positive current inventory, removes
    /// each one in turn, and re-solves the reduced active set. Recovery is
    /// accepted only when the reduced solve succeeds and the removed phase's
    /// accepted minimum TPD exceeds the keep threshold. Returns `None` when no
    /// phase satisfies all recovery preconditions.
    fn recover_active_boundary<F>(
        &mut self,
        iteration: usize,
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
            // The prepared problem retains the original input composition,
            // but a retargeted temperature sweep may have activated this
            // phase at an earlier accepted point. Boundary recovery must use
            // that accepted continuation state, not stale construction-time
            // inventory, otherwise a newly created phase can never be
            // considered for later disappearance.
            let continuation_phase_moles = seed
                .iter()
                .zip(species_phase.iter())
                .filter(|(_, phase)| **phase == phase_index)
                .map(|(log_moles, _)| log_moles.exp())
                .sum::<f64>();
            if !continuation_phase_moles.is_finite()
                || continuation_phase_moles <= self.phase_manager.phase_eps
            {
                continue;
            }

            let phase = PhaseIndex::new(phase_index, active.len())?;
            let mut reduced_phase_set = phase_set.clone();
            reduced_phase_set.deactivate(phase);
            reduced_phase_set.settle_transitions();
            let reduced_active = reduced_phase_set.active_mask();
            self.diagnostics.record(
                EquilibriumDiagnosticsMode::Detailed,
                EquilibriumDiagnosticEvent::RecoveryProbeStarted {
                    iteration,
                    removed_phase_index: phase_index,
                    attempted_phase_set: reduced_phase_set.clone(),
                },
            );
            let candidate = match solve_active_set(
                self,
                &reduced_active,
                seed,
                species_phase,
                full_element_totals,
            ) {
                Ok(candidate) => candidate,
                Err(error) => {
                    self.record_recovery_probe_rejected(iteration, phase_index, error.to_string());
                    continue;
                }
            };
            let stability = match compute_phase_stability_reports_with_geometry_cache(
                &candidate.log_moles,
                &candidate.stability_gibbs,
                self.prepared.problem().phases(),
                species_phase,
                self.prepared.problem().element_composition(),
                candidate.conditions.temperature(),
                candidate.conditions.pressure(),
                candidate.conditions.reference_pressure(),
                &reduced_phase_set,
                &mut self.phase_stability_geometry_cache,
            ) {
                Ok(stability) => stability,
                Err(error) => {
                    self.record_recovery_probe_rejected(iteration, phase_index, error.to_string());
                    continue;
                }
            };
            let phase_total =
                compute_phase_totals(&candidate.log_moles, species_phase)[phase_index];
            let Some(minimum_tpd) = stability[phase_index].minimum_tpd else {
                self.record_recovery_probe_rejected(
                    iteration,
                    phase_index,
                    "the removed phase has no finite TPD minimum".to_string(),
                );
                continue;
            };
            let (_, dg_keep) = self
                .phase_manager
                .thresholds_at(candidate.conditions.temperature())?;
            if phase_total < self.phase_manager.phase_eps && minimum_tpd > dg_keep {
                return Ok(Some(PreparedBoundaryRecovery {
                    phase,
                    candidate,
                    phase_set: reduced_phase_set,
                    stability,
                }));
            }
            self.record_recovery_probe_rejected(
                iteration,
                phase_index,
                format!(
                    "boundary contract failed: total={phase_total:.3e} mol, minimum_tpd={minimum_tpd:.3e} J/mol, keep_threshold={dg_keep:.3e} J/mol"
                ),
            );
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
        &mut self,
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
        self.diagnostics.record(
            EquilibriumDiagnosticsMode::Summary,
            EquilibriumDiagnosticEvent::SolveAccepted {
                final_phase_set: phase_set.clone(),
                outer_iterations: iterations,
                transition_count: phase_control_report.transitions.len(),
            },
        );
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
            diagnostics: self.diagnostics.finish(),
        })
    }

    fn record_stability_diagnostics(
        &mut self,
        iteration: usize,
        dg_create: f64,
        dg_keep: f64,
        phase_totals: &[f64],
        stability: &[PhaseStabilityReport],
    ) {
        if !self.diagnostics.mode().captures_lifecycle() {
            return;
        }
        let phases = stability
            .iter()
            .map(|report| PhaseStabilityDiagnostic {
                phase_index: report.phase.index(),
                active: report.active,
                phase_total_moles: phase_totals[report.phase.index()],
                minimum_tpd_j_per_mol: report.minimum_tpd,
                incipient_composition: report.incipient_composition.clone(),
            })
            .collect();
        self.diagnostics.record(
            EquilibriumDiagnosticsMode::PhaseLifecycle,
            EquilibriumDiagnosticEvent::StabilityEvaluated {
                iteration,
                dg_create_j_per_mol: dg_create,
                dg_keep_j_per_mol: dg_keep,
                phases,
            },
        );
    }

    fn record_transition_diagnostics(
        &mut self,
        iteration: usize,
        activated: &[PhaseIndex],
        deactivated: &[PhaseIndex],
        reason: PhaseTransitionReason,
        previous_phase_set: PhaseSet,
        new_phase_set: PhaseSet,
        incipient_composition: Option<Vec<f64>>,
    ) {
        self.diagnostics.record(
            EquilibriumDiagnosticsMode::PhaseLifecycle,
            EquilibriumDiagnosticEvent::TransitionAccepted {
                iteration,
                activated_phase_indices: activated.iter().map(|phase| phase.index()).collect(),
                deactivated_phase_indices: deactivated.iter().map(|phase| phase.index()).collect(),
                reason,
                previous_phase_set,
                new_phase_set,
                incipient_composition,
            },
        );
    }

    fn record_recovery_probe_rejected(
        &mut self,
        iteration: usize,
        phase_index: usize,
        message: String,
    ) {
        self.diagnostics.record(
            EquilibriumDiagnosticsMode::Detailed,
            EquilibriumDiagnosticEvent::RecoveryProbeRejected {
                iteration,
                removed_phase_index: phase_index,
                message,
            },
        );
    }

    fn reject_repeated_phase_set(
        &mut self,
        visited: &mut HashSet<PhaseSet>,
        phase_set: &PhaseSet,
        iteration: usize,
    ) -> Result<(), ReactionExtentError> {
        match reject_repeated_phase_set(visited, phase_set, iteration) {
            Ok(()) => Ok(()),
            Err(error) => {
                self.diagnostics.record(
                    EquilibriumDiagnosticsMode::PhaseLifecycle,
                    EquilibriumDiagnosticEvent::PhaseControlCycleDetected {
                        iteration,
                        repeated_phase_set: phase_set.clone(),
                    },
                );
                Err(error)
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_diagnostics::{
        EquilibriumDiagnosticEvent, EquilibriumDiagnosticsMode, EquilibriumDiagnosticsOptions,
    };
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_log_moles::{
        Phase, PhaseKind, Solvers,
    };
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_problem::EquilibriumConditions;
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_solver_policy::{
        SolverAttemptOutcome, SolverAttemptReport, SolverBackend, SolverPolicy,
    };
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_workflows::InitialPhaseSet;
    use nalgebra::DMatrix;
    use std::collections::HashSet;
    use std::rc::Rc;
    use std::sync::mpsc;

    fn runner_with_accepted_continuation() -> PreparedPhaseControlRunner {
        let initial_moles = vec![1.0, 1e-20];
        let problem = EquilibriumProblem::new(
            vec!["O2".to_string(), "O".to_string()],
            initial_moles.clone(),
            LogMolesInitialGuess::from_moles(&initial_moles, 1e-30).unwrap(),
            DMatrix::from_row_slice(2, 1, &[2.0, 1.0]),
            vec![Rc::new(|_| 0.0), Rc::new(|_| 0.0)],
            vec![Phase {
                kind: PhaseKind::IdealGas,
                species: vec![0, 1],
            }],
            EquilibriumConditions::new(1000.0, 101_325.0, 101_325.0).unwrap(),
        )
        .unwrap();
        let mut runner = PreparedPhaseControlRunner::new(problem, Vec::new(), false).unwrap();
        let seed = LogMolesInitialGuess::from_moles(&[0.75, 0.25], 1e-30).unwrap();
        runner
            .retarget_numeric(
                EquilibriumConditions::new(1100.0, 101_325.0, 101_325.0).unwrap(),
                seed,
                vec![Rc::new(|_| 0.0), Rc::new(|_| 0.0)],
            )
            .unwrap();
        let phase_set =
            PhaseSet::from_policy(&InitialPhaseSet::AllCandidatePhases, &[true]).unwrap();
        runner.set_continuation_phase_set(phase_set).unwrap();
        runner
    }

    #[test]
    fn failed_lifecycle_attempt_restores_accepted_continuation_state() {
        let mut runner = runner_with_accepted_continuation();
        let expected_seed = runner.continuation_seed.clone();
        let expected_phase_set = runner.continuation_phase_set.clone();

        let result = runner.solve_with_fixed_active_solver(|_, _, _, _, _| {
            Err(ReactionExtentError::InvalidProblem {
                field: "injected_phase_control_failure",
                message: "test failure after continuation was consumed".to_string(),
            })
        });

        assert!(result.is_err());
        assert_eq!(runner.continuation_seed, expected_seed);
        assert_eq!(runner.continuation_phase_set, expected_phase_set);
    }

    #[test]
    fn failed_lifecycle_attempt_streams_rejection_and_rollback_evidence() {
        let (sender, receiver) = mpsc::channel();
        let mut runner = runner_with_accepted_continuation();
        runner.set_diagnostics_options(
            EquilibriumDiagnosticsOptions::enabled(EquilibriumDiagnosticsMode::Detailed)
                .with_sink(move |event| sender.send(event).unwrap()),
        );

        let result = runner.solve_with_fixed_active_solver(|_, _, _, _, _| {
            Err(ReactionExtentError::InvalidProblem {
                field: "injected_phase_control_failure",
                message: "test lifecycle failure".to_string(),
            })
        });

        assert!(result.is_err());
        let events = receiver.try_iter().collect::<Vec<_>>();
        assert!(events.iter().any(|event| matches!(
            event,
            EquilibriumDiagnosticEvent::ActiveSetCandidateRejected { .. }
        )));
        assert!(events.iter().any(|event| matches!(
            event,
            EquilibriumDiagnosticEvent::ContinuationRestored {
                retained_seed: true,
                ..
            }
        )));
        assert!(events.iter().any(|event| matches!(
            event,
            EquilibriumDiagnosticEvent::SolveFailed {
                continuation_restored: true,
                ..
            }
        )));
    }

    #[test]
    fn repeated_phase_set_streams_typed_cycle_evidence() {
        let (sender, receiver) = mpsc::channel();
        let mut runner = runner_with_accepted_continuation();
        runner.set_diagnostics_options(
            EquilibriumDiagnosticsOptions::enabled(EquilibriumDiagnosticsMode::PhaseLifecycle)
                .with_sink(move |event| sender.send(event).unwrap()),
        );
        let phase_set = runner
            .continuation_phase_set
            .clone()
            .expect("fixture must retain one accepted phase set");
        let mut visited = HashSet::from([phase_set.clone()]);

        let error = runner
            .reject_repeated_phase_set(&mut visited, &phase_set, 3)
            .expect_err("a repeated settled phase set must stop the outer loop");

        assert!(matches!(
            error,
            ReactionExtentError::PhaseControlCycleDetected { iteration: 3, .. }
        ));
        assert!(receiver.try_iter().any(|event| matches!(
            event,
            EquilibriumDiagnosticEvent::PhaseControlCycleDetected { iteration: 3, .. }
        )));
    }

    fn accepted_validation(min_moles: f64) -> EquilibriumCandidateReport {
        EquilibriumCandidateReport {
            residual_l2_norm: 0.0,
            residual_rms: 0.0,
            max_abs_residual: 0.0,
            raw_residual_l2_norm: 0.0,
            raw_residual_rms: 0.0,
            raw_max_abs_residual: 0.0,
            max_abs_element_balance_error: 0.0,
            reaction_affinity_l2_norm: 0.0,
            max_abs_reaction_affinity: 0.0,
            min_moles,
        }
    }

    fn accepted_report() -> EquilibriumSolveReport {
        let backend = SolverBackend::Legacy(Solvers::NR);
        EquilibriumSolveReport {
            policy: SolverPolicy::Single(backend),
            attempts: vec![SolverAttemptReport {
                backend,
                outcome: SolverAttemptOutcome::Accepted,
                metrics: None,
            }],
            accepted_backend: backend,
        }
    }

    fn candidate(
        log_moles: Vec<f64>,
        solved_active_mask: Vec<bool>,
        gibbs: &[GibbsFn],
        temperature: f64,
    ) -> PreparedActiveSetCandidate {
        PreparedActiveSetCandidate {
            log_moles,
            solved_active_mask,
            validation_report: accepted_validation(1.0e-12),
            solve_report: accepted_report(),
            keq_validation_status: None,
            stability_gibbs: gibbs.to_vec(),
            conditions: EquilibriumConditions::new(temperature, 101_325.0, 101_325.0)
                .expect("synthetic P,H candidate conditions must be valid"),
            projection_build: Duration::ZERO,
            formulation_build: Duration::ZERO,
            validation_duration: Duration::ZERO,
            rst_symbolic_reused: false,
        }
    }

    /// A minimal phase topology for deterministic P,H lifecycle tests.
    ///
    /// The gas `A` fixes a one-element reference assemblage. `B` and `C`
    /// form an initially absent ideal solution with the same elemental
    /// direction. The caller controls only their standard Gibbs values, which
    /// gives a wide separation between stable and unstable TPD fixtures.
    fn runner_with_inactive_ideal_solution(gibbs: Vec<GibbsFn>) -> PreparedPhaseControlRunner {
        let physical_initial_moles = vec![1.0, 0.0, 0.0];
        let numerical_seed = LogMolesInitialGuess::from_moles(&[1.0, 1.0e-30, 1.0e-30], 1.0e-30)
            .expect("trace numerical coordinates must be valid");
        let problem = EquilibriumProblem::new(
            vec!["A".to_string(), "B".to_string(), "C".to_string()],
            physical_initial_moles,
            numerical_seed,
            DMatrix::from_row_slice(3, 1, &[1.0, 1.0, 1.0]),
            gibbs,
            vec![
                Phase {
                    kind: PhaseKind::IdealGas,
                    species: vec![0],
                },
                Phase {
                    kind: PhaseKind::IdealSolution,
                    species: vec![1, 2],
                },
            ],
            EquilibriumConditions::new(300.0, 101_325.0, 101_325.0)
                .expect("synthetic P,H conditions must be valid"),
        )
        .expect("synthetic P,H lifecycle topology must be valid");
        PreparedPhaseControlRunner::new(problem, Vec::new(), false)
            .expect("synthetic P,H lifecycle runner must prepare")
    }

    fn stable_solution_gibbs() -> Vec<GibbsFn> {
        vec![
            Rc::new(|_| 0.0),
            Rc::new(|_| 10_000.0),
            Rc::new(|_| 10_000.0),
        ]
    }

    fn unstable_solution_gibbs() -> Vec<GibbsFn> {
        // At 450 K this produces an interior minimizer close to [0.79, 0.21],
        // far enough from a neutral [0.5, 0.5] probe seed to catch leakage.
        vec![Rc::new(|_| 0.0), Rc::new(|_| -5_000.0), Rc::new(|_| 0.0)]
    }

    #[test]
    fn reduced_lifecycle_candidate_is_accepted_without_probe_or_transition() {
        let gibbs = stable_solution_gibbs();
        let mut runner = runner_with_inactive_ideal_solution(gibbs.clone());
        let mut calls = Vec::new();

        let outcome = runner
            .solve_with_fixed_active_solver(|_, active, seed, _, _| {
                calls.push((active.to_vec(), seed.to_vec()));
                Ok(candidate(
                    vec![0.0, 1.0e-40_f64.ln(), 1.0e-40_f64.ln()],
                    active.to_vec(),
                    &gibbs,
                    450.0,
                ))
            })
            .expect("a valid reduced candidate must not require a recovery probe");

        assert_eq!(calls.len(), 1);
        assert_eq!(calls[0].0, vec![true, false]);
        assert!(outcome.phase_control_report.transitions.is_empty());
        assert_eq!(
            outcome.phase_control_report.final_phase_set.active_mask(),
            vec![true, false]
        );
        assert_eq!(outcome.solution.conditions().temperature(), 450.0);
        assert_eq!(
            outcome.solution.validation().max_abs_element_balance_error,
            0.0
        );
    }

    #[test]
    fn wider_probe_with_stable_tpd_is_rejected_without_leaking_a_phase_transition() {
        let gibbs = stable_solution_gibbs();
        let mut runner = runner_with_inactive_ideal_solution(gibbs.clone());
        let mut calls = Vec::new();

        let error = runner
            .solve_with_fixed_active_solver(|_, active, seed, _, _| {
                calls.push((active.to_vec(), seed.to_vec()));
                Ok(candidate(
                    vec![0.0, 0.5_f64.ln(), 0.5_f64.ln()],
                    vec![true, true],
                    &gibbs,
                    450.0,
                ))
            })
            .expect_err("a wider numerical probe without negative TPD must not publish");

        assert_eq!(
            calls.len(),
            1,
            "the rejected probe must not restart a solve"
        );
        assert_eq!(calls[0].0, vec![true, false]);
        assert!(
            error.to_string().contains("TPD creation threshold"),
            "the physical rejection must name the TPD boundary: {error}"
        );
        assert_eq!(
            runner.continuation_phase_set, None,
            "a rejected wider probe must not publish its phase set as continuation"
        );
    }

    #[test]
    fn wider_probe_activation_uses_tpd_minimizer_not_neutral_probe_seed() {
        let gibbs = unstable_solution_gibbs();
        let mut runner = runner_with_inactive_ideal_solution(gibbs.clone());
        let mut calls = Vec::new();

        let outcome = runner
            .solve_with_fixed_active_solver(|_, active, seed, _, _| {
                calls.push((active.to_vec(), seed.to_vec()));
                match calls.len() {
                    1 => Ok(candidate(
                        // The wider numerical branch deliberately uses the
                        // neutral probe composition [0.5, 0.5].
                        vec![0.0, 0.5_f64.ln(), 0.5_f64.ln()],
                        vec![true, true],
                        &gibbs,
                        450.0,
                    )),
                    2 => Ok(candidate(
                        vec![0.0, 0.79_f64.ln(), 0.21_f64.ln()],
                        active.to_vec(),
                        &gibbs,
                        450.0,
                    )),
                    _ => panic!("one TPD activation must require exactly one fixed-set restart"),
                }
            })
            .expect("negative TPD must activate the ideal solution through a restart");

        assert_eq!(calls.len(), 2);
        assert_eq!(calls[0].0, vec![true, false]);
        assert_eq!(calls[1].0, vec![true, true]);
        let restart_solution_total = calls[1].1[1].exp() + calls[1].1[2].exp();
        let restart_fraction_b = calls[1].1[1].exp() / restart_solution_total;
        assert!(
            (restart_fraction_b - 0.5).abs() > 0.1,
            "the physical restart must not retain the neutral probe composition"
        );

        let transition = outcome
            .phase_control_report
            .transitions
            .first()
            .expect("TPD activation must retain transition evidence");
        let incipient = transition
            .incipient_composition
            .as_ref()
            .expect("TPD activation must retain its minimizer composition");
        assert_eq!(transition.activated.len(), 1);
        assert!(matches!(
            transition.reason,
            PhaseTransitionReason::UnstableInactivePhase { minimum_tpd } if minimum_tpd < 0.0
        ));
        assert!((incipient.iter().sum::<f64>() - 1.0).abs() < 1.0e-12);
        assert!(incipient[0] > 0.7 && incipient[1] < 0.3);
        assert!(
            (restart_fraction_b - incipient[0]).abs() < 1.0e-10,
            "restart seed ratios must be the TPD minimizer, not the probe composition"
        );
        assert_eq!(
            outcome.phase_control_report.final_phase_set.active_mask(),
            vec![true, true]
        );
    }

    fn runner_with_accepted_unstable_solution_continuation(
        gibbs: Vec<GibbsFn>,
    ) -> PreparedPhaseControlRunner {
        let mut runner = runner_with_inactive_ideal_solution(gibbs.clone());
        let accepted_seed = LogMolesInitialGuess::from_moles(&[0.7, 1.0e-30, 1.0e-30], 1.0e-30)
            .expect("accepted continuation seed must be valid");
        runner
            .retarget_numeric(
                EquilibriumConditions::new(425.0, 101_325.0, 101_325.0)
                    .expect("accepted continuation conditions must be valid"),
                accepted_seed,
                gibbs,
            )
            .expect("accepted continuation must retarget");
        runner
            .set_continuation_phase_set(
                PhaseSet::from_policy(
                    &InitialPhaseSet::Explicit {
                        active: vec![PhaseIndex::new(0, 2).unwrap()],
                        excluded: Vec::new(),
                    },
                    &[true, false],
                )
                .expect("accepted continuation phase set must be valid"),
            )
            .expect("accepted continuation phase set must install");
        runner
    }

    #[test]
    fn transition_failpoint_rolls_back_after_transition_diagnostics() {
        let gibbs = unstable_solution_gibbs();
        let mut runner = runner_with_accepted_unstable_solution_continuation(gibbs.clone());
        let expected_seed = runner.continuation_seed.clone();
        let expected_phase_set = runner.continuation_phase_set.clone();
        let (sender, receiver) = mpsc::channel();
        runner.set_diagnostics_options(
            EquilibriumDiagnosticsOptions::enabled(EquilibriumDiagnosticsMode::PhaseLifecycle)
                .with_sink(move |event| sender.send(event).unwrap()),
        );
        runner.arm_test_failpoint_once(TestFailPoint::AfterPhaseTransitionBeforeCommit);

        let error = runner
            .solve_with_fixed_active_solver(|_, _active, _, _, _| {
                Ok(candidate(
                    vec![0.0, 0.5_f64.ln(), 0.5_f64.ln()],
                    vec![true, true],
                    &gibbs,
                    450.0,
                ))
            })
            .expect_err("the transition boundary failpoint must abort the attempt");

        assert!(
            error
                .to_string()
                .contains("after phase transition before commit")
        );
        assert_eq!(runner.continuation_seed, expected_seed);
        assert_eq!(runner.continuation_phase_set, expected_phase_set);

        let events = receiver.try_iter().collect::<Vec<_>>();
        assert!(
            events.iter().any(|event| matches!(
                event,
                EquilibriumDiagnosticEvent::TransitionAccepted { .. }
            ))
        );
        assert!(events.iter().any(|event| matches!(
            event,
            EquilibriumDiagnosticEvent::ContinuationRestored {
                retained_seed: true,
                retained_phase_set: Some(_),
                ..
            }
        )));
        assert!(events.iter().any(|event| matches!(
            event,
            EquilibriumDiagnosticEvent::SolveFailed {
                continuation_restored: true,
                ..
            }
        )));
    }

    #[test]
    fn failed_restart_after_tpd_probe_restores_the_previous_continuation() {
        let gibbs = unstable_solution_gibbs();
        let mut runner = runner_with_accepted_unstable_solution_continuation(gibbs.clone());
        let expected_seed = runner.continuation_seed.clone();
        let expected_phase_set = runner.continuation_phase_set.clone();
        let mut calls = 0usize;

        let error = runner
            .solve_with_fixed_active_solver(|_, _active, _, _, _| {
                calls += 1;
                if calls == 1 {
                    Ok(candidate(
                        vec![0.0, 0.5_f64.ln(), 0.5_f64.ln()],
                        vec![true, true],
                        &gibbs,
                        450.0,
                    ))
                } else {
                    Err(ReactionExtentError::InvalidProblem {
                        field: "injected_restart_failure",
                        message: "the TPD-selected fixed set failed before acceptance".to_string(),
                    })
                }
            })
            .expect_err("a failed restart after a probe must not publish a phase transition");

        assert!(
            calls >= 2,
            "the TPD-selected fixed set must be attempted before bounded boundary recovery"
        );
        assert!(error.to_string().contains("injected_restart_failure"));
        assert_eq!(runner.continuation_seed, expected_seed);
        assert_eq!(runner.continuation_phase_set, expected_phase_set);
    }
}
