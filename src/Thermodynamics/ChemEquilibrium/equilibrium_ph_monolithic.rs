//! Fixed-active-set monolithic `P,H` nonlinear runner.
//!
//! The runner solves `[ln(n_1), ..., ln(n_m), theta_T]` in one nonlinear
//! system. It deliberately owns no phase lifecycle: activation, destruction,
//! hysteresis, and transactional publication stay at the phase-control
//! boundary. This separation lets the numerical formulation be tested without
//! mutating a broader multiphase workflow.
//!
//! RustedSciThe backends receive a dedicated symbolic `(ln(n), theta_T)`
//! payload. The retained legacy solvers consume the independently implemented
//! analytical Jacobian from `equilibrium_ph_formulation`, which remains a
//! valuable validation and fallback route.

use std::time::{Duration, Instant};

use nalgebra::DMatrix;

use crate::Thermodynamics::ChemEquilibrium::equilibrium_active_set::ActiveSetProjection;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_backend_adapter::{
    EquilibriumNonlinearBackend, solve_backend_cascade_with_control,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_constraints::{
    EnthalpyScale, TemperatureBounds,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_ids::PhaseIndex;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_log_moles::EquilibriumSolverSettings;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_nonlinear::ReactionExtentError;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_ph_formulation::{
    PhCandidateSnapshot, PreparedPhFormulation,
};
pub(crate) use crate::Thermodynamics::ChemEquilibrium::equilibrium_ph_options::PhAcceptanceOptions as MonolithicPhAcceptanceOptions;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_ph_thermochemistry::ResolvedThermochemistry;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_problem::{
    EquilibriumConditions, EquilibriumProblem, LogMolesInitialGuess, PreparedEquilibriumProblem,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_rst_backend::prepare_rst_symbolic_ph_problem;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_solver_policy::{
    EquilibriumSolveReport, SolverBackend, SolverCascadeBudget, SolverPolicy,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_validation::{
    EquilibriumAcceptanceCriteria, EquilibriumCandidateReport, EquilibriumCandidateResiduals,
    validate_equilibrium_candidate,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_workflows::{
    PhaseTotalSeedPolicy, seed_activated_phase_with_composition,
};
use crate::Thermodynamics::ChemEquilibrium::prepared_phase_control_runner::{
    PreparedActiveSetCandidate, PreparedPhaseControlRunner,
};

/// Accepted fixed-active-set monolithic P,H candidate.
// This is intentionally wider than the first runner consumer. The facade
// will need the original unknowns and P,T validation evidence when it starts
// publishing a unified nested/monolithic diagnostic report.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub(crate) struct MonolithicPhSolveOutcome {
    /// Full nonlinear coordinate `[ln(n), theta_T]` accepted by the backend cascade.
    pub(crate) unknowns: Vec<f64>,
    /// Physical values and residuals evaluated at the accepted coordinate.
    pub(crate) snapshot: PhCandidateSnapshot,
    /// Common fixed-P,T conservation and affinity validation report.
    pub(crate) pt_validation: EquilibriumCandidateReport,
    /// Deterministic record of attempted numerical backends.
    pub(crate) solve_report: EquilibriumSolveReport,
}

/// Immutable numerical runner for one already fixed active phase set.
#[derive(Clone)]
pub(crate) struct PreparedMonolithicPhRunner {
    formulation: PreparedPhFormulation,
    settings: EquilibriumSolverSettings,
    acceptance_options: MonolithicPhAcceptanceOptions,
}

impl PreparedMonolithicPhRunner {
    /// Creates a runner after validating common solver and enthalpy controls.
    pub(crate) fn new(
        formulation: PreparedPhFormulation,
        settings: EquilibriumSolverSettings,
        acceptance_options: MonolithicPhAcceptanceOptions,
    ) -> Result<Self, ReactionExtentError> {
        settings.validate()?;
        Ok(Self {
            formulation,
            settings,
            acceptance_options,
        })
    }

    /// Solves the coupled composition-temperature equations from one physical
    /// temperature seed. No state is published on failure.
    pub(crate) fn solve_from_temperature_seed(
        &self,
        temperature_seed: f64,
    ) -> Result<MonolithicPhSolveOutcome, ReactionExtentError> {
        let initial_unknowns = self.formulation.initial_unknowns(temperature_seed)?;
        self.solve_from_initial_unknowns(initial_unknowns, None)
    }

    /// Solves from an accepted `[log-moles, temperature]` continuation state.
    /// An optional prepared RST problem is reused when supplied; legacy
    /// backends continue to use the same analytical formulation.
    pub(crate) fn solve_from_log_moles_and_temperature_seed(
        &self,
        log_moles: &crate::Thermodynamics::ChemEquilibrium::equilibrium_problem::
            LogMolesInitialGuess,
        temperature_seed: f64,
        rst_problem: Option<
            &crate::Thermodynamics::ChemEquilibrium::equilibrium_rst_backend::RstPreparedProblem,
        >,
    ) -> Result<MonolithicPhSolveOutcome, ReactionExtentError> {
        let initial_unknowns = self
            .formulation
            .initial_unknowns_from_log_moles(log_moles, temperature_seed)?;
        self.solve_from_initial_unknowns(initial_unknowns, rst_problem)
    }

    /// Core dispatch for the monolithic `[log-moles, theta_T]` unknown vector.
    ///
    /// Resolves the solver policy and cascade budget, then delegates to the
    /// common backend adapter. The optional RST symbolic problem is forwarded
    /// when the policy selects a RustedSciThe backend; legacy backends receive
    /// the analytical formulation from [`PhFormulation`].
    fn solve_from_initial_unknowns(
        &self,
        initial_unknowns: Vec<f64>,
        prepared_rst_problem: Option<
            &crate::Thermodynamics::ChemEquilibrium::equilibrium_rst_backend::RstPreparedProblem,
        >,
    ) -> Result<MonolithicPhSolveOutcome, ReactionExtentError> {
        let policy = self.monolithic_policy()?;
        let backends = policy.ordered_backends();
        let budget = self.settings.solver_budget.unwrap_or_else(|| {
            SolverCascadeBudget::new(
                backends.len(),
                self.settings.solver_params.max_iter,
                self.settings
                    .solver_params
                    .max_iter
                    .saturating_mul(backends.len()),
            )
        });
        let backend_refs: Vec<&dyn EquilibriumNonlinearBackend> = backends
            .iter()
            .map(|backend| backend as &dyn EquilibriumNonlinearBackend)
            .collect();
        let owned_rst_problem = if prepared_rst_problem.is_none()
            && backends
                .iter()
                .any(|backend| matches!(backend, SolverBackend::RustedSciThe(_)))
        {
            Some(prepare_rst_symbolic_ph_problem(&self.formulation)?)
        } else {
            None
        };
        let rst_problem = prepared_rst_problem.or(owned_rst_problem.as_ref());
        let residual = |unknowns: &[f64]| self.formulation.residual(unknowns);
        let jacobian = |unknowns: &[f64]| self.formulation.jacobian(unknowns);
        let feasible = |unknowns: &[f64]| self.formulation.is_feasible(unknowns);
        let solver_tolerance = self.settings.solver_params.tol;
        let element_tolerance = 10.0 * solver_tolerance.max(1.0e-8);
        let validate_candidate = |unknowns: &[f64]| {
            let snapshot = self.formulation.candidate_snapshot(unknowns)?;
            let validation = validate_equilibrium_candidate(
                EquilibriumCandidateResiduals {
                    log_moles: &snapshot.log_moles,
                    raw_residual: &snapshot.raw_pt_residual,
                    acceptance_residual: &snapshot.scaled_pt_residual,
                },
                EquilibriumAcceptanceCriteria::new(
                    solver_tolerance,
                    element_tolerance,
                    solver_tolerance,
                )?
                .with_element_balance_relative_tolerance(element_tolerance)?,
                self.formulation
                    .prepared_pt()
                    .problem()
                    .element_composition(),
                self.formulation.prepared_pt().element_totals(),
            )?;
            if !self
                .acceptance_options
                .accepts_enthalpy_error(snapshot.enthalpy_error, self.formulation.enthalpy_scale())
            {
                return Err(ReactionExtentError::InvalidCandidate {
                    field: "candidate_enthalpy",
                    message: format!(
                        "enthalpy error {} J exceeds the P,H acceptance limit",
                        snapshot.enthalpy_error
                    ),
                });
            }
            Ok(validation)
        };
        let (unknowns, pt_validation, solve_report) = solve_backend_cascade_with_control(
            &backend_refs,
            initial_unknowns,
            &residual,
            Some(&jacobian as &dyn Fn(&[f64]) -> Result<DMatrix<f64>, ReactionExtentError>),
            &feasible,
            &validate_candidate,
            policy,
            budget,
            &self.settings.solver_params,
            rst_problem,
            None,
        )?;
        let snapshot = self.formulation.candidate_snapshot(&unknowns)?;
        Ok(MonolithicPhSolveOutcome {
            unknowns,
            snapshot,
            pt_validation,
            solve_report,
        })
    }

    /// Resolves the effective solver policy for the monolithic `P,H` runner.
    ///
    /// Returns the explicit policy when one has been installed, or falls back
    /// to [`SolverPolicy::legacy_default`] using the preferred legacy backend.
    /// This keeps the monolithic path consistent with the fixed-`P,T` cascade
    /// contract without requiring every caller to supply a policy.
    fn monolithic_policy(&self) -> Result<SolverPolicy, ReactionExtentError> {
        let policy = self
            .settings
            .solver_policy
            .clone()
            .unwrap_or_else(|| SolverPolicy::legacy_default(self.settings.solver));
        Ok(policy)
    }
}

/// Solves one active mask for the shared phase-control lifecycle.
///
/// The callback boundary deliberately carries only a fixed active mask and a
/// log-mole restart seed. It projects the numerical problem and the resolved
/// thermochemistry together, so component order cannot drift between the
/// conservation matrix, Gibbs functions, enthalpy, or Cp. The full candidate
/// returned to the lifecycle contains Gibbs capabilities at the solved
/// temperature for phase-stability decisions.
pub(crate) fn solve_monolithic_active_set_candidate(
    runner: &mut PreparedPhaseControlRunner,
    active: &[bool],
    full_seed: &[f64],
    species_phase: &[usize],
    full_element_totals: &[f64],
    thermochemistry: &ResolvedThermochemistry,
    temperature_bounds: TemperatureBounds,
    target_enthalpy: f64,
    enthalpy_scale: EnthalpyScale,
    temperature_seed: &mut f64,
    acceptance_options: &MonolithicPhAcceptanceOptions,
) -> Result<PreparedActiveSetCandidate, ReactionExtentError> {
    // A P,H root may lie on a branch where a phase absent from the initial
    // inventory becomes thermodynamically required. The ordinary active-set
    // lifecycle starts such a phase inactive, but a fixed monolithic system
    // cannot discover that branch if the current reduced problem has no
    // enthalpy root at all. Probe the full trace-seeded layout as a recovery
    // candidate; the phase manager still owns the subsequent activation or
    // rejection decision and the final physical publication remains
    // transactional.
    if active.iter().any(|is_active| !is_active) {
        let all_active = vec![true; active.len()];
        match solve_monolithic_active_set_candidate_inner(
            runner,
            active,
            full_seed,
            species_phase,
            full_element_totals,
            thermochemistry,
            temperature_bounds,
            target_enthalpy,
            enthalpy_scale,
            temperature_seed,
            acceptance_options,
        ) {
            Ok(candidate) => return Ok(candidate),
            Err(primary_error) => {
                // A single trace floor can be too far from a condensed-phase
                // branch for a legacy Newton-like backend. Keep this bounded
                // and deterministic: a few physically positive continuation
                // seeds are cheaper and more honest than pretending that one
                // failed nonlinear start proves the phase set impossible.
                for fraction in [1.0e-8, 1.0e-4, 1.0e-2, 1.0e-1] {
                    let mut probe_seed = full_seed.to_vec();
                    for phase_index in 0..all_active.len() {
                        if !active[phase_index] {
                            let component_count = species_phase
                                .iter()
                                .filter(|&&phase| phase == phase_index)
                                .count();
                            if component_count == 0 {
                                return Err(ReactionExtentError::InvalidProblem {
                                    field: "monolithic_ph_probe_seed",
                                    message: format!(
                                        "inactive phase {phase_index} has no declared components"
                                    ),
                                });
                            }
                            // This is a numerical all-active recovery probe,
                            // not a physical phase-appearance decision. The
                            // later TPD lifecycle replaces this neutral seed
                            // with its accepted minimizer composition.
                            let neutral_composition =
                                vec![1.0 / component_count as f64; component_count];
                            seed_activated_phase_with_composition(
                                &mut probe_seed,
                                PhaseIndex::new(phase_index, all_active.len())?,
                                species_phase,
                                &neutral_composition,
                                PhaseTotalSeedPolicy::RelativeToSystemTotal {
                                    fraction,
                                    minimum: 1.0e-12,
                                },
                            )?;
                        }
                    }
                    if let Ok(mut candidate) = solve_monolithic_active_set_candidate_inner(
                        runner,
                        &all_active,
                        &probe_seed,
                        species_phase,
                        full_element_totals,
                        thermochemistry,
                        temperature_bounds,
                        target_enthalpy,
                        enthalpy_scale,
                        temperature_seed,
                        acceptance_options,
                    ) {
                        candidate.solved_active_mask = all_active.clone();
                        return Ok(candidate);
                    }
                }
                return Err(primary_error);
            }
        }
    }

    solve_monolithic_active_set_candidate_inner(
        runner,
        active,
        full_seed,
        species_phase,
        full_element_totals,
        thermochemistry,
        temperature_bounds,
        target_enthalpy,
        enthalpy_scale,
        temperature_seed,
        acceptance_options,
    )
}

fn solve_monolithic_active_set_candidate_inner(
    runner: &mut PreparedPhaseControlRunner,
    active: &[bool],
    full_seed: &[f64],
    species_phase: &[usize],
    full_element_totals: &[f64],
    thermochemistry: &ResolvedThermochemistry,
    temperature_bounds: TemperatureBounds,
    target_enthalpy: f64,
    enthalpy_scale: EnthalpyScale,
    temperature_seed: &mut f64,
    acceptance_options: &MonolithicPhAcceptanceOptions,
) -> Result<PreparedActiveSetCandidate, ReactionExtentError> {
    let started = Instant::now();
    let problem = runner.prepared_problem();
    let settings = runner.solver_settings();
    let projection = ActiveSetProjection::build(
        problem.problem().phases(),
        species_phase,
        problem.problem().element_composition(),
        active,
        settings.solver_params.tol,
    )?;
    projection
        .validate_element_totals_representable(full_element_totals, settings.solver_params.tol)?;

    let reduced_log_seed = projection.project_log_moles(full_seed)?;
    let indices = projection
        .active_species
        .iter()
        .map(|species| species.index())
        .collect::<Vec<_>>();
    let reduced_components = indices
        .iter()
        .map(|&index| problem.problem().components()[index].clone())
        .collect::<Vec<_>>();
    let reduced_moles = reduced_log_seed
        .iter()
        .map(|value| value.exp())
        .collect::<Vec<_>>();
    let reduced_thermochemistry = thermochemistry.subset(&indices)?;
    // `PreparedEquilibriumProblem` still carries the historical infallible
    // Gibbs closure even though the monolithic formulation evaluates the
    // typed thermochemistry bundle directly.  Cross that boundary only with
    // a fully validated finite snapshot; never encode a source error as NaN.
    let reduced_gibbs =
        reduced_thermochemistry.gibbs_snapshot_for_legacy_boundary(*temperature_seed)?;
    let conditions = EquilibriumConditions::new(
        *temperature_seed,
        problem.problem().conditions().pressure(),
        problem.problem().conditions().reference_pressure(),
    )?;
    let reduced_totals = projection.reduced_element_totals(full_element_totals)?;
    let reduced_problem = EquilibriumProblem::new(
        reduced_components,
        reduced_moles,
        LogMolesInitialGuess::new(reduced_log_seed.clone())?,
        projection.element_composition.clone(),
        reduced_gibbs,
        projection.phases.clone(),
        conditions,
    )?;
    let prepared =
        PreparedEquilibriumProblem::new_with_element_totals(reduced_problem, Some(reduced_totals))?;
    let formulation = PreparedPhFormulation::new(
        prepared,
        reduced_thermochemistry,
        temperature_bounds,
        target_enthalpy,
        enthalpy_scale,
    )?;
    let monolithic = PreparedMonolithicPhRunner::new(formulation, settings, *acceptance_options)?;
    let outcome = monolithic.solve_from_temperature_seed(*temperature_seed)?;
    *temperature_seed = outcome.snapshot.temperature;
    let solved_conditions = EquilibriumConditions::new(
        outcome.snapshot.temperature,
        problem.problem().conditions().pressure(),
        problem.problem().conditions().reference_pressure(),
    )?;
    let log_moles = projection.scatter_log_moles(
        &outcome.snapshot.log_moles,
        crate::Thermodynamics::ChemEquilibrium::equilibrium_workflows::
            PHASE_CONTROL_TRACE_MOLE_FLOOR
            .ln(),
    )?;
    // Stability is evaluated at `outcome.snapshot.temperature` by the phase
    // lifecycle.  The snapshot adapter validates that exact point before the
    // legacy infallible callback is published to phase-control.
    let stability_gibbs =
        thermochemistry.gibbs_snapshot_for_legacy_boundary(outcome.snapshot.temperature)?;
    Ok(PreparedActiveSetCandidate {
        log_moles,
        solved_active_mask: active.to_vec(),
        validation_report: outcome.pt_validation,
        solve_report: outcome.solve_report,
        keq_validation_status: None,
        stability_gibbs,
        conditions: solved_conditions,
        projection_build: started.elapsed(),
        formulation_build: Duration::ZERO,
        validation_duration: Duration::ZERO,
        rst_symbolic_reused: false,
    })
}

#[cfg(test)]
mod tests {
    use std::rc::Rc;
    use std::sync::Arc;

    use RustedSciThe::symbolic::symbolic_engine::Expr;
    use nalgebra::DMatrix;

    use super::*;
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_constraints::{
        EnthalpyScale, TemperatureBounds,
    };
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_log_moles::{
        GibbsFn, Phase, PhaseKind, Solvers,
    };
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_ph_workflow::{
        MolarThermoFunction, ResolvedThermochemistry, ThermochemistryProvenance,
    };
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_problem::{
        EquilibriumConditions, EquilibriumProblem, LogMolesInitialGuess, PreparedEquilibriumProblem,
    };
    use crate::Thermodynamics::phase_layout::{PhaseComponentId, PhaseId};

    fn one_species_formulation(with_symbolic_capability: bool) -> PreparedPhFormulation {
        let initial_moles = vec![1.0];
        let problem = EquilibriumProblem::new(
            vec!["A".to_string()],
            initial_moles.clone(),
            LogMolesInitialGuess::from_initial_moles(&initial_moles).unwrap(),
            DMatrix::from_row_slice(1, 1, &[1.0]),
            vec![Rc::new(|_| 0.0) as GibbsFn],
            vec![Phase {
                kind: PhaseKind::IdealGas,
                species: vec![0],
            }],
            EquilibriumConditions::new(900.0, 101_325.0, 101_325.0).unwrap(),
        )
        .unwrap();
        let phase = PhaseId::new(None);
        let thermochemistry = ResolvedThermochemistry::from_functions(
            vec![ThermochemistryProvenance::new(
                PhaseComponentId::new(phase, "A"),
                "synthetic",
                "A",
                "gas",
            )],
            TemperatureBounds::new(300.0, 2_000.0).unwrap(),
            vec![Arc::new(|_| Ok(0.0)) as MolarThermoFunction],
            vec![Arc::new(|temperature| Ok(10.0 * temperature)) as MolarThermoFunction],
            vec![Some(Arc::new(|_| Ok(10.0)) as MolarThermoFunction)],
        )
        .unwrap();
        let thermochemistry = if with_symbolic_capability {
            thermochemistry
                .with_symbolic_expressions(
                    vec![Expr::Const(0.0)],
                    vec![Expr::Const(10.0) * Expr::Var("T".to_string())],
                )
                .unwrap()
        } else {
            thermochemistry
        };
        PreparedPhFormulation::new(
            PreparedEquilibriumProblem::new(problem).unwrap(),
            thermochemistry,
            TemperatureBounds::new(500.0, 1_500.0).unwrap(),
            10_000.0,
            EnthalpyScale::new(10_000.0).unwrap(),
        )
        .unwrap()
    }

    #[test]
    fn monolithic_runner_solves_composition_and_temperature_together() {
        let mut settings = EquilibriumSolverSettings::default();
        settings.solver = Solvers::NR;
        let runner = PreparedMonolithicPhRunner::new(
            one_species_formulation(false),
            settings,
            MonolithicPhAcceptanceOptions::new(1.0e-8, 1.0e-6).unwrap(),
        )
        .unwrap();

        let outcome = runner.solve_from_temperature_seed(900.0).unwrap();
        assert!((outcome.snapshot.temperature - 1_000.0).abs() < 1.0e-5);
        assert!((outcome.snapshot.moles[0] - 1.0).abs() < 1.0e-10);
        assert!(outcome.snapshot.enthalpy_error.abs() < 1.0e-5);
        assert!(matches!(
            outcome.solve_report.accepted_backend,
            SolverBackend::Legacy(_)
        ));
    }

    #[test]
    fn monolithic_runner_rejects_rst_without_symbolic_ph_capabilities() {
        let mut settings = EquilibriumSolverSettings::default();
        settings.solver_policy = Some(SolverPolicy::rusted_scithe_default());
        let runner = PreparedMonolithicPhRunner::new(
            one_species_formulation(false),
            settings,
            MonolithicPhAcceptanceOptions::new(1.0e-8, 1.0e-6).unwrap(),
        )
        .unwrap();

        assert!(matches!(
            runner.solve_from_temperature_seed(900.0),
            Err(ReactionExtentError::UnsupportedBackendCapability {
                backend,
                capability: "monolithic_p_h_symbolic_thermochemistry",
                ..
            }) if backend == "RustedSciThe"
        ));
    }

    #[test]
    fn monolithic_runner_lets_rst_differentiate_the_coupled_ph_system() {
        let mut settings = EquilibriumSolverSettings::default();
        settings.solver_policy = Some(SolverPolicy::rusted_scithe_default());
        let runner = PreparedMonolithicPhRunner::new(
            one_species_formulation(true),
            settings,
            // RST's own default termination is looser than the legacy test
            // fixture's 1e-8 scaled energy gate. This test exercises the
            // symbolic route under an explicitly compatible physical gate.
            MonolithicPhAcceptanceOptions::new(1.0e-7, 1.0e-6).unwrap(),
        )
        .unwrap();

        let outcome = runner.solve_from_temperature_seed(900.0).unwrap();
        assert!((outcome.snapshot.temperature - 1_000.0).abs() < 1.0e-3);
        assert!(outcome.snapshot.enthalpy_error.abs() < 1.0e-3);
        assert!(matches!(
            outcome.solve_report.accepted_backend,
            SolverBackend::RustedSciThe(_)
        ));
    }
}
