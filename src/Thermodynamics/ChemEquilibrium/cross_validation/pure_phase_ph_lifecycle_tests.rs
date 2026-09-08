//! I3 synthetic lifecycle tests for the independent pure-phase `P,H` work.
//!
//! The thermochemistry remains analytic and local, but the canonical side is
//! intentionally the real monolithic P,H active-set adapter plus the shared
//! phase-control runner. This proves TPD-driven publication separately from
//! the I1-I2 scalar reference route.

use std::rc::Rc;
use std::sync::Arc;

use nalgebra::DMatrix;

use crate::Thermodynamics::ChemEquilibrium::equilibrium_activity::PhaseActivityModel;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_constraints::{
    EnthalpyScale, TemperatureBounds,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_ids::PhaseIndex;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_log_moles::{
    EquilibriumSolverSettings, GibbsFn, Phase, Solvers,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_ph_monolithic::solve_monolithic_active_set_candidate;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_ph_options::PhMonolithicOptions;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_ph_thermochemistry::{
    MolarThermoFunction, ResolvedThermochemistry, ThermochemistryProvenance,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_problem::{
    EquilibriumConditions, EquilibriumProblem, LogMolesInitialGuess,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_solver_policy::SolverPolicy;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_workflows::{
    PhaseStatus, PhaseTransitionReason,
};
use crate::Thermodynamics::ChemEquilibrium::phase_boundary_validation::MOLAR_GAS_CONSTANT;
use crate::Thermodynamics::ChemEquilibrium::prepared_phase_control_runner::PreparedPhaseControlRunner;
use crate::Thermodynamics::ChemEquilibrium::pure_phase_ph_validation::{
    PurePhasePhConditions, PurePhasePhProblem, PurePhasePhValidator,
};
use crate::Thermodynamics::phase_layout::{PhaseComponentId, PhaseId};

const T_STAR: f64 = 800.0;
const REACTION_ENTHALPY: f64 = 2_000.0;
const GAS_HEAT_CAPACITY: f64 = 10.0;
const FAVORABLE_CANDIDATE_MOLES: f64 = 1.0e-2;
// The fixed-topology I2 reference is only meaningful where the favorable
// branch has a positive chemical extent.  Below roughly 787 K this analytic
// fixture's equilibrium would lie at the zero-candidate boundary instead.
const TEMPERATURE_LOWER: f64 = 795.0;
const TEMPERATURE_UPPER: f64 = 820.0;

fn thermo_function<F>(callback: F) -> MolarThermoFunction
where
    F: Fn(
            f64,
        ) -> Result<
            f64,
            crate::Thermodynamics::ChemEquilibrium::equilibrium_nonlinear::ReactionExtentError,
        > + Send
        + Sync
        + 'static,
{
    Arc::new(callback)
}

fn gas_standard_gibbs(temperature: f64) -> f64 {
    GAS_HEAT_CAPACITY * temperature * (1.0 - (temperature / T_STAR).ln())
}

fn ln_q(extent: f64) -> f64 {
    let a = 2.0 - 2.0 * extent;
    let b = 2.0 - extent;
    let total = a + b;
    -2.0 * (a / total).ln() - (b / total).ln()
}

fn favorable_ln_k(temperature: f64) -> f64 {
    ln_q(FAVORABLE_CANDIDATE_MOLES)
        + REACTION_ENTHALPY / MOLAR_GAS_CONSTANT * (1.0 / T_STAR - 1.0 / temperature)
}

fn favorable_candidate_gibbs(temperature: f64) -> f64 {
    3.0 * gas_standard_gibbs(temperature)
        - MOLAR_GAS_CONSTANT * temperature * favorable_ln_k(temperature)
}

fn favorable_candidate_enthalpy(temperature: f64) -> f64 {
    3.0 * GAS_HEAT_CAPACITY * temperature + REACTION_ENTHALPY
}

/// Returns the positive reaction extent satisfying the analytic chemical
/// equation at `T_STAR` after a constant pure-candidate Gibbs offset.
fn candidate_equilibrium_moles_at_t_star(candidate_offset: f64) -> f64 {
    let target_ln_k = favorable_ln_k(T_STAR) - candidate_offset / (MOLAR_GAS_CONSTANT * T_STAR);
    let mut lower = 0.0;
    let mut upper = 1.0 - 1.0e-12;
    assert!(
        ln_q(lower) < target_ln_k && target_ln_k < ln_q(upper),
        "synthetic in-band fixture requires a positive interior chemical root"
    );
    for _ in 0..96 {
        let midpoint = 0.5 * (lower + upper);
        if ln_q(midpoint) < target_ln_k {
            lower = midpoint;
        } else {
            upper = midpoint;
        }
    }
    0.5 * (lower + upper)
}

/// Shared synthetic P,H fixture. The candidate starts physically absent, so
/// the first fixed set is gas-only and phase control must decide its fate.
struct SyntheticLifecycleFixture {
    problem: EquilibriumProblem,
    thermochemistry: ResolvedThermochemistry,
    target_enthalpy: f64,
    independent_problem: PurePhasePhProblem,
}

fn synthetic_fixture_with_candidate_inventory(
    favorable_candidate: bool,
    initial_candidate_moles: f64,
) -> SyntheticLifecycleFixture {
    let (candidate_offset, target_candidate_moles) = if favorable_candidate {
        (0.0, Some(FAVORABLE_CANDIDATE_MOLES))
    } else {
        (20_000.0, None)
    };
    synthetic_fixture_with_candidate_offset(
        initial_candidate_moles,
        candidate_offset,
        target_candidate_moles,
    )
}

/// Builds one analytic lifecycle fixture with an explicitly controlled pure
/// candidate offset. `target_candidate_moles = None` selects the gas-only
/// boundary; otherwise the target is constructed from that positive phase
/// amount at `T_STAR`.
fn synthetic_fixture_with_candidate_offset(
    initial_candidate_moles: f64,
    candidate_offset: f64,
    target_candidate_moles: Option<f64>,
) -> SyntheticLifecycleFixture {
    assert!(
        (0.0..1.0).contains(&initial_candidate_moles),
        "the synthetic condensed inventory must preserve positive gas amounts"
    );
    let candidate_gibbs =
        move |temperature: f64| favorable_candidate_gibbs(temperature) + candidate_offset;
    let candidate_enthalpy =
        move |temperature: f64| favorable_candidate_enthalpy(temperature) + candidate_offset;
    let target_enthalpy = target_candidate_moles.map_or_else(
        || {
            // The gas-only state at T=800 K; a stable inactive candidate must
            // not be needed to satisfy the P,H constraint.
            4.0 * GAS_HEAT_CAPACITY * T_STAR
        },
        |candidate_moles| {
            // Construct an active reference point near the physical boundary.
            // It is still strongly TPD-favorable, while being a realistic
            // target for the lifecycle's deliberate trace-mole restart seed.
            (2.0 - 2.0 * candidate_moles) * GAS_HEAT_CAPACITY * T_STAR
                + (2.0 - candidate_moles) * GAS_HEAT_CAPACITY * T_STAR
                + candidate_moles * favorable_candidate_enthalpy(T_STAR)
                + candidate_moles * candidate_offset
        },
    );
    // Keep elemental totals fixed at A=2 and B=2 when an initially active
    // condensed phase is requested: A(g) + 2 S and B(g) + S remain constant.
    let initial_moles = vec![
        2.0 - 2.0 * initial_candidate_moles,
        2.0 - initial_candidate_moles,
        initial_candidate_moles,
    ];
    let numerical_seed = LogMolesInitialGuess::from_moles(
        &[
            initial_moles[0],
            initial_moles[1],
            initial_moles[2].max(1.0e-30),
        ],
        1.0e-30,
    )
    .unwrap();
    let components: Vec<String> = vec!["A(g)".into(), "B(g)".into(), "S(cond)".into()];
    let problem = EquilibriumProblem::new(
        components,
        initial_moles.clone(),
        numerical_seed,
        DMatrix::from_row_slice(3, 2, &[1.0, 0.0, 0.0, 1.0, 2.0, 1.0]),
        vec![
            Rc::new(gas_standard_gibbs) as GibbsFn,
            Rc::new(gas_standard_gibbs) as GibbsFn,
            Rc::new(candidate_gibbs) as GibbsFn,
        ],
        vec![
            Phase {
                kind: PhaseActivityModel::IdealGas,
                species: vec![0, 1],
            },
            Phase {
                kind: PhaseActivityModel::IdealSolution,
                species: vec![2],
            },
        ],
        EquilibriumConditions::new(T_STAR, 101_325.0, 101_325.0).unwrap(),
    )
    .unwrap();
    let gas = PhaseId::new(Some("gas".into()));
    let condensed = PhaseId::new(Some("condensed".into()));
    let thermochemistry = ResolvedThermochemistry::from_functions(
        vec![
            ThermochemistryProvenance::new(
                PhaseComponentId::new(gas.clone(), "A(g)"),
                "synthetic",
                "A",
                "gas",
            ),
            ThermochemistryProvenance::new(
                PhaseComponentId::new(gas, "B(g)"),
                "synthetic",
                "B",
                "gas",
            ),
            ThermochemistryProvenance::new(
                PhaseComponentId::new(condensed, "S(cond)"),
                "synthetic",
                "S",
                "condensed",
            ),
        ],
        TemperatureBounds::new(TEMPERATURE_LOWER, TEMPERATURE_UPPER).unwrap(),
        vec![
            thermo_function(|temperature| Ok(gas_standard_gibbs(temperature))),
            thermo_function(|temperature| Ok(gas_standard_gibbs(temperature))),
            thermo_function(move |temperature| Ok(candidate_gibbs(temperature))),
        ],
        vec![
            thermo_function(|temperature| Ok(GAS_HEAT_CAPACITY * temperature)),
            thermo_function(|temperature| Ok(GAS_HEAT_CAPACITY * temperature)),
            thermo_function(move |temperature| Ok(candidate_enthalpy(temperature))),
        ],
        vec![
            Some(thermo_function(|_| Ok(GAS_HEAT_CAPACITY))),
            Some(thermo_function(|_| Ok(GAS_HEAT_CAPACITY))),
            Some(thermo_function(|_| Ok(3.0 * GAS_HEAT_CAPACITY))),
        ],
    )
    .unwrap();
    let independent_problem = PurePhasePhProblem::new(
        vec!["A(g)".into(), "B(g)".into()],
        initial_moles[..2].to_vec(),
        vec![-2.0, -1.0],
        initial_candidate_moles,
        1.0,
        "S(cond)",
        vec![
            thermo_function(|temperature| Ok(gas_standard_gibbs(temperature))),
            thermo_function(|temperature| Ok(gas_standard_gibbs(temperature))),
        ],
        thermo_function(move |temperature| Ok(candidate_gibbs(temperature))),
        vec![
            thermo_function(|temperature| Ok(GAS_HEAT_CAPACITY * temperature)),
            thermo_function(|temperature| Ok(GAS_HEAT_CAPACITY * temperature)),
        ],
        thermo_function(move |temperature| Ok(candidate_enthalpy(temperature))),
        PurePhasePhConditions::new(
            101_325.0,
            101_325.0,
            target_enthalpy,
            TEMPERATURE_LOWER,
            TEMPERATURE_UPPER,
        )
        .unwrap(),
    )
    .unwrap();
    SyntheticLifecycleFixture {
        problem,
        thermochemistry,
        target_enthalpy,
        independent_problem,
    }
}

fn synthetic_fixture(favorable_candidate: bool) -> SyntheticLifecycleFixture {
    synthetic_fixture_with_candidate_inventory(favorable_candidate, 0.0)
}

fn solve_phase_controlled_ph(
    fixture: SyntheticLifecycleFixture,
) -> (
    crate::Thermodynamics::ChemEquilibrium::prepared_phase_control_runner::PreparedPhaseControlOutcome,
    ResolvedThermochemistry,
    PurePhasePhProblem,
){
    let initial_enthalpies = fixture.thermochemistry.evaluate_enthalpy(T_STAR).unwrap();
    let scale = EnthalpyScale::from_magnitudes(
        fixture.target_enthalpy,
        fixture.problem.initial_moles(),
        &initial_enthalpies,
    )
    .unwrap();
    let mut runner = PreparedPhaseControlRunner::new(fixture.problem, Vec::new(), false).unwrap();
    let settings: &mut EquilibriumSolverSettings = runner.configure_solver();
    settings.solver = Solvers::LM;
    settings.solver_policy = Some(SolverPolicy::legacy_default(Solvers::LM));
    settings.solver_params.tol = 1.0e-10;
    settings.solver_params.max_iter = 300;
    runner
        .configure_phase_control()
        .set_explicit_hysteresis(-10.0, 10.0);
    let thermochemistry = fixture.thermochemistry;
    let independent_problem = fixture.independent_problem;
    let mut temperature_seed = T_STAR;
    let bounds = TemperatureBounds::new(TEMPERATURE_LOWER, TEMPERATURE_UPPER).unwrap();
    let acceptance = PhMonolithicOptions::default();
    let outcome = runner
        .solve_with_fixed_active_solver(|runner, active, seed, species_phase, totals| {
            solve_monolithic_active_set_candidate(
                runner,
                active,
                seed,
                species_phase,
                totals,
                &thermochemistry,
                bounds,
                fixture.target_enthalpy,
                scale,
                &mut temperature_seed,
                &acceptance,
            )
        })
        .unwrap();
    (outcome, thermochemistry, independent_problem)
}

fn configure_lifecycle_runner(problem: EquilibriumProblem) -> PreparedPhaseControlRunner {
    let mut runner = PreparedPhaseControlRunner::new(problem, Vec::new(), false).unwrap();
    let settings: &mut EquilibriumSolverSettings = runner.configure_solver();
    settings.solver = Solvers::LM;
    settings.solver_policy = Some(SolverPolicy::legacy_default(Solvers::LM));
    settings.solver_params.tol = 1.0e-10;
    settings.solver_params.max_iter = 300;
    runner
        .configure_phase_control()
        .set_explicit_hysteresis(-10.0, 10.0);
    runner
}

fn assert_matches_independent_i2(
    outcome: &crate::Thermodynamics::ChemEquilibrium::prepared_phase_control_runner::PreparedPhaseControlOutcome,
    thermochemistry: &ResolvedThermochemistry,
    independent_problem: &PurePhasePhProblem,
) {
    let independent = PurePhasePhValidator::default()
        .solve(independent_problem)
        .expect("I2 scalar reference must solve the accepted P,H point");
    let temperature = outcome.solution.conditions().temperature();
    let moles = outcome.solution.moles();
    let total_enthalpy = thermochemistry
        .enthalpy_model()
        .evaluate_total(moles, temperature)
        .unwrap();
    assert!((temperature - independent.temperature).abs() <= 1.0e-4);
    assert!((moles[0] - independent.gas_moles[0]).abs() <= 1.0e-5);
    assert!((moles[1] - independent.gas_moles[1]).abs() <= 1.0e-5);
    assert!((moles[2] - independent.candidate_moles).abs() <= 1.0e-5);
    assert!((total_enthalpy - independent.total_enthalpy).abs() <= 1.0e-3);
}

fn install_accepted_continuation(
    runner: &mut PreparedPhaseControlRunner,
    outcome: &crate::Thermodynamics::ChemEquilibrium::prepared_phase_control_runner::PreparedPhaseControlOutcome,
) {
    // `PreparedPhaseControlRunner` deliberately does not invent the next
    // request's conditions. The owning P,H range workflow publishes this
    // exact pair only after acceptance; mirror that narrow ownership boundary
    // here instead of making the runner a second range orchestrator.
    let seed = LogMolesInitialGuess::from_moles(outcome.solution.moles(), 1.0e-30).unwrap();
    let gibbs = runner.prepared_problem().problem().gibbs().to_vec();
    runner
        .retarget_numeric(outcome.solution.conditions(), seed, gibbs)
        .unwrap();
    runner
        .set_continuation_phase_set(outcome.phase_control_report.final_phase_set.clone())
        .unwrap();
}

#[test]
fn i3_ph_favorable_inactive_pure_phase_is_activated_by_tpd_and_matches_i2() {
    let (outcome, thermochemistry, independent_problem) =
        solve_phase_controlled_ph(synthetic_fixture(true));
    let report = &outcome.phase_control_report;

    let gas = PhaseIndex::new(0, 2).unwrap();
    let condensed = PhaseIndex::new(1, 2).unwrap();
    assert_eq!(report.initial_phase_set.status(gas), PhaseStatus::Active);
    assert_eq!(
        report.initial_phase_set.status(condensed),
        PhaseStatus::Inactive
    );
    assert_eq!(report.final_phase_set.status(gas), PhaseStatus::Active);
    assert_eq!(
        report.final_phase_set.status(condensed),
        PhaseStatus::Active
    );
    let transition = report
        .transitions
        .iter()
        .find(|transition| !transition.activated.is_empty())
        .expect("a favorable inactive pure phase must retain activation evidence");
    assert!(matches!(
        transition.reason,
        PhaseTransitionReason::UnstableInactivePhase { minimum_tpd } if minimum_tpd < -10.0
    ));
    assert_eq!(
        report.transitions.len(),
        1,
        "the deliberately simple favorable fixture must restart exactly once"
    );
    assert!(transition.transition_duration > std::time::Duration::ZERO);
    assert!(outcome.solution.moles()[2] > 1.0e-3);

    assert_matches_independent_i2(&outcome, &thermochemistry, &independent_problem);
}

#[test]
fn i3_ph_stable_inactive_pure_phase_stays_inactive_without_transition() {
    let (outcome, _thermochemistry, independent_problem) =
        solve_phase_controlled_ph(synthetic_fixture(false));
    let report = &outcome.phase_control_report;

    let gas = PhaseIndex::new(0, 2).unwrap();
    let condensed = PhaseIndex::new(1, 2).unwrap();
    assert_eq!(report.initial_phase_set.status(gas), PhaseStatus::Active);
    assert_eq!(
        report.initial_phase_set.status(condensed),
        PhaseStatus::Inactive
    );
    assert_eq!(report.final_phase_set.status(gas), PhaseStatus::Active);
    assert_eq!(
        report.final_phase_set.status(condensed),
        PhaseStatus::Inactive
    );
    assert!(report.transitions.is_empty());
    assert!(outcome.solution.moles()[2] <= 1.0e-25);

    // I1 <-> I3 stability comparison for an absent pure phase.  The scalar
    // route defines `ln(Q)-ln(K)`, while the production TPD has energy units;
    // for this one-reaction, ideal-gas/pure-condensed fixture they must carry
    // the same sign and satisfy `TPD = R*T*(ln(Q)-ln(K))`.
    let temperature = outcome.solution.conditions().temperature();
    let independent_log_driving_force = independent_problem
        .chemical_log_residual_for_gas_moles(&outcome.solution.moles()[..2], temperature)
        .unwrap();
    let canonical_tpd = outcome.acceptance_report.phase_stability[1]
        .minimum_tpd
        .expect("the inactive pure candidate must retain TPD evidence");
    assert!(independent_log_driving_force > 0.0);
    assert!(canonical_tpd > 0.0);
    assert!(
        (canonical_tpd - MOLAR_GAS_CONSTANT * temperature * independent_log_driving_force).abs()
            <= 1.0e-6,
        "independent log driving force and canonical TPD diverged"
    );
}

#[test]
fn i3_ph_stable_active_pure_phase_recovers_the_gas_only_boundary() {
    // I3: the condensed candidate is initially present but deliberately
    // unfavourable.  The accepted state is the gas-only P,H boundary, so the
    // runner must publish an explicit phase removal rather than retain a
    // numerical trace merely because log-moles cannot represent zero.
    let (outcome, _thermochemistry, independent_problem) =
        solve_phase_controlled_ph(synthetic_fixture_with_candidate_inventory(false, 1.0e-1));
    let report = &outcome.phase_control_report;

    let gas = PhaseIndex::new(0, 2).unwrap();
    let condensed = PhaseIndex::new(1, 2).unwrap();
    assert_eq!(report.initial_phase_set.status(gas), PhaseStatus::Active);
    assert_eq!(
        report.initial_phase_set.status(condensed),
        PhaseStatus::Active
    );
    assert_eq!(report.final_phase_set.status(gas), PhaseStatus::Active);
    assert_eq!(
        report.final_phase_set.status(condensed),
        PhaseStatus::Inactive
    );
    let transition = report
        .transitions
        .iter()
        .find(|transition| !transition.deactivated.is_empty())
        .expect("stable active pure phase must retain boundary-removal evidence");
    assert!(matches!(
        transition.reason,
        // This fixture has no positive interior equilibrium for the strongly
        // unfavourable pure phase. Its physical route is therefore boundary
        // recovery, not numerical vanishing of an otherwise valid interior
        // active solution. Other production P,H states may legitimately use
        // `VanishingUnstableActivePhase`; this test does not narrow that API.
        PhaseTransitionReason::BoundaryUnstableActivePhase { .. }
    ));
    assert_eq!(
        report.transitions.len(),
        1,
        "the controlled boundary fixture must publish exactly one deactivation"
    );
    assert!(transition.transition_duration > std::time::Duration::ZERO);
    assert!(outcome.solution.moles()[2] <= 1.0e-20);
    assert!(outcome.acceptance_report.complementarity.satisfied);

    // I1 <-> I3 boundary evidence: after removing the phase, the reduced gas
    // state must make recreating it thermodynamically unfavourable. The scalar
    // coordinate is normalized by candidate stoichiometry, hence the explicit
    // division even though this fixture happens to use nu_candidate = 1.
    let temperature = outcome.solution.conditions().temperature();
    let independent_log_driving_force = independent_problem
        .chemical_log_residual_for_gas_moles(&outcome.solution.moles()[..2], temperature)
        .unwrap();
    let canonical_tpd = outcome.acceptance_report.phase_stability[1]
        .minimum_tpd
        .expect("boundary-recovered inactive candidate must retain TPD evidence");
    let expected_tpd = MOLAR_GAS_CONSTANT * temperature * independent_log_driving_force
        / independent_problem.candidate_stoichiometry();
    assert!(independent_log_driving_force > 0.0);
    assert!(canonical_tpd > 0.0);
    assert!(
        (canonical_tpd - expected_tpd).abs() <= 1.0e-6,
        "reduced-boundary independent driving force and canonical TPD diverged"
    );
}

#[test]
fn i3_ph_continuation_uses_only_accepted_states_and_matches_i2_at_each_target() {
    // I3: neighbouring enthalpy targets share one P,H lifecycle runner.  The
    // deliberately impossible middle request must roll back, after which the
    // final accepted point still starts from the preceding accepted phase set
    // and agrees with the independent I2 scalar route.
    let fixture = synthetic_fixture(true);
    let initial_enthalpies = fixture.thermochemistry.evaluate_enthalpy(T_STAR).unwrap();
    let mut runner = configure_lifecycle_runner(fixture.problem);
    let thermochemistry = fixture.thermochemistry;
    let independent_problem = fixture.independent_problem;
    let bounds = TemperatureBounds::new(TEMPERATURE_LOWER, TEMPERATURE_UPPER).unwrap();
    let acceptance = PhMonolithicOptions::default();
    let mut temperature_seed = T_STAR;
    let targets = [
        independent_problem.conditions().target_enthalpy() - 40.0,
        independent_problem.conditions().target_enthalpy() + 40.0,
    ];

    let mut accepted = Vec::new();
    for target_enthalpy in targets {
        let scale = EnthalpyScale::from_magnitudes(
            target_enthalpy,
            runner.prepared_problem().problem().initial_moles(),
            &initial_enthalpies,
        )
        .unwrap();
        let outcome = runner
            .solve_with_fixed_active_solver(|runner, active, seed, species_phase, totals| {
                solve_monolithic_active_set_candidate(
                    runner,
                    active,
                    seed,
                    species_phase,
                    totals,
                    &thermochemistry,
                    bounds,
                    target_enthalpy,
                    scale,
                    &mut temperature_seed,
                    &acceptance,
                )
            })
            .unwrap();
        let reference = independent_problem
            .with_target_enthalpy(target_enthalpy)
            .unwrap();
        assert_matches_independent_i2(&outcome, &thermochemistry, &reference);
        install_accepted_continuation(&mut runner, &outcome);
        accepted.push(outcome);
    }
    assert_eq!(
        accepted[0]
            .phase_control_report
            .final_phase_set
            .active_mask(),
        accepted[1]
            .phase_control_report
            .initial_phase_set
            .active_mask(),
        "the second target must start from the preceding accepted phase state"
    );

    let rejected_target = independent_problem.conditions().target_enthalpy() - 1.0e6;
    let rejected_scale = EnthalpyScale::from_magnitudes(
        rejected_target,
        runner.prepared_problem().problem().initial_moles(),
        &initial_enthalpies,
    )
    .unwrap();
    assert!(
        runner
            .solve_with_fixed_active_solver(|runner, active, seed, species_phase, totals| {
                solve_monolithic_active_set_candidate(
                    runner,
                    active,
                    seed,
                    species_phase,
                    totals,
                    &thermochemistry,
                    bounds,
                    rejected_target,
                    rejected_scale,
                    &mut temperature_seed,
                    &acceptance,
                )
            })
            .is_err()
    );

    let final_target = independent_problem.conditions().target_enthalpy();
    let final_scale = EnthalpyScale::from_magnitudes(
        final_target,
        runner.prepared_problem().problem().initial_moles(),
        &initial_enthalpies,
    )
    .unwrap();
    let recovered = runner
        .solve_with_fixed_active_solver(|runner, active, seed, species_phase, totals| {
            solve_monolithic_active_set_candidate(
                runner,
                active,
                seed,
                species_phase,
                totals,
                &thermochemistry,
                bounds,
                final_target,
                final_scale,
                &mut temperature_seed,
                &acceptance,
            )
        })
        .unwrap();
    assert_eq!(
        accepted[1]
            .phase_control_report
            .final_phase_set
            .active_mask(),
        recovered
            .phase_control_report
            .initial_phase_set
            .active_mask(),
        "a failed target must not replace the accepted continuation phase set"
    );
    assert_matches_independent_i2(&recovered, &thermochemistry, &independent_problem);
}

#[test]
fn i3_ph_in_band_hysteresis_keeps_only_the_accepted_active_history() {
    // I3: this offset makes the gas-only candidate's TPD fall inside the
    // explicit [-10, +10] J/mol band. The same target must therefore keep an
    // already accepted active phase, but cannot create that phase from a
    // fresh gas-only history.
    const IN_BAND_CANDIDATE_OFFSET: f64 = 12.0;
    let in_band_candidate_moles = candidate_equilibrium_moles_at_t_star(IN_BAND_CANDIDATE_OFFSET);
    let in_band_reference = synthetic_fixture_with_candidate_offset(
        0.0,
        IN_BAND_CANDIDATE_OFFSET,
        Some(in_band_candidate_moles),
    );
    let target_enthalpy = in_band_reference.target_enthalpy;
    let in_band_thermochemistry = in_band_reference.thermochemistry;
    let active_fixture = synthetic_fixture_with_candidate_offset(
        in_band_candidate_moles,
        0.0,
        Some(FAVORABLE_CANDIDATE_MOLES),
    );
    let activation_target_enthalpy = active_fixture.target_enthalpy;
    let initial_enthalpies = active_fixture
        .thermochemistry
        .evaluate_enthalpy(T_STAR)
        .unwrap();
    let activation_scale = EnthalpyScale::from_magnitudes(
        activation_target_enthalpy,
        active_fixture.problem.initial_moles(),
        &initial_enthalpies,
    )
    .unwrap();
    let mut runner = configure_lifecycle_runner(active_fixture.problem);
    let activation_thermochemistry = active_fixture.thermochemistry;
    let bounds = TemperatureBounds::new(TEMPERATURE_LOWER, TEMPERATURE_UPPER).unwrap();
    let acceptance = PhMonolithicOptions::default();
    let mut temperature_seed = T_STAR;

    let established = runner
        .solve_with_fixed_active_solver(|runner, active, seed, species_phase, totals| {
            solve_monolithic_active_set_candidate(
                runner,
                active,
                seed,
                species_phase,
                totals,
                &activation_thermochemistry,
                bounds,
                activation_target_enthalpy,
                activation_scale,
                &mut temperature_seed,
                &acceptance,
            )
        })
        .unwrap();
    assert_eq!(
        established
            .phase_control_report
            .final_phase_set
            .active_mask(),
        &[true, true]
    );
    install_accepted_continuation(&mut runner, &established);

    let in_band_enthalpies = in_band_thermochemistry.evaluate_enthalpy(T_STAR).unwrap();
    let in_band_scale = EnthalpyScale::from_magnitudes(
        target_enthalpy,
        runner.prepared_problem().problem().initial_moles(),
        &in_band_enthalpies,
    )
    .unwrap();

    let continued = runner
        .solve_with_fixed_active_solver(|runner, active, seed, species_phase, totals| {
            solve_monolithic_active_set_candidate(
                runner,
                active,
                seed,
                species_phase,
                totals,
                &in_band_thermochemistry,
                bounds,
                target_enthalpy,
                in_band_scale,
                &mut temperature_seed,
                &acceptance,
            )
        })
        .unwrap();
    assert_eq!(
        continued
            .phase_control_report
            .initial_phase_set
            .active_mask(),
        &[true, true]
    );
    assert_eq!(
        continued.phase_control_report.final_phase_set.active_mask(),
        &[true, true],
        "the in-band keep threshold must retain the accepted active phase"
    );
    assert!(continued.phase_control_report.transitions.is_empty());

    let fresh_fixture = synthetic_fixture_with_candidate_offset(
        0.0,
        IN_BAND_CANDIDATE_OFFSET,
        Some(in_band_candidate_moles),
    );
    assert!((fresh_fixture.target_enthalpy - target_enthalpy).abs() <= f64::EPSILON);
    let (fresh, _fresh_thermochemistry, _fresh_independent_problem) =
        solve_phase_controlled_ph(fresh_fixture);
    assert_eq!(
        fresh.phase_control_report.initial_phase_set.active_mask(),
        &[true, false]
    );
    assert_eq!(
        fresh.phase_control_report.final_phase_set.active_mask(),
        &[true, false],
        "an in-band gas-only history must not activate the candidate"
    );
    assert!(fresh.phase_control_report.transitions.is_empty());
    let fresh_tpd = fresh.acceptance_report.phase_stability[1]
        .minimum_tpd
        .expect("the fresh inactive candidate must retain TPD evidence");
    assert!(
        (-10.0..10.0).contains(&fresh_tpd),
        "fresh candidate TPD {fresh_tpd:e} must lie inside the explicit hysteresis band"
    );
}
