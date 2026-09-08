//! Layered synthetic regression tests for pure-phase boundary validation.
//!
//! This module deliberately keeps three questions separate:
//!
//! 1. Does independent `ln(Q)-ln(K)` agree with canonical TPD?
//! 2. Does a fixed two-phase canonical Gibbs solve agree with the independent
//!    scalar `K_eq` extent result?
//! 3. Does the production phase-control lifecycle make the history-dependent
//!    topology decision correctly?
//!
//! The fixtures are abstract and synthetic. Real thermochemical data, P,H,
//! external databases, and multi-component candidate phases are intentionally
//! outside this regression layer.

use crate::Thermodynamics::ChemEquilibrium::equilibrium_activity::PhaseActivityModel;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_ids::PhaseIndex;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_log_moles::{GibbsFn, Phase, Solvers};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_nonlinear::ReactionExtentError;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_phase_stability::{
    PhaseStabilityConditions, PhaseStabilityLayout, PhaseStabilityReport, PhaseStabilityStatus,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_prepared_runner::PreparedEquilibriumRunner;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_problem::{
    EquilibriumConditions, EquilibriumProblem, LogMolesInitialGuess, PreparedEquilibriumProblem,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_solver_policy::SolverPolicy;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_workflows::{
    InitialPhaseSet, PHASE_CONTROL_TRACE_MOLE_FLOOR, PhaseManager, PhaseSet, PhaseStabilityReason,
    PhaseTransitionPlan, PhaseTransitionReason, compute_phase_stability_reports,
};
use crate::Thermodynamics::ChemEquilibrium::phase_boundary_validation::{
    CanonicalPurePhaseEvidence, MOLAR_GAS_CONSTANT, PurePhaseBoundaryElementComposition,
    PurePhaseBoundaryPrediction, PurePhaseBoundaryProblem, PurePhaseBoundaryStructuralTolerances,
    PurePhaseBoundaryTemperatureSearchSettings, PurePhaseBoundaryTolerances,
    PurePhaseBoundaryValidator, PurePhaseCrossValidationStatus, PurePhaseCrossValidationTolerances,
    PurePhaseTopologyExpectation, bisect_pure_phase_boundary_temperature,
    compare_pure_phase_validation, cross_validate_pure_phase, evaluate_pure_phase_boundary,
};
use crate::Thermodynamics::ChemEquilibrium::prepared_phase_control_runner::{
    PreparedPhaseControlOutcome, PreparedPhaseControlRunner,
};
use nalgebra::DMatrix;
use std::rc::Rc;

fn constant_gibbs(value: f64) -> GibbsFn {
    Rc::new(move |_| value)
}

fn analytic_boundary_problem(
    temperature: f64,
    boundary_temperature: f64,
    slope_per_kelvin: f64,
) -> PurePhaseBoundaryProblem {
    let pressure = 101_325.0;
    let ln_q_at_absence = 2.0_f64.ln();
    let candidate_gibbs: GibbsFn = Rc::new(move |evaluated_temperature| {
        let ln_k =
            ln_q_at_absence + slope_per_kelvin * (boundary_temperature - evaluated_temperature);
        -MOLAR_GAS_CONSTANT * evaluated_temperature * ln_k
    });
    PurePhaseBoundaryProblem::new(
        vec!["A".to_string(), "B".to_string()],
        vec![1.0, 1.0],
        vec![-2.0, 1.0],
        1.0,
        vec![constant_gibbs(0.0), constant_gibbs(0.0)],
        candidate_gibbs,
        EquilibriumConditions::new(temperature, pressure, pressure).unwrap(),
        "S",
    )
    .unwrap()
}

fn canonical_pure_candidate_tpd(
    temperature: f64,
    boundary_temperature: f64,
    slope_per_kelvin: f64,
) -> f64 {
    canonical_pure_candidate_tpd_at_conditions(
        temperature,
        boundary_temperature,
        slope_per_kelvin,
        101_325.0,
        101_325.0,
    )
}

fn canonical_pure_candidate_tpd_at_conditions(
    temperature: f64,
    boundary_temperature: f64,
    slope_per_kelvin: f64,
    pressure: f64,
    reference_pressure: f64,
) -> f64 {
    let ln_q_at_absence = 2.0_f64.ln();
    let ln_k = ln_q_at_absence + slope_per_kelvin * (boundary_temperature - temperature);
    canonical_tpd_for_target_k_at_conditions(ln_k.exp(), temperature, pressure, reference_pressure)
}

fn canonical_tpd_for_target_k_at_conditions(
    target_k: f64,
    temperature: f64,
    pressure: f64,
    reference_pressure: f64,
) -> f64 {
    let candidate_gibbs = -MOLAR_GAS_CONSTANT * temperature * target_k.ln();
    let phases = vec![
        Phase {
            kind: PhaseActivityModel::IdealGas,
            species: vec![0, 1],
        },
        Phase {
            kind: PhaseActivityModel::IdealSolution,
            species: vec![2],
        },
    ];
    let phase_set =
        PhaseSet::from_policy(&InitialPhaseSet::FromInitialMoles, &[true, false]).unwrap();
    let reports = compute_phase_stability_reports(
        &[0.0, 0.0, PHASE_CONTROL_TRACE_MOLE_FLOOR.ln()],
        &[
            constant_gibbs(0.0),
            constant_gibbs(0.0),
            constant_gibbs(candidate_gibbs),
        ],
        &phases,
        &[0, 0, 1],
        &DMatrix::from_row_slice(3, 2, &[1.0, 1.0, 1.0, 0.0, 1.0, 2.0]),
        temperature,
        pressure,
        reference_pressure,
        &phase_set,
    )
    .expect("canonical TPD accepts the synthetic gas-plus-pure layout");
    reports[1]
        .minimum_tpd
        .expect("inactive pure candidate has an evaluated TPD minimum")
}

/// Test-local bisection over canonical TPD.  It stays in the regression module
/// because production phase control uses hysteresis rather than a mathematical
/// root finder to select a discrete topology.
fn bisect_canonical_tpd_temperature(
    lower_temperature: f64,
    upper_temperature: f64,
    boundary_temperature: f64,
    slope_per_kelvin: f64,
    max_abs_tpd: f64,
    max_iterations: usize,
) -> Result<f64, String> {
    let mut lower = lower_temperature;
    let mut upper = upper_temperature;
    let mut lower_tpd = canonical_pure_candidate_tpd(lower, boundary_temperature, slope_per_kelvin);
    let upper_tpd = canonical_pure_candidate_tpd(upper, boundary_temperature, slope_per_kelvin);
    if lower_tpd.signum() == upper_tpd.signum() {
        return Err(format!(
            "canonical TPD boundary root is not bracketed: lower={lower:.16e}, upper={upper:.16e}, tpd_lower={lower_tpd:.16e}, tpd_upper={upper_tpd:.16e}"
        ));
    }

    for _ in 0..max_iterations {
        let midpoint = 0.5 * (lower + upper);
        let midpoint_tpd =
            canonical_pure_candidate_tpd(midpoint, boundary_temperature, slope_per_kelvin);
        if midpoint_tpd.abs() <= max_abs_tpd {
            return Ok(midpoint);
        }
        if lower_tpd.signum() != midpoint_tpd.signum() {
            upper = midpoint;
        } else {
            lower = midpoint;
            lower_tpd = midpoint_tpd;
        }
    }
    let upper_tpd = canonical_pure_candidate_tpd(upper, boundary_temperature, slope_per_kelvin);
    Err(format!(
        "canonical TPD boundary bisection did not converge: lower={lower:.16e}, upper={upper:.16e}, tpd_lower={lower_tpd:.16e}, tpd_upper={upper_tpd:.16e}, max_abs_tpd={max_abs_tpd:.16e}, max_iterations={max_iterations}"
    ))
}

fn boundary_problem_for_target_k(target_k: f64, scale: f64) -> PurePhaseBoundaryProblem {
    boundary_problem_for_target_k_at_temperature(target_k, scale, 1_000.0)
}

fn boundary_problem_for_target_k_at_temperature(
    target_k: f64,
    scale: f64,
    temperature: f64,
) -> PurePhaseBoundaryProblem {
    boundary_problem_for_target_k_at_conditions(target_k, scale, temperature, 101_325.0, 101_325.0)
}

fn boundary_problem_for_target_k_at_conditions(
    target_k: f64,
    scale: f64,
    temperature: f64,
    pressure: f64,
    reference_pressure: f64,
) -> PurePhaseBoundaryProblem {
    PurePhaseBoundaryProblem::new(
        vec!["A".to_string(), "B".to_string()],
        vec![1.0, 1.0],
        vec![-2.0 * scale, scale],
        scale,
        vec![constant_gibbs(0.0), constant_gibbs(0.0)],
        constant_gibbs(-MOLAR_GAS_CONSTANT * temperature * target_k.ln()),
        EquilibriumConditions::new(temperature, pressure, reference_pressure).unwrap(),
        "S",
    )
    .unwrap()
}

fn scaled_boundary_problem(scale: f64) -> PurePhaseBoundaryProblem {
    boundary_problem_for_target_k(4.0, scale)
}

/// Builds the same A/B/S equilibrium with an inert gas I carrying its own
/// conserved element.  I has no reaction stoichiometry, but it must still
/// affect the gas activities through the total gas mole number.
fn inert_dilution_boundary_problem(inert_moles: f64) -> PurePhaseBoundaryProblem {
    let temperature = 1_000.0;
    let pressure = 101_325.0;
    PurePhaseBoundaryProblem::new(
        vec!["A".to_string(), "B".to_string(), "I".to_string()],
        vec![1.0, 1.0, inert_moles],
        vec![-2.0, 1.0, 0.0],
        1.0,
        vec![
            constant_gibbs(0.0),
            constant_gibbs(0.0),
            constant_gibbs(0.0),
        ],
        constant_gibbs(-MOLAR_GAS_CONSTANT * temperature * 4.0_f64.ln()),
        EquilibriumConditions::new(temperature, pressure, pressure).unwrap(),
        "S",
    )
    .unwrap()
    .with_element_composition(
        PurePhaseBoundaryElementComposition::new(
            vec!["X".to_string(), "Y".to_string(), "Inert".to_string()],
            DMatrix::from_row_slice(
                3,
                3,
                &[
                    1.0, 1.0, 0.0, // A
                    2.0, 0.0, 0.0, // B
                    0.0, 0.0, 1.0, // I
                ],
            ),
            vec![0.0, 2.0, 0.0],
        )
        .unwrap(),
        PurePhaseBoundaryStructuralTolerances::default(),
    )
    .unwrap()
}

/// Mirrors the canonical inactive-pure-phase TPD calculation for the inert
/// dilution fixture.  It verifies that the independent reaction expression
/// and the production TPD convention react identically to dilution.
fn canonical_inert_dilution_tpd(inert_moles: f64) -> f64 {
    let temperature = 1_000.0;
    let pressure = 101_325.0;
    let phases = vec![
        Phase {
            kind: PhaseActivityModel::IdealGas,
            species: vec![0, 1, 2],
        },
        Phase {
            kind: PhaseActivityModel::IdealSolution,
            species: vec![3],
        },
    ];
    let reports = compute_phase_stability_reports(
        &[
            0.0,
            0.0,
            inert_moles.ln(),
            PHASE_CONTROL_TRACE_MOLE_FLOOR.ln(),
        ],
        &[
            constant_gibbs(0.0),
            constant_gibbs(0.0),
            constant_gibbs(0.0),
            constant_gibbs(-MOLAR_GAS_CONSTANT * temperature * 4.0_f64.ln()),
        ],
        &phases,
        &[0, 0, 0, 1],
        &DMatrix::from_row_slice(
            4,
            3,
            &[
                1.0, 1.0, 0.0, // A
                2.0, 0.0, 0.0, // B
                0.0, 0.0, 1.0, // I
                0.0, 2.0, 0.0, // S
            ],
        ),
        temperature,
        pressure,
        pressure,
        &phase_set(&[true, false]),
    )
    .expect("canonical TPD accepts the inert-dilution layout");
    reports[1]
        .minimum_tpd
        .expect("inactive pure candidate has an evaluated TPD minimum")
}

fn permuted_boundary_problem_for_target_k(target_k: f64, scale: f64) -> PurePhaseBoundaryProblem {
    let temperature = 1_000.0;
    let pressure = 101_325.0;
    PurePhaseBoundaryProblem::new(
        vec!["B".to_string(), "A".to_string()],
        vec![1.0, 1.0],
        vec![scale, -2.0 * scale],
        scale,
        vec![constant_gibbs(0.0), constant_gibbs(0.0)],
        constant_gibbs(-MOLAR_GAS_CONSTANT * temperature * target_k.ln()),
        EquilibriumConditions::new(temperature, pressure, pressure).unwrap(),
        "S",
    )
    .unwrap()
    .with_element_composition(
        // Species are B,A here and the element columns are intentionally Y,X.
        // This is physically identical to A,B with X,Y, but neither rows nor
        // columns retain their original positions.
        PurePhaseBoundaryElementComposition::new(
            vec!["Y".to_string(), "X".to_string()],
            DMatrix::from_row_slice(2, 2, &[0.0, 1.0, 1.0, 1.0]),
            vec![2.0, 1.0],
        )
        .unwrap(),
        PurePhaseBoundaryStructuralTolerances::default(),
    )
    .unwrap()
}

fn phase_set(active: &[bool]) -> PhaseSet {
    PhaseSet::from_policy(&InitialPhaseSet::FromInitialMoles, active).unwrap()
}

/// Builds the production runner for the same strict synthetic family used by
/// the independent validator:
///
///     2 A(g) <=> B(g) + S(pure)
///
/// The candidate begins physically absent.  Consequently, the first fixed
/// active-set solve contains only the gas phase; any later appearance of `S`
/// must come from the production TPD lifecycle, not from its initial seed.
fn production_phase_control_runner(target_k: f64) -> PreparedPhaseControlRunner {
    production_phase_control_runner_with_candidate_inventory(target_k, 0.0)
}

/// Creates the same production topology with a controllable initial inventory
/// of the pure phase.  A positive amount is needed only by the disappearance
/// story: it makes the phase genuinely active before the nonlinear solve
/// drives it below the phase-destruction threshold.
fn production_phase_control_runner_with_candidate_inventory(
    target_k: f64,
    candidate_initial_moles: f64,
) -> PreparedPhaseControlRunner {
    let temperature = 1_000.0;
    let pressure = 101_325.0;
    let initial_moles = vec![1.0, 1.0, candidate_initial_moles];
    let numeric_seed = LogMolesInitialGuess::from_moles(
        &[
            1.0,
            1.0,
            candidate_initial_moles.max(PHASE_CONTROL_TRACE_MOLE_FLOOR),
        ],
        PHASE_CONTROL_TRACE_MOLE_FLOOR,
    )
    .expect("trace seed must be valid even though the candidate is physically absent");
    let candidate_gibbs = -MOLAR_GAS_CONSTANT * temperature * target_k.ln();
    let problem = EquilibriumProblem::new(
        vec!["A".to_string(), "B".to_string(), "S".to_string()],
        initial_moles,
        numeric_seed,
        DMatrix::from_row_slice(3, 2, &[1.0, 1.0, 1.0, 0.0, 1.0, 2.0]),
        vec![
            constant_gibbs(0.0),
            constant_gibbs(0.0),
            constant_gibbs(candidate_gibbs),
        ],
        vec![
            Phase {
                kind: PhaseActivityModel::IdealGas,
                species: vec![0, 1],
            },
            Phase {
                // A one-component ideal solution is a pure condensed phase:
                // its phase activity is identically one.
                kind: PhaseActivityModel::IdealSolution,
                species: vec![2],
            },
        ],
        EquilibriumConditions::new(temperature, pressure, pressure)
            .expect("synthetic P,T conditions must be valid"),
    )
    .expect("strict synthetic topology must produce a valid canonical problem");

    let mut runner = PreparedPhaseControlRunner::new(problem, Vec::new(), false)
        .expect("production phase-control runner must prepare the synthetic problem");
    let settings = runner.configure_solver();
    settings.solver_policy = Some(SolverPolicy::legacy_default(Solvers::LM));
    // This is test-local precision for comparison with the independent scalar
    // bisection route.  It does not alter production defaults.
    settings.solver_params.tol = 1e-12;
    settings.solver_params.max_iter = 200;
    runner
}

/// Retargets the synthetic prepared runner without rebuilding its layout.
/// The candidate Gibbs closure is point-local in this fixture, while the
/// runner retains only accepted continuation state across the two points.
fn retarget_production_runner(
    runner: &mut PreparedPhaseControlRunner,
    temperature: f64,
    target_k: f64,
    seed: LogMolesInitialGuess,
) {
    let pressure = 101_325.0;
    runner
        .retarget_numeric(
            EquilibriumConditions::new(temperature, pressure, pressure)
                .expect("synthetic retarget conditions must be valid"),
            seed,
            vec![
                constant_gibbs(0.0),
                constant_gibbs(0.0),
                constant_gibbs(-MOLAR_GAS_CONSTANT * temperature * target_k.ln()),
            ],
        )
        .expect("prepared synthetic runner must accept a numeric retarget");
}

/// For the gas-only A/B state in this fixture, `ln(Q) = ln(2)`.  Choosing K
/// from a requested TPD makes the temperature sweep analytically controlled:
/// `TPD = R*T*(ln(2)-ln(K))` per mole of pure candidate.
fn target_k_for_absence_tpd(temperature: f64, tpd: f64) -> f64 {
    2.0 * (-tpd / (MOLAR_GAS_CONSTANT * temperature)).exp()
}

/// Solves the fixed two-phase topology through canonical Gibbs minimization.
///
/// This is deliberately distinct from `PreparedPhaseControlRunner`: the
/// candidate phase is supplied as active from the beginning, so the test
/// isolates fixed-topology thermodynamics from outer-loop behavior.
fn canonical_fixed_phase_solution_for_target_k(target_k: f64) -> Vec<f64> {
    let temperature = 1_000.0;
    let pressure = 101_325.0;
    let initial_moles = vec![1.0, 1.0, 1e-20];
    let candidate_gibbs = -MOLAR_GAS_CONSTANT * temperature * target_k.ln();
    let problem = EquilibriumProblem::new(
        vec!["A".to_string(), "B".to_string(), "S".to_string()],
        initial_moles.clone(),
        LogMolesInitialGuess::from_moles(&initial_moles, 1e-30)
            .expect("fixed-topology fixture must have a valid log seed"),
        DMatrix::from_row_slice(3, 2, &[1.0, 1.0, 1.0, 0.0, 1.0, 2.0]),
        vec![
            constant_gibbs(0.0),
            constant_gibbs(0.0),
            constant_gibbs(candidate_gibbs),
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
        EquilibriumConditions::new(temperature, pressure, pressure)
            .expect("fixed-topology conditions must be valid"),
    )
    .expect("strict synthetic family must build a canonical problem");
    let prepared = PreparedEquilibriumProblem::new(problem)
        .expect("canonical path must derive its reaction basis");
    let mut runner = PreparedEquilibriumRunner::new(prepared, Vec::new())
        .expect("numeric canonical runner must prepare");
    let settings = runner.configure();
    settings.solver_policy = Some(SolverPolicy::legacy_default(Solvers::LM));
    // Test-local accuracy matches the independent scalar bisection reference.
    settings.solver_params.tol = 1e-12;
    settings.solver_params.max_iter = 200;
    runner
        .solve()
        .expect("canonical fixed-topology fixture must solve")
        .solution
        .moles()
        .to_vec()
}

/// Independent non-collinear synthetic family:
///
///     A(g) + B(g) <=> C(g) + 2 S(pure)
///
/// The three gas rows are independent over abstract elements X/Y/Z, while
/// the full four-species system has exactly one reaction direction.  `S` has
/// a coefficient of two, so this fixture specifically checks that TPD is
/// compared with `Delta_rG / nu_S`, not raw reaction Gibbs energy.
fn second_family_boundary_problem(target_k: f64) -> PurePhaseBoundaryProblem {
    let temperature = 900.0;
    let pressure = 101_325.0;
    let composition = PurePhaseBoundaryElementComposition::new(
        vec!["X".to_string(), "Y".to_string(), "Z".to_string()],
        DMatrix::from_row_slice(3, 3, &[2.0, 0.0, 2.0, 0.0, 2.0, 2.0, 0.0, 0.0, 2.0]),
        vec![1.0, 1.0, 1.0],
    )
    .expect("second-family element matrix is valid");
    let candidate_gibbs = -MOLAR_GAS_CONSTANT * temperature * target_k.ln() / 2.0;
    PurePhaseBoundaryProblem::new(
        vec!["A".to_string(), "B".to_string(), "C".to_string()],
        vec![1.0, 1.0, 1.0],
        vec![-1.0, -1.0, 1.0],
        2.0,
        vec![
            constant_gibbs(0.0),
            constant_gibbs(0.0),
            constant_gibbs(0.0),
        ],
        constant_gibbs(candidate_gibbs),
        EquilibriumConditions::new(temperature, pressure, pressure)
            .expect("second-family P,T conditions are valid"),
        "S",
    )
    .expect("second-family boundary problem is valid")
    .with_element_composition(
        composition,
        PurePhaseBoundaryStructuralTolerances::default(),
    )
    .expect("second-family reaction conserves every abstract element")
}

fn second_family_canonical_tpd(target_k: f64) -> f64 {
    let problem = second_family_boundary_problem(target_k);
    let temperature = problem.conditions().temperature();
    let candidate_gibbs = -MOLAR_GAS_CONSTANT * temperature * target_k.ln() / 2.0;
    let phase_set = phase_set(&[true, false]);
    let reports = compute_phase_stability_reports(
        &[0.0, 0.0, 0.0, PHASE_CONTROL_TRACE_MOLE_FLOOR.ln()],
        &[
            constant_gibbs(0.0),
            constant_gibbs(0.0),
            constant_gibbs(0.0),
            constant_gibbs(candidate_gibbs),
        ],
        &[
            Phase {
                kind: PhaseActivityModel::IdealGas,
                species: vec![0, 1, 2],
            },
            Phase {
                kind: PhaseActivityModel::IdealSolution,
                species: vec![3],
            },
        ],
        &[0, 0, 0, 1],
        &DMatrix::from_row_slice(
            4,
            3,
            &[2.0, 0.0, 2.0, 0.0, 2.0, 2.0, 0.0, 0.0, 2.0, 1.0, 1.0, 1.0],
        ),
        temperature,
        problem.conditions().pressure(),
        problem.conditions().reference_pressure(),
        &phase_set,
    )
    .expect("second-family canonical TPD must evaluate");
    reports[1]
        .minimum_tpd
        .expect("inactive pure second-family candidate has a TPD minimum")
}

fn second_family_production_runner(target_k: f64) -> PreparedPhaseControlRunner {
    let problem = second_family_boundary_problem(target_k);
    let temperature = problem.conditions().temperature();
    let pressure = problem.conditions().pressure();
    let initial_moles = vec![1.0, 1.0, 1.0, 0.0];
    let candidate_gibbs = -MOLAR_GAS_CONSTANT * temperature * target_k.ln() / 2.0;
    let numeric_seed = LogMolesInitialGuess::from_moles(
        &[1.0, 1.0, 1.0, PHASE_CONTROL_TRACE_MOLE_FLOOR],
        PHASE_CONTROL_TRACE_MOLE_FLOOR,
    )
    .expect("second-family trace seed is valid");
    let canonical_problem = EquilibriumProblem::new(
        vec![
            "A".to_string(),
            "B".to_string(),
            "C".to_string(),
            "S".to_string(),
        ],
        initial_moles,
        numeric_seed,
        DMatrix::from_row_slice(
            4,
            3,
            &[2.0, 0.0, 2.0, 0.0, 2.0, 2.0, 0.0, 0.0, 2.0, 1.0, 1.0, 1.0],
        ),
        vec![
            constant_gibbs(0.0),
            constant_gibbs(0.0),
            constant_gibbs(0.0),
            constant_gibbs(candidate_gibbs),
        ],
        vec![
            Phase {
                kind: PhaseActivityModel::IdealGas,
                species: vec![0, 1, 2],
            },
            Phase {
                kind: PhaseActivityModel::IdealSolution,
                species: vec![3],
            },
        ],
        EquilibriumConditions::new(temperature, pressure, pressure)
            .expect("second-family canonical conditions are valid"),
    )
    .expect("second-family canonical problem is valid");
    let mut runner = PreparedPhaseControlRunner::new(canonical_problem, Vec::new(), false)
        .expect("second-family runner must prepare");
    let settings = runner.configure_solver();
    settings.solver_policy = Some(SolverPolicy::legacy_default(Solvers::LM));
    settings.solver_params.tol = 1e-12;
    settings.solver_params.max_iter = 200;
    runner
}

fn lifecycle_stability_reports(active: &[bool], candidate_tpd: f64) -> Vec<PhaseStabilityReport> {
    active
        .iter()
        .enumerate()
        .map(|(phase, &is_active)| PhaseStabilityReport {
            phase: PhaseIndex::new(phase, active.len()).unwrap(),
            active: is_active,
            status: PhaseStabilityStatus::Evaluated,
            conditions: PhaseStabilityConditions {
                temperature: 600.0,
                pressure: 101_325.0,
                reference_pressure: 101_325.0,
            },
            layout: PhaseStabilityLayout {
                system_species_count: 2,
                system_phase_count: active.len(),
                element_count: 1,
                phase_component_indices: vec![phase],
            },
            minimum_tpd: if phase == 1 {
                Some(candidate_tpd)
            } else {
                None
            },
            incipient_composition: None,
            element_potentials: None,
            elemental_feasibility: None,
            minimizer: None,
        })
        .collect()
}

#[test]
fn analytic_temperature_sweep_matches_independent_keq_and_canonical_tpd() {
    let boundary_temperature = 600.0;
    let slope_per_kelvin = 0.02;
    let tolerances = PurePhaseBoundaryTolerances {
        max_abs_boundary_log_residual: 1e-11,
        ..PurePhaseBoundaryTolerances::default()
    };

    for temperature in [500.0, boundary_temperature, 700.0] {
        let problem =
            analytic_boundary_problem(temperature, boundary_temperature, slope_per_kelvin);
        let independent = evaluate_pure_phase_boundary(&problem, tolerances).unwrap();
        let canonical_tpd =
            canonical_pure_candidate_tpd(temperature, boundary_temperature, slope_per_kelvin);
        let expected_tpd =
            independent.reaction_gibbs_at_absence / problem.candidate_stoichiometry();

        assert!(
            (canonical_tpd - expected_tpd).abs() < 1e-8,
            "T={temperature}: canonical TPD and independent driving force diverged"
        );
        if independent.log_residual_at_absence.abs() <= tolerances.max_abs_boundary_log_residual {
            assert!(
                canonical_tpd.abs() < 1e-8,
                "T={temperature} must be the boundary"
            );
            assert_eq!(
                independent.prediction,
                PurePhaseBoundaryPrediction::Boundary
            );
        } else {
            assert_eq!(
                canonical_tpd.signum(),
                independent.log_residual_at_absence.signum(),
                "T={temperature}: TPD and ln(Q)-ln(K) must have the same sign"
            );
        }
    }
}

#[test]
fn canonical_tpd_boundary_temperature_matches_independent_keq_root() {
    let known_boundary_temperature = 600.0;
    let slope_per_kelvin = 0.02;
    let log_tolerances = PurePhaseBoundaryTolerances {
        max_abs_boundary_log_residual: 1e-11,
        ..PurePhaseBoundaryTolerances::default()
    };
    let independent = bisect_pure_phase_boundary_temperature(
        500.0,
        700.0,
        PurePhaseBoundaryTemperatureSearchSettings::default(),
        log_tolerances,
        |temperature| {
            Ok(analytic_boundary_problem(
                temperature,
                known_boundary_temperature,
                slope_per_kelvin,
            ))
        },
    )
    .expect("analytic independent boundary must be bracketed");
    let canonical = bisect_canonical_tpd_temperature(
        500.0,
        700.0,
        known_boundary_temperature,
        slope_per_kelvin,
        1e-9,
        96,
    )
    .expect("test-local canonical TPD root must converge inside its bracket");

    assert!((independent.temperature - known_boundary_temperature).abs() < 1e-8);
    assert!((canonical - known_boundary_temperature).abs() < 1e-8);
    assert!(
        (canonical - independent.temperature).abs() < 1e-8,
        "canonical TPD and independent K_eq roots must be the same mathematical boundary"
    );
}

#[test]
fn canonical_tpd_boundary_bisection_rejects_unconverged_bracket() {
    let error = bisect_canonical_tpd_temperature(500.0, 700.0, 611.125, 0.02, 0.0, 1)
        .expect_err("test-local root finder must not publish an unverified midpoint");

    assert!(error.contains("did not converge"));
    assert!(error.contains("tpd_lower"));
    assert!(error.contains("tpd_upper"));
}

#[test]
fn canonical_pure_phase_tpd_matches_independent_boundary_driving_force() {
    let temperature = 500.0;
    let problem = analytic_boundary_problem(temperature, 600.0, 0.02);
    let independent =
        evaluate_pure_phase_boundary(&problem, PurePhaseBoundaryTolerances::default())
            .expect("independent boundary equation must be well defined");
    let canonical_tpd = canonical_pure_candidate_tpd(temperature, 600.0, 0.02);
    let independent_per_candidate =
        independent.reaction_gibbs_at_absence / problem.candidate_stoichiometry();

    assert!(
        (canonical_tpd - independent_per_candidate).abs() < 1e-8,
        "TPD must equal RT*(lnQ-lnK)/nu_candidate for the controlled fixture"
    );
}

#[test]
fn exact_boundary_reports_history_dependent_topology_not_a_false_mismatch() {
    let problem = analytic_boundary_problem(600.0, 600.0, 0.02);
    let independent = PurePhaseBoundaryValidator::default()
        .validate(&problem)
        .unwrap();
    let comparison = compare_pure_phase_validation(
        &problem,
        &CanonicalPurePhaseEvidence {
            gas_species: problem.gas_species().to_vec(),
            candidate_name: problem.candidate_name().to_string(),
            // Either state is permitted on an exact thermodynamic boundary.
            candidate_active: true,
            gas_moles: vec![1.0, 1.0],
            boundary_gas_moles: vec![1.0, 1.0],
            candidate_moles: 0.25,
            boundary_minimum_tpd: Some(0.0),
        },
        &independent,
        PurePhaseCrossValidationTolerances::default(),
    )
    .unwrap();

    assert_eq!(
        comparison.topology_expectation,
        PurePhaseTopologyExpectation::HysteresisDependent
    );
    assert_eq!(comparison.topology_agreement, None);
    assert_eq!(comparison.composition_agreement, None);
    assert_eq!(comparison.thermodynamic_agreement, Some(true));
    assert_eq!(comparison.status, PurePhaseCrossValidationStatus::Complete);
    assert!(comparison.accepted);
}

#[test]
fn cross_validation_requires_named_complete_evidence_for_the_same_case() {
    let problem = scaled_boundary_problem(1.0);
    let validator = PurePhaseBoundaryValidator::default();
    let independent = validator
        .validate(&problem)
        .expect("controlled independent fixture must solve");
    let equilibrium = independent
        .equilibrium
        .as_ref()
        .expect("favorable candidate must have a finite reference state");
    let expected_tpd =
        independent.boundary.reaction_gibbs_at_absence / problem.candidate_stoichiometry();
    let evidence = CanonicalPurePhaseEvidence {
        gas_species: problem.gas_species().to_vec(),
        candidate_name: problem.candidate_name().to_string(),
        candidate_active: true,
        gas_moles: equilibrium.gas_moles.clone(),
        boundary_gas_moles: problem.gas_moles_at_absence().to_vec(),
        candidate_moles: equilibrium.candidate_moles,
        boundary_minimum_tpd: Some(expected_tpd),
    };

    let comparison = cross_validate_pure_phase(
        &problem,
        &evidence,
        validator,
        PurePhaseCrossValidationTolerances::default(),
    )
    .expect("high-level cross-validation must keep problem/result identity aligned");
    assert_eq!(comparison.status, PurePhaseCrossValidationStatus::Complete);
    assert!(comparison.accepted);

    let different_state = scaled_boundary_problem(2.0);
    let wrong_case = validator
        .validate(&different_state)
        .expect("second controlled fixture must solve");
    let error = compare_pure_phase_validation(
        &problem,
        &evidence,
        &wrong_case,
        PurePhaseCrossValidationTolerances::default(),
    )
    .expect_err("independent evidence from another state must be rejected");
    assert!(format!("{error:?}").contains("different pure-phase boundary case"));

    let different_thermochemistry = boundary_problem_for_target_k(8.0, 1.0);
    let wrong_thermochemistry = validator
        .validate(&different_thermochemistry)
        .expect("same-layout fixture with different Gibbs data must solve");
    let error = compare_pure_phase_validation(
        &problem,
        &evidence,
        &wrong_thermochemistry,
        PurePhaseCrossValidationTolerances::default(),
    )
    .expect_err("independent evidence with different Gibbs closures must be rejected");
    assert!(format!("{error:?}").contains("boundary does not match"));

    let another_temperature = boundary_problem_for_target_k_at_temperature(4.0, 1.0, 1_001.0);
    let wrong_temperature = validator
        .validate(&another_temperature)
        .expect("same physical fixture at another temperature must solve");
    let error = compare_pure_phase_validation(
        &problem,
        &evidence,
        &wrong_temperature,
        PurePhaseCrossValidationTolerances::default(),
    )
    .expect_err("independent evidence from another temperature must be rejected");
    assert!(format!("{error:?}").contains("different pure-phase boundary case"));
}

#[test]
fn species_permutation_preserves_boundary_root_and_final_equilibrium() {
    let validator = PurePhaseBoundaryValidator::default();
    let original_problem = scaled_boundary_problem(1.0)
        .with_element_composition(
            PurePhaseBoundaryElementComposition::new(
                vec!["X".to_string(), "Y".to_string()],
                DMatrix::from_row_slice(2, 2, &[1.0, 1.0, 1.0, 0.0]),
                vec![1.0, 2.0],
            )
            .unwrap(),
            PurePhaseBoundaryStructuralTolerances::default(),
        )
        .unwrap();
    let permuted_problem = permuted_boundary_problem_for_target_k(4.0, 1.0);
    let original_space = original_problem
        .validate_strict_independent_family(PurePhaseBoundaryStructuralTolerances::default())
        .unwrap();
    let permuted_space = permuted_problem
        .validate_strict_independent_family(PurePhaseBoundaryStructuralTolerances::default())
        .unwrap();
    assert_eq!(original_space.full_reaction_dimension, 1);
    assert_eq!(permuted_space.full_reaction_dimension, 1);
    assert_eq!(original_space.gas_only_reaction_dimension, 0);
    assert_eq!(permuted_space.gas_only_reaction_dimension, 0);
    let original = validator.validate(&original_problem).unwrap();
    let permuted = validator.validate(&permuted_problem).unwrap();
    let original_equilibrium = original.equilibrium.as_ref().unwrap();
    let permuted_equilibrium = permuted.equilibrium.as_ref().unwrap();

    let original_evidence = CanonicalPurePhaseEvidence {
        gas_species: original_problem.gas_species().to_vec(),
        candidate_name: original_problem.candidate_name().to_string(),
        candidate_active: true,
        gas_moles: original_equilibrium.gas_moles.clone(),
        boundary_gas_moles: original_problem.gas_moles_at_absence().to_vec(),
        candidate_moles: original_equilibrium.candidate_moles,
        boundary_minimum_tpd: Some(
            original.boundary.reaction_gibbs_at_absence
                / original_problem.candidate_stoichiometry(),
        ),
    };
    let permuted_evidence = CanonicalPurePhaseEvidence {
        gas_species: permuted_problem.gas_species().to_vec(),
        candidate_name: permuted_problem.candidate_name().to_string(),
        candidate_active: true,
        gas_moles: permuted_equilibrium.gas_moles.clone(),
        boundary_gas_moles: permuted_problem.gas_moles_at_absence().to_vec(),
        candidate_moles: permuted_equilibrium.candidate_moles,
        boundary_minimum_tpd: Some(
            permuted.boundary.reaction_gibbs_at_absence
                / permuted_problem.candidate_stoichiometry(),
        ),
    };

    let original_report = cross_validate_pure_phase(
        &original_problem,
        &original_evidence,
        validator,
        PurePhaseCrossValidationTolerances::default(),
    )
    .unwrap();
    let permuted_report = cross_validate_pure_phase(
        &permuted_problem,
        &permuted_evidence,
        validator,
        PurePhaseCrossValidationTolerances::default(),
    )
    .unwrap();

    assert_eq!(
        original_report.status,
        PurePhaseCrossValidationStatus::Complete
    );
    assert_eq!(
        permuted_report.status,
        PurePhaseCrossValidationStatus::Complete
    );
    assert!((original_equilibrium.gas_moles[0] - permuted_equilibrium.gas_moles[1]).abs() < 1e-12);
    assert!((original_equilibrium.gas_moles[1] - permuted_equilibrium.gas_moles[0]).abs() < 1e-12);
    assert!(
        (original_equilibrium.candidate_moles - permuted_equilibrium.candidate_moles).abs() < 1e-12
    );
}

#[test]
fn cross_validation_rejects_permuted_components_and_marks_missing_axes_incomplete() {
    let problem = analytic_boundary_problem(600.0, 600.0, 0.02);
    let validator = PurePhaseBoundaryValidator::default();
    let independent = validator
        .validate(&problem)
        .expect("boundary fixture must validate");
    let incomplete = CanonicalPurePhaseEvidence {
        gas_species: problem.gas_species().to_vec(),
        candidate_name: problem.candidate_name().to_string(),
        candidate_active: true,
        gas_moles: vec![1.0, 1.0],
        boundary_gas_moles: vec![1.0, 1.0],
        candidate_moles: 0.25,
        boundary_minimum_tpd: None,
    };
    let comparison = compare_pure_phase_validation(
        &problem,
        &incomplete,
        &independent,
        PurePhaseCrossValidationTolerances::default(),
    )
    .expect("missing optional evidence must report a status, not fail comparison");
    assert_eq!(
        comparison.status,
        PurePhaseCrossValidationStatus::InsufficientEvidence
    );
    assert!(!comparison.accepted);

    let mut permuted = incomplete;
    permuted.gas_species.reverse();
    let error = compare_pure_phase_validation(
        &problem,
        &permuted,
        &independent,
        PurePhaseCrossValidationTolerances::default(),
    )
    .expect_err("equal-length but permuted component evidence must be rejected");
    assert!(format!("{error:?}").contains("component identities"));
}

#[test]
fn cross_validation_report_localizes_thermodynamic_topology_and_composition_failures() {
    let problem = scaled_boundary_problem(1.0);
    let independent = PurePhaseBoundaryValidator::default()
        .validate(&problem)
        .expect("favorable independent fixture must solve");
    let equilibrium = independent
        .equilibrium
        .as_ref()
        .expect("favorable fixture must have a finite two-phase state");
    let expected_tpd =
        independent.boundary.reaction_gibbs_at_absence / problem.candidate_stoichiometry();
    let tolerances = PurePhaseCrossValidationTolerances::default();

    let wrong_tpd = compare_pure_phase_validation(
        &problem,
        &CanonicalPurePhaseEvidence {
            gas_species: problem.gas_species().to_vec(),
            candidate_name: problem.candidate_name().to_string(),
            candidate_active: true,
            gas_moles: equilibrium.gas_moles.clone(),
            boundary_gas_moles: problem.gas_moles_at_absence().to_vec(),
            candidate_moles: equilibrium.candidate_moles,
            boundary_minimum_tpd: Some(expected_tpd + 1.0),
        },
        &independent,
        tolerances,
    )
    .expect("diagnostic fixture must compare");
    assert_eq!(wrong_tpd.thermodynamic_agreement, Some(false));
    assert_eq!(wrong_tpd.topology_agreement, Some(true));
    assert_eq!(wrong_tpd.composition_agreement, Some(true));
    assert_eq!(wrong_tpd.status, PurePhaseCrossValidationStatus::Disagreed);
    assert!(!wrong_tpd.accepted);

    let wrong_topology = compare_pure_phase_validation(
        &problem,
        &CanonicalPurePhaseEvidence {
            gas_species: problem.gas_species().to_vec(),
            candidate_name: problem.candidate_name().to_string(),
            candidate_active: false,
            gas_moles: equilibrium.gas_moles.clone(),
            boundary_gas_moles: problem.gas_moles_at_absence().to_vec(),
            candidate_moles: equilibrium.candidate_moles,
            boundary_minimum_tpd: Some(expected_tpd),
        },
        &independent,
        tolerances,
    )
    .expect("diagnostic fixture must compare");
    assert_eq!(wrong_topology.thermodynamic_agreement, Some(true));
    assert_eq!(wrong_topology.topology_agreement, Some(false));
    assert_eq!(wrong_topology.composition_agreement, Some(true));
    assert_eq!(
        wrong_topology.status,
        PurePhaseCrossValidationStatus::Disagreed
    );
    assert!(!wrong_topology.accepted);

    let mut wrong_gas_moles = equilibrium.gas_moles.clone();
    wrong_gas_moles[0] += 1.0;
    let wrong_composition = compare_pure_phase_validation(
        &problem,
        &CanonicalPurePhaseEvidence {
            gas_species: problem.gas_species().to_vec(),
            candidate_name: problem.candidate_name().to_string(),
            candidate_active: true,
            gas_moles: wrong_gas_moles,
            boundary_gas_moles: problem.gas_moles_at_absence().to_vec(),
            candidate_moles: equilibrium.candidate_moles,
            boundary_minimum_tpd: Some(expected_tpd),
        },
        &independent,
        tolerances,
    )
    .expect("diagnostic fixture must compare");
    assert_eq!(wrong_composition.thermodynamic_agreement, Some(true));
    assert_eq!(wrong_composition.topology_agreement, Some(true));
    assert_eq!(wrong_composition.composition_agreement, Some(false));
    assert_eq!(
        wrong_composition.status,
        PurePhaseCrossValidationStatus::Disagreed
    );
    assert!(!wrong_composition.accepted);
}

#[test]
fn reaction_coordinate_scaling_preserves_physical_and_per_candidate_quantities() {
    let validator = PurePhaseBoundaryValidator::default();
    let reference_problem = scaled_boundary_problem(1.0);
    let reference_boundary =
        evaluate_pure_phase_boundary(&reference_problem, validator.tolerances).unwrap();
    let reference_equilibrium = validator
        .validate(&reference_problem)
        .unwrap()
        .equilibrium
        .expect("reference candidate formation is favorable");
    let reference_per_candidate =
        reference_boundary.reaction_gibbs_at_absence / reference_problem.candidate_stoichiometry();
    let reference_gas_total = reference_equilibrium.gas_moles.iter().sum::<f64>();
    let reference_mole_fractions = reference_equilibrium
        .gas_moles
        .iter()
        .map(|moles| moles / reference_gas_total)
        .collect::<Vec<_>>();

    for scale in [0.5, 1.0, 2.0] {
        let problem = scaled_boundary_problem(scale);
        let boundary = evaluate_pure_phase_boundary(&problem, validator.tolerances).unwrap();
        let equilibrium = validator
            .validate(&problem)
            .unwrap()
            .equilibrium
            .expect("rescaled candidate formation is favorable");

        assert!(
            (boundary.ln_q_at_absence - scale * reference_boundary.ln_q_at_absence).abs() < 1e-12
        );
        assert!((boundary.ln_k - scale * reference_boundary.ln_k).abs() < 1e-12);
        assert!(
            (boundary.log_residual_at_absence - scale * reference_boundary.log_residual_at_absence)
                .abs()
                < 1e-12
        );
        assert!(
            (boundary.reaction_gibbs_at_absence
                - scale * reference_boundary.reaction_gibbs_at_absence)
                .abs()
                < 1e-8
        );
        assert!(
            (boundary.reaction_gibbs_at_absence / problem.candidate_stoichiometry()
                - reference_per_candidate)
                .abs()
                < 1e-8
        );
        assert!(
            (equilibrium.extent - reference_equilibrium.extent / scale).abs() < 1e-10,
            "extent is a coordinate and must scale inversely with nu"
        );
        for (actual, reference) in equilibrium
            .gas_moles
            .iter()
            .zip(&reference_equilibrium.gas_moles)
        {
            assert!((actual - reference).abs() < 1e-10);
        }
        assert!(
            (equilibrium.candidate_moles - reference_equilibrium.candidate_moles).abs() < 1e-10
        );
        let gas_total = equilibrium.gas_moles.iter().sum::<f64>();
        for (fraction, reference_fraction) in equilibrium
            .gas_moles
            .iter()
            .map(|moles| moles / gas_total)
            .zip(&reference_mole_fractions)
        {
            assert!(
                (fraction - reference_fraction).abs() < 1e-12,
                "gas mole fractions are physical quantities and must not depend on reaction-coordinate scale"
            );
        }
        assert_eq!(boundary.prediction, reference_boundary.prediction);
    }
}

#[test]
fn pressure_ratio_controls_lnq_and_joint_pressure_reference_scaling_is_invariant() {
    let temperature = 1_000.0;
    let reference_pressure = 101_325.0;
    let baseline = boundary_problem_for_target_k_at_conditions(
        4.0,
        1.0,
        temperature,
        reference_pressure,
        reference_pressure,
    );
    let doubled_pressure = boundary_problem_for_target_k_at_conditions(
        4.0,
        1.0,
        temperature,
        2.0 * reference_pressure,
        reference_pressure,
    );
    let jointly_scaled = boundary_problem_for_target_k_at_conditions(
        4.0,
        1.0,
        temperature,
        2.0 * reference_pressure,
        2.0 * reference_pressure,
    );
    let tolerances = PurePhaseBoundaryTolerances::default();
    let baseline_report = evaluate_pure_phase_boundary(&baseline, tolerances).unwrap();
    let doubled_report = evaluate_pure_phase_boundary(&doubled_pressure, tolerances).unwrap();
    let jointly_scaled_report = evaluate_pure_phase_boundary(&jointly_scaled, tolerances).unwrap();

    let delta_nu_gas = baseline.gas_stoichiometry().iter().sum::<f64>();
    let expected_delta_lnq = delta_nu_gas * 2.0_f64.ln();
    assert!(
        (doubled_report.ln_q_at_absence - baseline_report.ln_q_at_absence - expected_delta_lnq)
            .abs()
            < 1e-12
    );
    assert!(
        (doubled_report.reaction_gibbs_at_absence
            - baseline_report.reaction_gibbs_at_absence
            - MOLAR_GAS_CONSTANT * temperature * expected_delta_lnq)
            .abs()
            < 1e-8
    );

    let canonical_baseline_tpd = canonical_tpd_for_target_k_at_conditions(
        4.0,
        temperature,
        reference_pressure,
        reference_pressure,
    );
    let canonical_doubled_tpd = canonical_tpd_for_target_k_at_conditions(
        4.0,
        temperature,
        2.0 * reference_pressure,
        reference_pressure,
    );
    assert!((canonical_baseline_tpd - baseline_report.reaction_gibbs_at_absence).abs() < 1e-8);
    assert!((canonical_doubled_tpd - doubled_report.reaction_gibbs_at_absence).abs() < 1e-8);

    assert!(
        (jointly_scaled_report.ln_q_at_absence - baseline_report.ln_q_at_absence).abs() < 1e-12
    );
    assert!((jointly_scaled_report.ln_k - baseline_report.ln_k).abs() < 1e-12);
    assert_eq!(jointly_scaled_report.prediction, baseline_report.prediction);
}

#[test]
fn inert_dilution_changes_q_and_matches_canonical_tpd() {
    let validator = PurePhaseBoundaryValidator::default();
    let low_inert = 1.0;
    let high_inert = 3.0;
    let low_problem = inert_dilution_boundary_problem(low_inert);
    let high_problem = inert_dilution_boundary_problem(high_inert);
    let low = validator.validate(&low_problem).unwrap();
    let high = validator.validate(&high_problem).unwrap();

    let low_boundary = &low.boundary;
    let high_boundary = &high.boundary;
    let delta_nu_gas = low_problem.gas_stoichiometry().iter().sum::<f64>();
    let low_total = 2.0 + low_inert;
    let high_total = 2.0 + high_inert;
    let expected_delta_ln_q = -delta_nu_gas * (high_total / low_total).ln();
    assert!(
        (high_boundary.ln_q_at_absence - low_boundary.ln_q_at_absence - expected_delta_ln_q).abs()
            < 1e-12,
        "an inert gas changes reacting-species activities through dilution"
    );

    let low_tpd = canonical_inert_dilution_tpd(low_inert);
    let high_tpd = canonical_inert_dilution_tpd(high_inert);
    assert!(
        (low_tpd - low_boundary.reaction_gibbs_at_absence).abs() < 1e-8,
        "independent delta-G must use the same diluted activities as canonical TPD"
    );
    assert!(
        (high_tpd - high_boundary.reaction_gibbs_at_absence).abs() < 1e-8,
        "independent delta-G must use the same diluted activities as canonical TPD"
    );
    assert!((high_tpd - low_tpd - MOLAR_GAS_CONSTANT * 1_000.0 * expected_delta_ln_q).abs() < 1e-8);
}

#[test]
fn global_inventory_scaling_preserves_boundary_and_scales_equilibrium() {
    let validator = PurePhaseBoundaryValidator::default();
    let reference_problem = scaled_boundary_problem(1.0);
    let reference = validator.validate(&reference_problem).unwrap();
    let reference_equilibrium = reference
        .equilibrium
        .as_ref()
        .expect("the synthetic candidate is favorable at the reference inventory");
    let reference_total = reference_equilibrium.gas_moles.iter().sum::<f64>();
    let reference_fractions = reference_equilibrium
        .gas_moles
        .iter()
        .map(|moles| moles / reference_total)
        .collect::<Vec<_>>();

    for inventory_scale in [1e-9, 1e-3, 1.0, 1e3, 1e9] {
        let problem = PurePhaseBoundaryProblem::new(
            vec!["A".to_string(), "B".to_string()],
            vec![inventory_scale, inventory_scale],
            vec![-2.0, 1.0],
            1.0,
            vec![constant_gibbs(0.0), constant_gibbs(0.0)],
            constant_gibbs(-MOLAR_GAS_CONSTANT * 1_000.0 * 4.0_f64.ln()),
            EquilibriumConditions::new(1_000.0, 101_325.0, 101_325.0).unwrap(),
            "S",
        )
        .unwrap();
        let actual = validator.validate(&problem).unwrap();
        let equilibrium = actual
            .equilibrium
            .as_ref()
            .expect("inventory scaling must not change the favorable topology");

        assert!(
            (actual.boundary.ln_q_at_absence - reference.boundary.ln_q_at_absence).abs() < 1e-12,
            "scale={inventory_scale:e}"
        );
        assert!(
            (actual.boundary.reaction_gibbs_at_absence
                - reference.boundary.reaction_gibbs_at_absence)
                .abs()
                < 1e-8,
            "scale={inventory_scale:e}"
        );
        assert_eq!(actual.boundary.prediction, reference.boundary.prediction);
        for ((actual_moles, reference_moles), fraction) in equilibrium
            .gas_moles
            .iter()
            .zip(&reference_equilibrium.gas_moles)
            .zip(&reference_fractions)
        {
            let expected_moles = inventory_scale * reference_moles;
            let scale_aware_tolerance = 1e-12_f64.max(expected_moles.abs() * 1e-10);
            assert!(
                (actual_moles - expected_moles).abs() <= scale_aware_tolerance,
                "scale={inventory_scale:e}, actual={actual_moles:e}, expected={expected_moles:e}"
            );
            let total = equilibrium.gas_moles.iter().sum::<f64>();
            assert!(
                (actual_moles / total - fraction).abs() < 1e-12,
                "gas composition must be homogeneous in the total inventory"
            );
        }
        let expected_candidate = inventory_scale * reference_equilibrium.candidate_moles;
        let candidate_tolerance = 1e-12_f64.max(expected_candidate.abs() * 1e-10);
        assert!(
            (equilibrium.candidate_moles - expected_candidate).abs() <= candidate_tolerance,
            "candidate inventory must scale with the whole system"
        );
    }
}

#[test]
fn deterministic_generated_boundary_matrix_preserves_independent_contracts() {
    // This is intentionally a fixed generated corpus rather than random fuzz:
    // its descriptor is sufficient to reproduce any future regression.
    let cases: [(&str, f64, f64, f64, f64); 5] = [
        ("appear-unit", 4.0, 1.0, 1.0, 1.0),
        ("appear-scaled", 4.0, 2.0, 1.0, 1e-6),
        ("boundary", 2.0, 0.5, 1.0, 1.0),
        ("stable-low-k", 1.0, 1.0, 1.0, 1e6),
        ("pressure-shifted", 4.0, 1.0, 2.0, 1.0),
    ];
    let validator = PurePhaseBoundaryValidator::default();

    for (fixture_id, target_k, stoichiometric_scale, pressure_factor, inventory_scale) in cases {
        let temperature = 1_000.0;
        let pressure = 101_325.0 * pressure_factor;
        let problem = PurePhaseBoundaryProblem::new(
            vec!["A".to_string(), "B".to_string()],
            vec![inventory_scale, inventory_scale],
            vec![-2.0 * stoichiometric_scale, stoichiometric_scale],
            stoichiometric_scale,
            vec![constant_gibbs(0.0), constant_gibbs(0.0)],
            constant_gibbs(-MOLAR_GAS_CONSTANT * temperature * target_k.ln()),
            EquilibriumConditions::new(temperature, pressure, 101_325.0).unwrap(),
            "S",
        )
        .unwrap()
        // A/B/S is deliberately supplied with an independent element matrix:
        // A = X + Y, B = X, and S = X + 2Y. Together with 2A -> B + S,
        // this keeps the generated corpus honest about rank and conservation
        // instead of treating a matching scalar residual as sufficient evidence.
        .with_element_composition(
            PurePhaseBoundaryElementComposition::new(
                vec!["X".to_string(), "Y".to_string()],
                DMatrix::from_row_slice(2, 2, &[1.0, 1.0, 1.0, 0.0]),
                vec![1.0, 2.0],
            )
            .unwrap(),
            PurePhaseBoundaryStructuralTolerances::default(),
        )
        .unwrap();
        let reaction_space = problem
            .validate_strict_independent_family(PurePhaseBoundaryStructuralTolerances::default())
            .unwrap_or_else(|error| panic!("fixture={fixture_id}: {error:?}"));
        assert_eq!(
            reaction_space.full_reaction_dimension, 1,
            "fixture={fixture_id}"
        );
        assert_eq!(
            reaction_space.gas_only_reaction_dimension, 0,
            "fixture={fixture_id}"
        );
        assert!(
            reaction_space
                .element_balance_residuals
                .iter()
                .all(|residual| residual.abs() <= 1e-12),
            "fixture={fixture_id}, residuals={:?}",
            reaction_space.element_balance_residuals
        );
        let result = validator.validate(&problem).unwrap_or_else(|error| {
            panic!(
                "fixture={fixture_id}, K={target_k:e}, nu_scale={stoichiometric_scale:e}, pressure_factor={pressure_factor:e}, inventory_scale={inventory_scale:e}: {error:?}"
            )
        });
        let expected_tpd =
            result.boundary.reaction_gibbs_at_absence / problem.candidate_stoichiometry();
        let canonical_tpd =
            canonical_tpd_for_target_k_at_conditions(target_k, temperature, pressure, 101_325.0);
        assert!(
            (canonical_tpd - expected_tpd).abs() < 1e-8,
            "fixture={fixture_id}, K={target_k:e}, nu_scale={stoichiometric_scale:e}, pressure_factor={pressure_factor:e}, inventory_scale={inventory_scale:e}"
        );
        if let Some(equilibrium) = result.equilibrium.as_ref() {
            assert!(equilibrium.gas_moles.iter().all(|moles| *moles > 0.0));
            assert!(
                equilibrium.log_residual.abs()
                    <= validator.tolerances.max_abs_equilibrium_log_residual,
                "fixture={fixture_id}"
            );
        }

        // The boundary TPD is analytically prescribed separately from the
        // reaction residual above. For active cases this synthetic corpus
        // deliberately reuses the scalar finite extent only to populate a
        // complete final-composition record. Production-vs-independent
        // composition evidence is covered by the dedicated runner stories;
        // this test owns the exhaustive comparator contract over generated
        // rank/conservation cases.
        let canonical = if let Some(equilibrium) = result.equilibrium.as_ref() {
            CanonicalPurePhaseEvidence {
                gas_species: problem.gas_species().to_vec(),
                candidate_name: problem.candidate_name().to_string(),
                candidate_active: true,
                gas_moles: equilibrium.gas_moles.clone(),
                boundary_gas_moles: problem.gas_moles_at_absence().to_vec(),
                candidate_moles: equilibrium.candidate_moles,
                boundary_minimum_tpd: Some(canonical_tpd),
            }
        } else {
            CanonicalPurePhaseEvidence {
                gas_species: problem.gas_species().to_vec(),
                candidate_name: problem.candidate_name().to_string(),
                candidate_active: false,
                gas_moles: problem.gas_moles_at_absence().to_vec(),
                boundary_gas_moles: problem.gas_moles_at_absence().to_vec(),
                candidate_moles: 0.0,
                boundary_minimum_tpd: Some(canonical_tpd),
            }
        };
        let comparison = compare_pure_phase_validation(
            &problem,
            &canonical,
            &result,
            PurePhaseCrossValidationTolerances::default(),
        )
        .unwrap_or_else(|error| panic!("fixture={fixture_id}: {error:?}"));
        assert_eq!(
            comparison.status,
            PurePhaseCrossValidationStatus::Complete,
            "fixture={fixture_id}: {comparison:?}"
        );
        assert!(comparison.accepted, "fixture={fixture_id}: {comparison:?}");
    }
}

#[test]
fn p10_1_pt_cross_validation_uses_scale_aware_mole_tolerances() {
    let problem = scaled_boundary_problem(1.0);
    let validator = PurePhaseBoundaryValidator::default();
    let independent = validator
        .validate(&problem)
        .expect("large controlled inventory must retain a finite equilibrium");
    let equilibrium = independent
        .equilibrium
        .as_ref()
        .expect("favorable controlled candidate must form");
    let tpd = independent.boundary.reaction_gibbs_at_absence / problem.candidate_stoichiometry();
    let scaled_difference = 5e-5;
    let scale_aware_tolerances = PurePhaseCrossValidationTolerances {
        // Deliberately make the absolute floor much smaller than the controlled
        // perturbation so this test proves that the relative branch is used.
        max_abs_gas_mole_delta: 1e-9,
        max_relative_gas_mole_delta: 1e-4,
        max_abs_candidate_mole_delta: 1e-9,
        max_relative_candidate_mole_delta: 1e-4,
        ..PurePhaseCrossValidationTolerances::default()
    };
    let evidence = CanonicalPurePhaseEvidence {
        gas_species: problem.gas_species().to_vec(),
        candidate_name: problem.candidate_name().to_string(),
        candidate_active: true,
        gas_moles: equilibrium
            .gas_moles
            .iter()
            .map(|moles| moles * (1.0 + scaled_difference))
            .collect(),
        boundary_gas_moles: problem.gas_moles_at_absence().to_vec(),
        candidate_moles: equilibrium.candidate_moles * (1.0 + scaled_difference),
        boundary_minimum_tpd: Some(tpd),
    };

    let scale_aware =
        compare_pure_phase_validation(&problem, &evidence, &independent, scale_aware_tolerances)
            .expect("finite evidence must compare");
    assert_eq!(scale_aware.status, PurePhaseCrossValidationStatus::Complete);
    assert!(scale_aware.accepted);
    assert!(
        scale_aware
            .max_abs_gas_mole_delta
            .expect("finite equilibrium has a gas comparison")
            > scale_aware_tolerances.max_abs_gas_mole_delta,
        "the relative tolerance, not the absolute floor, must accept this scaled case"
    );

    let absolute_only = compare_pure_phase_validation(
        &problem,
        &evidence,
        &independent,
        PurePhaseCrossValidationTolerances {
            max_relative_gas_mole_delta: 0.0,
            max_relative_candidate_mole_delta: 0.0,
            ..scale_aware_tolerances
        },
    )
    .expect("strict evidence must compare rather than fail structurally");
    assert_eq!(absolute_only.composition_agreement, Some(false));
    assert_eq!(
        absolute_only.status,
        PurePhaseCrossValidationStatus::Disagreed
    );
}

#[test]
fn p10_1_pt_cross_validation_rejects_invalid_tolerance_matrix() {
    let problem = scaled_boundary_problem(1.0);
    let independent = PurePhaseBoundaryValidator::default()
        .validate(&problem)
        .expect("controlled independent fixture must solve");
    let equilibrium = independent
        .equilibrium
        .as_ref()
        .expect("favorable fixture must have a finite reference state");
    let evidence = CanonicalPurePhaseEvidence {
        gas_species: problem.gas_species().to_vec(),
        candidate_name: problem.candidate_name().to_string(),
        candidate_active: true,
        gas_moles: equilibrium.gas_moles.clone(),
        boundary_gas_moles: problem.gas_moles_at_absence().to_vec(),
        candidate_moles: equilibrium.candidate_moles,
        boundary_minimum_tpd: Some(
            independent.boundary.reaction_gibbs_at_absence / problem.candidate_stoichiometry(),
        ),
    };
    let invalid_tolerances = [
        PurePhaseCrossValidationTolerances {
            max_abs_gas_mole_delta: 0.0,
            ..PurePhaseCrossValidationTolerances::default()
        },
        PurePhaseCrossValidationTolerances {
            max_abs_candidate_mole_delta: -1.0,
            ..PurePhaseCrossValidationTolerances::default()
        },
        PurePhaseCrossValidationTolerances {
            max_relative_gas_mole_delta: f64::NAN,
            ..PurePhaseCrossValidationTolerances::default()
        },
        PurePhaseCrossValidationTolerances {
            max_relative_candidate_mole_delta: f64::INFINITY,
            ..PurePhaseCrossValidationTolerances::default()
        },
    ];
    for tolerances in invalid_tolerances {
        let error = compare_pure_phase_validation(&problem, &evidence, &independent, tolerances)
            .expect_err("invalid P,T cross-validation tolerances must be typed errors");
        assert!(matches!(error, ReactionExtentError::InvalidProblem { .. }));
    }
}

#[test]
fn boundary_classification_respects_independent_tolerance_edges() {
    let tolerance = 1e-6;
    let validator = PurePhaseBoundaryValidator {
        tolerances: PurePhaseBoundaryTolerances {
            max_abs_boundary_log_residual: tolerance,
            ..PurePhaseBoundaryTolerances::default()
        },
        ..PurePhaseBoundaryValidator::default()
    };
    let ln_q = 2.0_f64.ln();

    for (residual, expected) in [
        (-1.01 * tolerance, PurePhaseBoundaryPrediction::ShouldAppear),
        (-0.99 * tolerance, PurePhaseBoundaryPrediction::Boundary),
        (0.0, PurePhaseBoundaryPrediction::Boundary),
        (0.99 * tolerance, PurePhaseBoundaryPrediction::Boundary),
        (
            1.01 * tolerance,
            PurePhaseBoundaryPrediction::StableInactive,
        ),
    ] {
        // ln(Q) - ln(K) is prescribed directly, independently of phase-control
        // hysteresis thresholds used by the production outer loop.
        let target_k = (ln_q - residual).exp();
        let problem = boundary_problem_for_target_k(target_k, 1.0);
        let actual = validator.validate(&problem).unwrap();
        assert_eq!(
            actual.boundary.prediction, expected,
            "residual={residual:.16e}, tolerance={tolerance:.16e}"
        );
    }
}

#[test]
fn fixed_topology_gibbs_solution_matches_independent_finite_equilibrium_under_scaling() {
    let canonical_moles = canonical_fixed_phase_solution_for_target_k(4.0);
    let validator = PurePhaseBoundaryValidator::default();

    for scale in [0.5, 1.0, 2.0] {
        let problem = scaled_boundary_problem(scale);
        let independent = validator
            .validate(&problem)
            .expect("scaled independent reaction must be solvable");
        let equilibrium = independent
            .equilibrium
            .as_ref()
            .expect("favorable candidate must form at finite extent");
        let comparison = compare_pure_phase_validation(
            &problem,
            &CanonicalPurePhaseEvidence {
                gas_species: problem.gas_species().to_vec(),
                candidate_name: problem.candidate_name().to_string(),
                candidate_active: canonical_moles[2] > 1e-12,
                gas_moles: canonical_moles[..2].to_vec(),
                boundary_gas_moles: problem.gas_moles_at_absence().to_vec(),
                candidate_moles: canonical_moles[2],
                boundary_minimum_tpd: None,
            },
            &independent,
            PurePhaseCrossValidationTolerances {
                max_abs_gas_mole_delta: 1e-8,
                max_abs_candidate_mole_delta: 1e-8,
                ..PurePhaseCrossValidationTolerances::default()
            },
        )
        .expect("strict family has directly comparable final states");

        assert_eq!(comparison.topology_agreement, Some(true), "scale={scale}");
        assert_eq!(
            comparison.composition_agreement,
            Some(true),
            "scale={scale}"
        );
        assert_eq!(comparison.thermodynamic_agreement, None);
        assert_eq!(
            comparison.status,
            PurePhaseCrossValidationStatus::ConsistentButPartial,
            "scale={scale}: {comparison:?}"
        );
        assert!(!comparison.accepted, "scale={scale}: {comparison:?}");
        assert!(equilibrium.candidate_moles > 0.0 && canonical_moles[2] > 0.0);
    }
}

#[test]
fn production_phase_manager_hysteresis_matrix_respects_history_aware_topology() {
    let manager = PhaseManager::new(1e-6, -1.0, 1.0);

    let inactive = phase_set(&[true, false]);
    assert!(matches!(
        manager
            .classify_phases(
                &[1.0, 0.0],
                &lifecycle_stability_reports(&[true, false], -2.0),
                &inactive,
            )
            .unwrap(),
        PhaseTransitionPlan::Activate { phase } if phase.index() == 1
    ));
    assert!(matches!(
        manager
            .classify_phases(
                &[1.0, 0.0],
                &lifecycle_stability_reports(&[true, false], 0.0),
                &inactive,
            )
            .unwrap(),
        PhaseTransitionPlan::NoTransition {
            reason: PhaseStabilityReason::NoPhaseTransitionNeeded
        }
    ));

    let active = phase_set(&[true, true]);
    assert!(matches!(
        manager
            .classify_phases(
                &[1.0, 1e-8],
                &lifecycle_stability_reports(&[true, true], 0.0),
                &active,
            )
            .unwrap(),
        PhaseTransitionPlan::Hold { phase } if phase.index() == 1
    ));
    assert!(matches!(
        manager
            .classify_phases(
                &[1.0, 1e-8],
                &lifecycle_stability_reports(&[true, true], 2.0),
                &active,
            )
            .unwrap(),
        PhaseTransitionPlan::Deactivate { phase } if phase.index() == 1
    ));
}

#[test]
fn production_phase_manager_hysteresis_edges_are_strict() {
    let manager = PhaseManager::new(1e-6, -1.0, 1.0);
    let inactive = phase_set(&[true, false]);
    assert!(matches!(
        manager
            .classify_phases(
                &[1.0, 0.0],
                &lifecycle_stability_reports(&[true, false], -1.0),
                &inactive,
            )
            .unwrap(),
        PhaseTransitionPlan::NoTransition {
            reason: PhaseStabilityReason::NoPhaseTransitionNeeded
        }
    ));

    let active = phase_set(&[true, true]);
    assert!(matches!(
        manager
            .classify_phases(
                &[1.0, 1e-8],
                &lifecycle_stability_reports(&[true, true], 1.0),
                &active,
            )
            .unwrap(),
        PhaseTransitionPlan::Hold { phase } if phase.index() == 1
    ));
}

#[test]
fn same_in_band_point_uses_previous_active_set_history() {
    let temperature = 1_001.0;
    let in_band_tpd = -0.5;
    let in_band_k = 2.0 * (-in_band_tpd / (MOLAR_GAS_CONSTANT * temperature)).exp();
    let policy = PhaseManager::new(1e-6, -1.0, 1.0);

    // First establish an accepted active history at a clearly unstable point.
    let mut continued_runner = production_phase_control_runner(4.0);
    *continued_runner.configure_phase_control() = policy.clone();
    let accepted_active = continued_runner
        .solve()
        .expect("the strongly favorable point must activate the pure phase");
    let accepted_seed = LogMolesInitialGuess::new(accepted_active.solution.log_moles().to_vec())
        .expect("an accepted solution is always a valid continuation seed");
    let accepted_phase_set = accepted_active.phase_control_report.final_phase_set.clone();
    assert_eq!(accepted_phase_set.active_mask(), &[true, true]);

    retarget_production_runner(&mut continued_runner, temperature, in_band_k, accepted_seed);
    continued_runner
        .set_continuation_phase_set(accepted_phase_set)
        .expect("accepted phase set matches the unchanged prepared layout");
    let continued = continued_runner
        .solve()
        .expect("an in-band continuation point must remain solvable");
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
        "the keep threshold retains an accepted active phase inside the hysteresis band"
    );
    assert!(continued.phase_control_report.transitions.is_empty());

    // The same thermodynamic point starts inactive when it has no accepted
    // history: its TPD lies above dg_create and therefore cannot activate.
    let mut fresh_runner = production_phase_control_runner(4.0);
    *fresh_runner.configure_phase_control() = policy;
    let fresh_seed = LogMolesInitialGuess::from_moles(
        &[1.0, 1.0, PHASE_CONTROL_TRACE_MOLE_FLOOR],
        PHASE_CONTROL_TRACE_MOLE_FLOOR,
    )
    .expect("fresh trace seed must be valid");
    retarget_production_runner(&mut fresh_runner, temperature, in_band_k, fresh_seed);
    let fresh = fresh_runner
        .solve()
        .expect("fresh gas-only in-band point must be stable");
    assert_eq!(
        fresh.phase_control_report.initial_phase_set.active_mask(),
        &[true, false]
    );
    assert_eq!(
        fresh.phase_control_report.final_phase_set.active_mask(),
        &[true, false]
    );
    assert!(fresh.phase_control_report.transitions.is_empty());
}

#[test]
fn ascending_and_descending_boundary_sweeps_expose_hysteresis() {
    let policy = PhaseManager::new(1e-6, -1.0, 1.0);
    let in_band_k = target_k_for_absence_tpd(1_000.0, -0.5);
    // The endpoints deliberately stay far from the nonlinear pure-phase
    // boundary. The middle point is inside the hysteresis band; it is the
    // only point whose topology is allowed to depend on the path history.
    let ascending_points = [(1_000.0, 4.0), (1_001.0, in_band_k), (1_002.0, 1e-4)];
    let descending_points = [(1_002.0, 1e-4), (1_001.0, in_band_k), (1_000.0, 4.0)];

    let mut ascending_runner = production_phase_control_runner(4.0);
    *ascending_runner.configure_phase_control() = policy.clone();
    let mut ascending_outcomes: Vec<PreparedPhaseControlOutcome> = Vec::new();
    for (index, (temperature, target_k)) in ascending_points.into_iter().enumerate() {
        let outcome = if let Some(previous) = ascending_outcomes.last() {
            let seed = LogMolesInitialGuess::new(previous.solution.log_moles().to_vec()).unwrap();
            retarget_production_runner(&mut ascending_runner, temperature, target_k, seed);
            ascending_runner
                .set_continuation_phase_set(previous.phase_control_report.final_phase_set.clone())
                .unwrap();
            ascending_runner
                .solve()
                .unwrap_or_else(|error| panic!("ascending point {index} failed: {error:?}"))
        } else {
            ascending_runner.solve().unwrap()
        };
        assert!(
            outcome.acceptance_report.complementarity.satisfied,
            "ascending point {index} must publish complementarity evidence"
        );
        ascending_outcomes.push(outcome);
    }

    assert_eq!(
        ascending_outcomes[0]
            .phase_control_report
            .final_phase_set
            .active_mask(),
        &[true, true]
    );
    assert_eq!(
        ascending_outcomes[1]
            .phase_control_report
            .initial_phase_set
            .active_mask(),
        &[true, true]
    );
    assert_eq!(
        ascending_outcomes[1]
            .phase_control_report
            .final_phase_set
            .active_mask(),
        &[true, true],
        "the ascending in-band point retains the previously active candidate"
    );
    assert!(
        ascending_outcomes[1]
            .phase_control_report
            .transitions
            .is_empty()
    );
    assert_eq!(
        ascending_outcomes[2]
            .phase_control_report
            .final_phase_set
            .active_mask(),
        &[true, false],
        "the positive TPD endpoint must deactivate the vanishing pure phase"
    );
    assert!(matches!(
        ascending_outcomes[2].phase_control_report.transitions.as_slice(),
        [transition] if matches!(transition.reason, PhaseTransitionReason::BoundaryUnstableActivePhase { minimum_tpd, .. } if minimum_tpd > 0.0)
    ));

    let mut descending_runner = production_phase_control_runner(1e-4);
    *descending_runner.configure_phase_control() = policy;
    let mut descending_outcomes: Vec<PreparedPhaseControlOutcome> = Vec::new();
    for (index, (temperature, target_k)) in descending_points.into_iter().enumerate() {
        let outcome = if let Some(previous) = descending_outcomes.last() {
            let seed = LogMolesInitialGuess::new(previous.solution.log_moles().to_vec()).unwrap();
            retarget_production_runner(&mut descending_runner, temperature, target_k, seed);
            descending_runner
                .set_continuation_phase_set(previous.phase_control_report.final_phase_set.clone())
                .unwrap();
            descending_runner
                .solve()
                .unwrap_or_else(|error| panic!("descending point {index} failed: {error:?}"))
        } else {
            descending_runner.solve().unwrap()
        };
        assert!(
            outcome.acceptance_report.complementarity.satisfied,
            "descending point {index} must publish complementarity evidence"
        );
        descending_outcomes.push(outcome);
    }

    assert_eq!(
        descending_outcomes[0]
            .phase_control_report
            .final_phase_set
            .active_mask(),
        &[true, false]
    );
    assert_eq!(
        descending_outcomes[1]
            .phase_control_report
            .final_phase_set
            .active_mask(),
        &[true, false],
        "the descending in-band point has no active-phase history to retain"
    );
    assert!(
        descending_outcomes[1]
            .phase_control_report
            .transitions
            .is_empty()
    );
    assert_eq!(
        descending_outcomes[2]
            .phase_control_report
            .final_phase_set
            .active_mask(),
        &[true, true]
    );
    assert!(matches!(
        descending_outcomes[2].phase_control_report.transitions.as_slice(),
        [transition] if matches!(transition.reason, PhaseTransitionReason::UnstableInactivePhase { minimum_tpd } if minimum_tpd < 0.0)
    ));
}

#[test]
fn prepared_phase_control_activates_pure_phase_from_negative_tpd_and_matches_keq() {
    let independent_problem = scaled_boundary_problem(1.0);
    let independent = PurePhaseBoundaryValidator::default()
        .validate(&independent_problem)
        .expect("independent K_eq reference must be solvable");
    let expected = independent
        .equilibrium
        .as_ref()
        .expect("K=4 makes the pure candidate favorable at the gas-only state");

    let outcome = production_phase_control_runner(4.0)
        .solve()
        .expect("production phase control must activate the TPD-unstable pure phase");
    let moles = outcome.solution.moles();
    assert_eq!(
        outcome.phase_control_report.transitions.len(),
        1,
        "simple pure-phase fixture must converge through exactly one activation transition"
    );
    let transition = outcome
        .phase_control_report
        .transitions
        .first()
        .expect("negative inactive-phase TPD must create one lifecycle transition");

    assert_eq!(
        outcome.phase_control_report.initial_phase_set.active_mask(),
        &[true, false]
    );
    assert_eq!(
        outcome.phase_control_report.final_phase_set.active_mask(),
        &[true, true]
    );
    assert_eq!(transition.activated.len(), 1);
    assert!(transition.deactivated.is_empty());
    assert!(matches!(
        transition.reason,
        PhaseTransitionReason::UnstableInactivePhase { minimum_tpd } if minimum_tpd < 0.0
    ));
    let transition_tpd = match transition.reason {
        PhaseTransitionReason::UnstableInactivePhase { minimum_tpd } => minimum_tpd,
        _ => unreachable!("the assertion above fixes the transition variant"),
    };
    let independent_tpd = independent.boundary.reaction_gibbs_at_absence
        / independent_problem.candidate_stoichiometry();
    assert!(
        (transition_tpd - independent_tpd).abs() < 1e-8,
        "activation TPD must equal the independent per-candidate K_eq driving force"
    );

    // The transition record must retain the actual TPD minimizer and the
    // restart seed derived from it.  For a pure phase the minimizer is [1],
    // so the candidate's seeded amount is exactly the recorded phase total.
    assert_eq!(
        transition.incipient_composition.as_deref(),
        Some(&[1.0][..])
    );
    let restart_candidate_moles = transition.restart_seed[2].exp();
    let restart_phase_total = transition.phase_totals[1];
    assert!(
        (restart_candidate_moles - restart_phase_total).abs()
            <= 1e-12 * restart_phase_total.max(1.0),
        "TPD restart seed {restart_candidate_moles:e} must retain recorded phase total {restart_phase_total:e}"
    );
    assert!(transition.phase_totals[1] > PHASE_CONTROL_TRACE_MOLE_FLOOR);

    for (actual, reference) in moles[..2].iter().zip(&expected.gas_moles) {
        assert!(
            (actual - reference).abs() < 1e-8,
            "canonical gas amount {actual:e} disagrees with independent K_eq {reference:e}"
        );
    }
    let actual_gas_total = moles[0] + moles[1];
    let expected_gas_total = expected.gas_moles.iter().sum::<f64>();
    for (actual, reference) in moles[..2].iter().zip(&expected.gas_moles) {
        assert!(
            (actual / actual_gas_total - reference / expected_gas_total).abs() < 1e-10,
            "canonical and independent gas mole fractions must agree"
        );
    }
    assert!(
        (moles[2] - expected.candidate_moles).abs() < 1e-8,
        "canonical pure-phase amount {} disagrees with independent K_eq {}",
        moles[2],
        expected.candidate_moles
    );
}

#[test]
fn prepared_phase_control_keeps_stable_pure_phase_inactive_without_restart() {
    let independent = PurePhaseBoundaryValidator::default()
        .validate(&boundary_problem_for_target_k(1.0, 1.0))
        .expect("stable independent boundary problem must validate");
    assert_eq!(
        independent.boundary.prediction,
        PurePhaseBoundaryPrediction::StableInactive
    );
    assert!(independent.equilibrium.is_none());

    let outcome = production_phase_control_runner(1.0)
        .solve()
        .expect("a stable inactive pure phase must not make phase control fail");
    let moles = outcome.solution.moles();

    assert_eq!(
        outcome.phase_control_report.initial_phase_set.active_mask(),
        &[true, false]
    );
    assert_eq!(
        outcome.phase_control_report.final_phase_set.active_mask(),
        &[true, false]
    );
    assert!(outcome.phase_control_report.transitions.is_empty());
    assert!((moles[0] - 1.0).abs() < 1e-10);
    assert!((moles[1] - 1.0).abs() < 1e-10);
    let canonical_tpd = outcome.acceptance_report.phase_stability[1]
        .minimum_tpd
        .expect("inactive pure candidate must retain canonical TPD evidence");
    let independent_tpd = independent.boundary.reaction_gibbs_at_absence
        / boundary_problem_for_target_k(1.0, 1.0).candidate_stoichiometry();
    assert!(
        canonical_tpd > 0.0,
        "the stable inactive case must be supported by positive canonical TPD"
    );
    assert!(
        (canonical_tpd - independent_tpd).abs() < 1e-8,
        "stable inactive canonical TPD must equal independent K_eq driving force"
    );
    assert!(
        moles[2] <= PHASE_CONTROL_TRACE_MOLE_FLOOR * 1.000_001,
        "inactive pure phase must remain a trace numerical coordinate, got {} mol",
        moles[2]
    );
}

#[test]
fn prepared_phase_control_deactivates_vanishing_pure_phase_and_restarts_gas_only() {
    // Start S above the test-local destruction threshold so it is a real
    // lifecycle-active phase. K << 1 makes the pure phase thermodynamically
    // unfavorable. A strict log-moles solve has no positive interior root at
    // that pure-phase boundary, so production boundary recovery must solve the
    // reduced gas set, verify its positive TPD, then deactivate and restart.
    let initial_candidate_moles = 1e-2;
    let phase_eps = 1e-3;
    let mut runner =
        production_phase_control_runner_with_candidate_inventory(1e-4, initial_candidate_moles);
    runner.configure_phase_control().phase_eps = phase_eps;

    let outcome = runner
        .solve()
        .expect("production phase control must recover the gas-only boundary state");
    let moles = outcome.solution.moles();
    assert_eq!(
        outcome.phase_control_report.initial_phase_set.active_mask(),
        &[true, true],
        "positive initial S inventory must enter the first fixed-set solve as active"
    );
    assert_eq!(
        outcome.phase_control_report.transitions.len(),
        1,
        "controlled disappearance must take exactly one production deactivation restart"
    );
    assert_eq!(
        outcome.phase_control_report.final_phase_set.active_mask(),
        &[true, false]
    );

    let transition = &outcome.phase_control_report.transitions[0];
    assert!(transition.activated.is_empty());
    assert_eq!(transition.deactivated.len(), 1);
    assert!(
        matches!(
            transition.reason,
            PhaseTransitionReason::BoundaryUnstableActivePhase {
                initial_phase_moles,
                minimum_tpd,
            } if (initial_phase_moles - initial_candidate_moles).abs() < 1e-15 && minimum_tpd > 0.0
        ),
        "pure-phase boundary disappearance must use validated boundary recovery, got {:?}",
        transition.reason
    );
    assert!(
        transition.minimum_tpds[1].is_some_and(|value| value > 0.0),
        "deactivation must retain positive canonical TPD evidence"
    );
    assert!(
        transition.restart_seed[2].exp() <= PHASE_CONTROL_TRACE_MOLE_FLOOR * 1.000_001,
        "restart seed must demote S to the named trace coordinate"
    );
    assert!(
        moles[2] <= PHASE_CONTROL_TRACE_MOLE_FLOOR * 1.000_001,
        "final inactive S amount must remain at trace level, got {} mol",
        moles[2]
    );

    // Removing S does not discard its conserved elements.  For the declared
    // composition A=(1,1), B=(1,0), S=(1,2), gas-only conservation yields
    // A=1+2*nS0 and B=1-nS0 exactly.
    assert!((moles[0] - (1.0 + 2.0 * initial_candidate_moles)).abs() < 1e-10);
    assert!((moles[1] - (1.0 - initial_candidate_moles)).abs() < 1e-10);
    assert!(
        outcome
            .acceptance_report
            .final_validation
            .max_abs_element_balance_error
            < 1e-10,
        "the accepted gas-only restart must retain elemental conservation"
    );
    assert!(outcome.acceptance_report.complementarity.satisfied);
}

#[test]
fn second_synthetic_family_matches_independent_tpd_and_full_lifecycle() {
    let problem = second_family_boundary_problem(9.0);
    let reaction_space = problem
        .validate_strict_independent_family(PurePhaseBoundaryStructuralTolerances::default())
        .expect("second family must have exactly one phase-forming direction");
    assert_eq!(reaction_space.full_reaction_dimension, 1);
    assert_eq!(reaction_space.gas_only_reaction_dimension, 0);
    assert_eq!(problem.candidate_stoichiometry(), 2.0);

    let independent = PurePhaseBoundaryValidator::default()
        .validate(&problem)
        .expect("second independent family must solve");
    assert_eq!(
        independent.boundary.prediction,
        PurePhaseBoundaryPrediction::ShouldAppear
    );
    let expected = independent
        .equilibrium
        .as_ref()
        .expect("K=9 must produce a finite second-family equilibrium");
    let independent_tpd =
        independent.boundary.reaction_gibbs_at_absence / problem.candidate_stoichiometry();
    let canonical_tpd = second_family_canonical_tpd(9.0);
    assert!(
        (canonical_tpd - independent_tpd).abs() < 1e-8,
        "TPD must use per-candidate driving force when nu_S != 1"
    );

    let outcome = second_family_production_runner(9.0)
        .solve()
        .expect("production phase control must activate the second-family pure phase");
    assert_eq!(
        outcome.phase_control_report.transitions.len(),
        1,
        "second-family fixture must converge through exactly one activation transition"
    );
    let transition = outcome
        .phase_control_report
        .transitions
        .first()
        .expect("negative second-family TPD must activate the candidate");
    let transition_tpd = match transition.reason {
        PhaseTransitionReason::UnstableInactivePhase { minimum_tpd } => minimum_tpd,
        _ => panic!("second-family appearance must be driven by inactive-phase TPD"),
    };
    assert!((transition_tpd - independent_tpd).abs() < 1e-8);
    assert_eq!(
        outcome.phase_control_report.final_phase_set.active_mask(),
        &[true, true]
    );

    let moles = outcome.solution.moles();
    for (actual, reference) in moles[..3].iter().zip(&expected.gas_moles) {
        assert!(
            (actual - reference).abs() < 1e-8,
            "second-family canonical gas amount {actual:e} disagrees with independent {reference:e}"
        );
    }
    assert!(
        (moles[3] - expected.candidate_moles).abs() < 1e-8,
        "second-family canonical candidate amount disagrees with independent K_eq"
    );
}

#[test]
fn second_synthetic_family_keeps_positive_tpd_candidate_inactive() {
    let problem = second_family_boundary_problem(1.0);
    let independent = PurePhaseBoundaryValidator::default()
        .validate(&problem)
        .expect("stable second-family independent boundary must validate");
    assert_eq!(
        independent.boundary.prediction,
        PurePhaseBoundaryPrediction::StableInactive
    );
    assert!(independent.equilibrium.is_none());
    let independent_tpd =
        independent.boundary.reaction_gibbs_at_absence / problem.candidate_stoichiometry();

    let outcome = second_family_production_runner(1.0)
        .solve()
        .expect("stable second-family candidate must not make phase control fail");
    assert!(outcome.phase_control_report.transitions.is_empty());
    assert_eq!(
        outcome.phase_control_report.final_phase_set.active_mask(),
        &[true, false]
    );
    let canonical_tpd = outcome.acceptance_report.phase_stability[1]
        .minimum_tpd
        .expect("inactive second-family candidate must retain TPD evidence");
    assert!(canonical_tpd > 0.0);
    assert!((canonical_tpd - independent_tpd).abs() < 1e-8);
}
