//! Full fixed-inventory production lifecycle for the frozen ternary VLE gauge.
//!
//! The flash and bubble/dew calculations in this module are independent
//! Raoult-oracle evidence. Every accepted thermodynamic state is produced by
//! the ordinary immutable phase-control runner, never by a fixture-specific
//! equilibrium routine.

use nalgebra::DMatrix;

use crate::Thermodynamics::ChemEquilibrium::equilibrium_activity::PhaseActivityModel;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_log_moles::{Phase, R, Solvers};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_problem::{
    EquilibriumConditions, EquilibriumProblem, LogMolesInitialGuess,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_solver_policy::SolverPolicy;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_workflows::PhaseManager;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_workflows::PhaseTransitionReason;
use crate::Thermodynamics::ChemEquilibrium::frozen_reference_nist_ternary_vle_antoine_gauge::{
    FROZEN_ANTOINE_REFERENCE_PRESSURE_PA, FrozenAntoineGaugeState, TernaryRaoultFlashPhaseClass,
    gauge_standard_gibbs_j_mol, load_nist_ternary_vle_antoine_gauge, molecular_element_matrix,
    solve_raoult_bubble_temperature, solve_raoult_dew_temperature, solve_raoult_flash,
    ternary_gauge_gibbs_functions,
};
use crate::Thermodynamics::ChemEquilibrium::prepared_phase_control_runner::{
    PreparedPhaseControlOutcome, PreparedPhaseControlRunner,
};

const PRESSURE_PA: f64 = 53_330.0;
const BULK: [f64; 3] = [0.334, 0.333, 0.333];
const TRACE_MOLES: f64 = 1.0e-30;
const COMPOSITION_TOLERANCE: f64 = 3.0e-6;
const CONSERVATION_TOLERANCE: f64 = 2.0e-8;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum ExpectedTopology {
    LiquidOnly,
    TwoPhase,
    GasOnly,
}

impl ExpectedTopology {
    fn phase_mask(self) -> Vec<bool> {
        match self {
            Self::LiquidOnly => vec![false, true],
            Self::TwoPhase => vec![true, true],
            Self::GasOnly => vec![true, false],
        }
    }
}

#[derive(Debug, Clone)]
struct LifecyclePoint {
    temperature_k: f64,
    used_accepted_continuation: bool,
    oracle: TernaryRaoultFlashPhaseClass,
    phase_mask: Vec<bool>,
    beta: f64,
    liquid: Option<[f64; 3]>,
    vapor: Option<[f64; 3]>,
    component_conservation_error: f64,
    chemical_potential_mismatch_j_mol: f64,
    transition_count: usize,
    activation_count: usize,
    deactivation_count: usize,
    gas_minimum_tpd: Option<f64>,
    liquid_minimum_tpd: Option<f64>,
    activation_minimum_tpd: Option<f64>,
    activation_incipient_composition: Option<Vec<f64>>,
}

fn boundaries_and_grid() -> (f64, f64, [f64; 7]) {
    let records = load_nist_ternary_vle_antoine_gauge().unwrap();
    let bubble = solve_raoult_bubble_temperature(records.rows(), BULK, PRESSURE_PA).unwrap();
    let dew = solve_raoult_dew_temperature(records.rows(), BULK, PRESSURE_PA).unwrap();
    assert!((bubble.temperature_k - 375.629_552_965).abs() <= 2.0e-7);
    assert!((dew.temperature_k - 379.074_339_452).abs() <= 2.0e-7);
    assert!((dew.temperature_k - bubble.temperature_k - 3.444_786_486).abs() <= 3.0e-7);
    (
        bubble.temperature_k,
        dew.temperature_k,
        [
            bubble.temperature_k - 3.0,
            bubble.temperature_k - 0.5,
            bubble.temperature_k + 0.5,
            0.5 * (bubble.temperature_k + dew.temperature_k),
            dew.temperature_k - 0.5,
            dew.temperature_k + 0.5,
            dew.temperature_k + 3.0,
        ],
    )
}

fn full_lifecycle_runner(
    temperature_k: f64,
    gas_initially_active: bool,
) -> PreparedPhaseControlRunner {
    let records = load_nist_ternary_vle_antoine_gauge().unwrap();
    let (_, molecular) = molecular_element_matrix(records.rows()).unwrap();
    let element_composition = DMatrix::from_fn(6, molecular.ncols(), |row, column| {
        molecular[(row % 3, column)]
    });
    let initial_moles = if gas_initially_active {
        vec![BULK[0], BULK[1], BULK[2], 0.0, 0.0, 0.0]
    } else {
        vec![0.0, 0.0, 0.0, BULK[0], BULK[1], BULK[2]]
    };
    let seed = initial_moles
        .iter()
        .map(|moles| moles.max(TRACE_MOLES))
        .collect::<Vec<_>>();
    let problem = EquilibriumProblem::new(
        vec![
            "toluene(g)".to_string(),
            "ethylbenzene(g)".to_string(),
            "chlorobenzene(g)".to_string(),
            "toluene(l)".to_string(),
            "ethylbenzene(l)".to_string(),
            "chlorobenzene(l)".to_string(),
        ],
        initial_moles,
        LogMolesInitialGuess::from_moles(&seed, TRACE_MOLES).unwrap(),
        element_composition,
        ternary_gauge_gibbs_functions(records.rows()).unwrap(),
        vec![
            Phase {
                kind: PhaseActivityModel::IdealGas,
                species: vec![0, 1, 2],
            },
            Phase {
                kind: PhaseActivityModel::IdealSolution,
                species: vec![3, 4, 5],
            },
        ],
        EquilibriumConditions::new(
            temperature_k,
            PRESSURE_PA,
            FROZEN_ANTOINE_REFERENCE_PRESSURE_PA,
        )
        .unwrap(),
    )
    .unwrap();
    let mut runner = PreparedPhaseControlRunner::new(problem, Vec::new(), false).unwrap();
    let settings = runner.configure_solver();
    settings.solver_policy = Some(SolverPolicy::legacy_default(Solvers::LM));
    settings.solver_params.tol = 1.0e-11;
    settings.solver_params.max_iter = 250;
    runner
}

fn continue_at(
    runner: &mut PreparedPhaseControlRunner,
    previous: &PreparedPhaseControlOutcome,
    temperature_k: f64,
) {
    let records = load_nist_ternary_vle_antoine_gauge().unwrap();
    runner
        .retarget_numeric(
            EquilibriumConditions::new(
                temperature_k,
                PRESSURE_PA,
                FROZEN_ANTOINE_REFERENCE_PRESSURE_PA,
            )
            .unwrap(),
            LogMolesInitialGuess::new(previous.solution.log_moles().to_vec()).unwrap(),
            ternary_gauge_gibbs_functions(records.rows()).unwrap(),
        )
        .unwrap();
    runner
        .set_continuation_phase_set(previous.phase_control_report.final_phase_set.clone())
        .unwrap();
}

fn accepted_point(
    outcome: &PreparedPhaseControlOutcome,
    temperature_k: f64,
    used_accepted_continuation: bool,
) -> LifecyclePoint {
    let records = load_nist_ternary_vle_antoine_gauge().unwrap();
    let flash = solve_raoult_flash(records.rows(), temperature_k, PRESSURE_PA, BULK).unwrap();
    let moles = outcome.solution.moles();
    let gas_total = moles[..3].iter().sum::<f64>();
    let liquid_total = moles[3..].iter().sum::<f64>();
    let total = gas_total + liquid_total;
    let beta = gas_total / total;
    let liquid = (liquid_total > 1.0e-15)
        .then(|| std::array::from_fn(|index| moles[index + 3] / liquid_total));
    let vapor =
        (gas_total > 1.0e-15).then(|| std::array::from_fn(|index| moles[index] / gas_total));
    let component_conservation_error = (0..3)
        .map(|index| (moles[index] + moles[index + 3] - BULK[index]).abs())
        .fold(0.0_f64, f64::max);
    let chemical_potential_mismatch_j_mol = match (liquid, vapor) {
        (Some(liquid), Some(vapor)) => (0..3)
            .map(|index| {
                let gas_mu = R
                    * temperature_k
                    * (vapor[index] * PRESSURE_PA / FROZEN_ANTOINE_REFERENCE_PRESSURE_PA).ln();
                let liquid_mu = gauge_standard_gibbs_j_mol(
                    &records.rows()[index],
                    FrozenAntoineGaugeState::Liquid,
                    temperature_k,
                )
                .unwrap()
                    + R * temperature_k * liquid[index].ln();
                (gas_mu - liquid_mu).abs()
            })
            .fold(0.0_f64, f64::max),
        _ => 0.0,
    };
    let tpds = outcome
        .acceptance_report
        .phase_stability
        .iter()
        .map(|report| report.minimum_tpd)
        .collect::<Vec<_>>();
    let activation = outcome
        .phase_control_report
        .transitions
        .iter()
        .find(|transition| !transition.activated.is_empty());
    let activation_minimum_tpd = activation.and_then(|transition| match transition.reason {
        PhaseTransitionReason::UnstableInactivePhase { minimum_tpd } => Some(minimum_tpd),
        _ => None,
    });
    let activation_incipient_composition =
        activation.and_then(|transition| transition.incipient_composition.clone());
    LifecyclePoint {
        temperature_k,
        used_accepted_continuation,
        oracle: flash.phase_class,
        phase_mask: outcome.phase_control_report.final_phase_set.active_mask(),
        beta,
        liquid,
        vapor,
        component_conservation_error,
        chemical_potential_mismatch_j_mol,
        transition_count: outcome.phase_control_report.transitions.len(),
        activation_count: outcome
            .phase_control_report
            .transitions
            .iter()
            .filter(|transition| !transition.activated.is_empty())
            .count(),
        deactivation_count: outcome
            .phase_control_report
            .transitions
            .iter()
            .filter(|transition| !transition.deactivated.is_empty())
            .count(),
        gas_minimum_tpd: tpds.first().copied().flatten(),
        liquid_minimum_tpd: tpds.get(1).copied().flatten(),
        activation_minimum_tpd,
        activation_incipient_composition,
    }
}

fn expected_topology(oracle: TernaryRaoultFlashPhaseClass) -> ExpectedTopology {
    match oracle {
        TernaryRaoultFlashPhaseClass::AllLiquid => ExpectedTopology::LiquidOnly,
        TernaryRaoultFlashPhaseClass::TwoPhase => ExpectedTopology::TwoPhase,
        TernaryRaoultFlashPhaseClass::AllVapor => ExpectedTopology::GasOnly,
    }
}

fn assert_point_matches_oracle(point: &LifecyclePoint) {
    assert_eq!(
        point.phase_mask,
        expected_topology(point.oracle).phase_mask()
    );
    assert!(point.component_conservation_error <= CONSERVATION_TOLERANCE);
    match point.oracle {
        TernaryRaoultFlashPhaseClass::AllLiquid => {
            assert!(point.beta <= 1.0e-14);
        }
        TernaryRaoultFlashPhaseClass::AllVapor => {
            assert!((point.beta - 1.0).abs() <= 1.0e-14);
            assert!(
                point
                    .liquid_minimum_tpd
                    .expect("inactive liquid must have TPD")
                    > 0.0
            );
        }
        TernaryRaoultFlashPhaseClass::TwoPhase => {
            let records = load_nist_ternary_vle_antoine_gauge().unwrap();
            let reference =
                solve_raoult_flash(records.rows(), point.temperature_k, PRESSURE_PA, BULK).unwrap();
            assert!(
                (point.beta - reference.vapor_fraction.unwrap()).abs() <= COMPOSITION_TOLERANCE
            );
            let liquid = point.liquid.unwrap();
            let vapor = point.vapor.unwrap();
            for index in 0..3 {
                assert!(
                    (liquid[index] - reference.liquid_mole_fractions.unwrap()[index]).abs()
                        <= COMPOSITION_TOLERANCE
                );
                assert!(
                    (vapor[index] - reference.vapor_mole_fractions.unwrap()[index]).abs()
                        <= COMPOSITION_TOLERANCE
                );
            }
            assert!(point.chemical_potential_mismatch_j_mol <= 2.0e-4);
        }
    }
}

fn run_forward_grid() -> Vec<LifecyclePoint> {
    let (_, _, temperatures) = boundaries_and_grid();
    let mut runner = full_lifecycle_runner(temperatures[0], false);
    let mut previous = runner.solve().unwrap();
    let mut points = vec![accepted_point(&previous, temperatures[0], false)];
    for temperature_k in temperatures.into_iter().skip(1) {
        continue_at(&mut runner, &previous, temperature_k);
        previous = runner.solve().unwrap();
        points.push(accepted_point(&previous, temperature_k, true));
    }
    points
}

/// Runs a new independent descending range. Its first gas-only point receives
/// no carried physical state from the forward sweep; only later points reuse
/// an accepted predecessor from this reverse history.
fn run_reverse_grid() -> Vec<LifecyclePoint> {
    let (_, _, temperatures) = boundaries_and_grid();
    let mut reverse_temperatures = temperatures.into_iter().rev();
    let first_temperature_k = reverse_temperatures.next().unwrap();
    let mut runner = full_lifecycle_runner(first_temperature_k, true);
    let mut previous = runner.solve().unwrap();
    let mut points = vec![accepted_point(&previous, first_temperature_k, false)];
    for temperature_k in reverse_temperatures {
        continue_at(&mut runner, &previous, temperature_k);
        previous = runner.solve().unwrap();
        points.push(accepted_point(&previous, temperature_k, true));
    }
    points
}

/// Evaluates the continuous gas-candidate TPD while intentionally suppressing
/// discrete transitions. This is a *probe only*: the actual history tests use
/// the untouched default `PhaseManager` policy. Keeping these concerns apart
/// prevents a bisection search from silently changing production hysteresis.
fn gas_candidate_tpd_probe(temperature_k: f64) -> f64 {
    let mut runner = full_lifecycle_runner(temperature_k, false);
    runner
        .configure_phase_control()
        .set_explicit_hysteresis(-1.0e9, 1.0e9);
    let outcome = runner.solve().unwrap();
    assert_eq!(
        outcome.phase_control_report.final_phase_set.active_mask(),
        ExpectedTopology::LiquidOnly.phase_mask(),
        "the wide-band probe must retain the liquid-only reference assemblage"
    );
    outcome.acceptance_report.phase_stability[0]
        .minimum_tpd
        .expect("inactive gas candidate must have a continuous TPD")
}

/// Locates a thermodynamic state inside the actual production hysteresis band
/// without assuming a temperature offset from the bubble point.
fn in_band_gas_temperature() -> (f64, f64, f64, f64) {
    let (bubble_k, _, _) = boundaries_and_grid();
    let manager = PhaseManager::default();
    let (dg_create, dg_keep) = manager.thresholds_at(bubble_k).unwrap();
    let target_tpd = 0.5 * (dg_create + dg_keep);
    let mut lower = bubble_k - 0.02;
    let mut upper = bubble_k + 0.02;
    let mut lower_value = gas_candidate_tpd_probe(lower) - target_tpd;
    let upper_value = gas_candidate_tpd_probe(upper) - target_tpd;
    assert!(
        lower_value.signum() != upper_value.signum(),
        "the local physical TPD probe must bracket the actual hysteresis band"
    );
    for _ in 0..48 {
        let middle = 0.5 * (lower + upper);
        let middle_value = gas_candidate_tpd_probe(middle) - target_tpd;
        if middle_value.abs() <= 1.0e-10 {
            return (middle, gas_candidate_tpd_probe(middle), dg_create, dg_keep);
        }
        if middle_value.signum() == lower_value.signum() {
            lower = middle;
            lower_value = middle_value;
        } else {
            upper = middle;
        }
    }
    let temperature_k = 0.5 * (lower + upper);
    (
        temperature_k,
        gas_candidate_tpd_probe(temperature_k),
        dg_create,
        dg_keep,
    )
}

#[test]
fn fixed_inventory_independent_bubble_dew_and_flash_classify_the_seven_point_grid() {
    let records = load_nist_ternary_vle_antoine_gauge().unwrap();
    let bubble = solve_raoult_bubble_temperature(records.rows(), BULK, PRESSURE_PA).unwrap();
    let dew = solve_raoult_dew_temperature(records.rows(), BULK, PRESSURE_PA).unwrap();
    for index in 0..3 {
        assert!(
            (bubble.vapor_mole_fractions[index]
                - [0.500_551_40, 0.232_212_53, 0.267_236_07][index])
                .abs()
                <= 2.0e-8
        );
        assert!(
            (dew.liquid_mole_fractions[index] - [0.201_292_07, 0.426_863_93, 0.371_843_99][index])
                .abs()
                <= 2.0e-8
        );
    }
    assert_eq!(bubble.pressure_pa, PRESSURE_PA);
    assert_eq!(dew.pressure_pa, PRESSURE_PA);
    let (_, _, temperatures) = boundaries_and_grid();
    let expected = [
        TernaryRaoultFlashPhaseClass::AllLiquid,
        TernaryRaoultFlashPhaseClass::AllLiquid,
        TernaryRaoultFlashPhaseClass::TwoPhase,
        TernaryRaoultFlashPhaseClass::TwoPhase,
        TernaryRaoultFlashPhaseClass::TwoPhase,
        TernaryRaoultFlashPhaseClass::AllVapor,
        TernaryRaoultFlashPhaseClass::AllVapor,
    ];
    for (temperature_k, expected) in temperatures.into_iter().zip(expected) {
        let flash = solve_raoult_flash(records.rows(), temperature_k, PRESSURE_PA, BULK).unwrap();
        assert_eq!(flash.phase_class, expected, "T={temperature_k}");
    }
}

#[test]
fn fixed_inventory_forward_production_lifecycle_matches_independent_flash() {
    let points = run_forward_grid();
    assert!(!points[0].used_accepted_continuation);
    assert!(
        points
            .iter()
            .skip(1)
            .all(|point| point.used_accepted_continuation)
    );
    for point in &points {
        assert_point_matches_oracle(point);
    }
    assert_eq!(
        points[2].activation_count, 1,
        "gas must appear between T2 and T3"
    );
    assert!(
        points[2]
            .activation_minimum_tpd
            .expect("gas activation must retain a physical TPD")
            < 0.0
    );
    let gas_incipient = points[2]
        .activation_incipient_composition
        .as_ref()
        .expect("gas activation must retain its three-component TPD minimizer");
    assert_eq!(gas_incipient.len(), 3);
    assert!(gas_incipient.iter().all(|fraction| *fraction > 0.0));
    assert!((gas_incipient.iter().sum::<f64>() - 1.0).abs() <= 1.0e-10);
    assert_eq!(
        points[5].deactivation_count, 1,
        "liquid must disappear between T5 and T6"
    );
    assert_eq!(
        points
            .iter()
            .map(|point| point.activation_count)
            .sum::<usize>(),
        1,
        "forward sweep must not chatter through repeated activation"
    );
    assert_eq!(
        points
            .iter()
            .map(|point| point.deactivation_count)
            .sum::<usize>(),
        1,
        "forward sweep must not chatter through repeated deactivation"
    );
    assert_eq!(
        points
            .iter()
            .map(|point| point.transition_count)
            .sum::<usize>(),
        2,
        "only gas appearance and liquid disappearance are accepted transitions"
    );
}

#[test]
fn fixed_inventory_reverse_production_lifecycle_matches_forward_interior_states() {
    let forward = run_forward_grid();
    let reverse = run_reverse_grid();
    assert!(!reverse[0].used_accepted_continuation);
    assert!(
        reverse
            .iter()
            .skip(1)
            .all(|point| point.used_accepted_continuation)
    );
    for point in &reverse {
        assert_point_matches_oracle(point);
    }
    assert_eq!(
        reverse[2].activation_count, 1,
        "liquid must appear between reverse T6 and T5"
    );
    assert!(
        reverse[2]
            .activation_minimum_tpd
            .expect("liquid activation must retain a physical TPD")
            < 0.0
    );
    let liquid_incipient = reverse[2]
        .activation_incipient_composition
        .as_ref()
        .expect("liquid activation must retain its three-component TPD minimizer");
    assert_eq!(liquid_incipient.len(), 3);
    assert!(liquid_incipient.iter().all(|fraction| *fraction > 0.0));
    assert!((liquid_incipient.iter().sum::<f64>() - 1.0).abs() <= 1.0e-10);
    assert_eq!(
        reverse[5].deactivation_count, 1,
        "gas must disappear between reverse T3 and T2"
    );
    assert_eq!(
        reverse
            .iter()
            .map(|point| point.activation_count)
            .sum::<usize>(),
        1,
        "reverse sweep must not chatter through repeated activation"
    );
    assert_eq!(
        reverse
            .iter()
            .map(|point| point.deactivation_count)
            .sum::<usize>(),
        1,
        "reverse sweep must not chatter through repeated deactivation"
    );

    // Forward is T1..T7 while reverse is T7..T1. The three strictly interior
    // points must agree physically even though their lifecycle history differs.
    for forward_index in 2..=4 {
        let reverse_index = 6 - forward_index;
        let forward_point = &forward[forward_index];
        let reverse_point = &reverse[reverse_index];
        assert!((forward_point.temperature_k - reverse_point.temperature_k).abs() <= 1.0e-12);
        assert!((forward_point.beta - reverse_point.beta).abs() <= COMPOSITION_TOLERANCE);
        for component in 0..3 {
            assert!(
                (forward_point.liquid.unwrap()[component]
                    - reverse_point.liquid.unwrap()[component])
                    .abs()
                    <= COMPOSITION_TOLERANCE
            );
            assert!(
                (forward_point.vapor.unwrap()[component] - reverse_point.vapor.unwrap()[component])
                    .abs()
                    <= COMPOSITION_TOLERANCE
            );
        }
        assert!(
            forward_point
                .component_conservation_error
                .max(reverse_point.component_conservation_error)
                <= CONSERVATION_TOLERANCE
        );
    }
}

#[test]
fn fixed_inventory_gas_hysteresis_has_history_dependent_in_band_and_unique_outside_band_topology() {
    let (in_band_k, probe_tpd, dg_create, dg_keep) = in_band_gas_temperature();
    assert!(
        probe_tpd > dg_create && probe_tpd < dg_keep,
        "the scalar criterion must be strictly inside the production hysteresis band"
    );

    // History A starts from the liquid branch at exactly the in-band target.
    // Its inactive gas candidate must remain inactive because it has not met
    // the stricter creation threshold.
    let mut inactive_history = full_lifecycle_runner(in_band_k, false);
    let inactive = inactive_history.solve().unwrap();
    assert_eq!(
        inactive.phase_control_report.final_phase_set.active_mask(),
        ExpectedTopology::LiquidOnly.phase_mask()
    );
    let inactive_tpd = inactive.acceptance_report.phase_stability[0]
        .minimum_tpd
        .expect("inactive in-band gas must retain TPD evidence");
    assert!((inactive_tpd - probe_tpd).abs() <= 1.0e-7);

    // History B first accepts the adjacent two-phase branch, then carries only
    // that accepted state into the identical in-band target. The gas remains
    // active because retention uses the distinct keep threshold.
    let (_, _, grid) = boundaries_and_grid();
    let mut active_history = full_lifecycle_runner(grid[2], false);
    let accepted_two_phase = active_history.solve().unwrap();
    assert_eq!(
        accepted_two_phase
            .phase_control_report
            .final_phase_set
            .active_mask(),
        ExpectedTopology::TwoPhase.phase_mask()
    );
    continue_at(&mut active_history, &accepted_two_phase, in_band_k);
    let retained = active_history.solve().unwrap();
    assert_eq!(
        retained.phase_control_report.final_phase_set.active_mask(),
        ExpectedTopology::TwoPhase.phase_mask(),
        "an accepted active history must be retained inside the keep band"
    );
    assert!(retained.phase_control_report.transitions.is_empty());

    // Outside the band, discrete history no longer changes the permitted
    // topology. Below the band both routes are liquid-only; above it both
    // routes are two-phase after gas creation.
    let below_band_k = in_band_k - 0.02;
    let above_band_k = in_band_k + 0.02;
    assert!(gas_candidate_tpd_probe(below_band_k) > dg_keep);
    assert!(gas_candidate_tpd_probe(above_band_k) < dg_create);
    let below_fresh = full_lifecycle_runner(below_band_k, false).solve().unwrap();
    assert_eq!(
        below_fresh
            .phase_control_report
            .final_phase_set
            .active_mask(),
        ExpectedTopology::LiquidOnly.phase_mask()
    );
    let above_fresh = full_lifecycle_runner(above_band_k, false).solve().unwrap();
    assert_eq!(
        above_fresh
            .phase_control_report
            .final_phase_set
            .active_mask(),
        ExpectedTopology::TwoPhase.phase_mask()
    );
}

/// Prints the full fixed-inventory lifecycle evidence for release review.
/// The strict unit tests own acceptance; this ignored story preserves the
/// numerical record needed to spot future branch or performance drift.
#[test]
#[ignore = "release-oriented full fixed-inventory ternary VLE lifecycle table"]
fn i5_nist_ternary_fixed_inventory_full_lifecycle_characterization() {
    let forward = run_forward_grid();
    let reverse = run_reverse_grid();
    for point in forward.iter().chain(&reverse) {
        assert_point_matches_oracle(point);
    }
    assert_eq!(
        forward
            .iter()
            .map(|point| point.transition_count)
            .sum::<usize>(),
        2,
        "forward story must retain one appearance and one disappearance"
    );
    assert_eq!(
        reverse
            .iter()
            .map(|point| point.transition_count)
            .sum::<usize>(),
        2,
        "reverse story must retain one appearance and one disappearance"
    );
    for forward_point in forward
        .iter()
        .filter(|point| point.oracle == TernaryRaoultFlashPhaseClass::TwoPhase)
    {
        let reverse_point = reverse
            .iter()
            .find(|point| (point.temperature_k - forward_point.temperature_k).abs() <= 1.0e-12)
            .expect("reverse story must contain every forward temperature");
        assert!((forward_point.beta - reverse_point.beta).abs() <= COMPOSITION_TOLERANCE);
    }
    println!("ternary Antoine-gauge fixed-inventory full lifecycle");
    println!(
        "P={PRESSURE_PA:.0} Pa, z=[0.334, 0.333, 0.333] mol, p0={FROZEN_ANTOINE_REFERENCE_PRESSURE_PA:.0} Pa"
    );
    print_lifecycle_table("forward", &forward);
    print_lifecycle_table("reverse", &reverse);
    print_hysteresis_characterization();

    let (max_beta_delta, max_liquid_delta, max_vapor_delta) = forward
        .iter()
        .chain(&reverse)
        .filter(|point| point.oracle == TernaryRaoultFlashPhaseClass::TwoPhase)
        .fold((0.0_f64, 0.0_f64, 0.0_f64), |metrics, point| {
            let records = load_nist_ternary_vle_antoine_gauge().unwrap();
            let reference =
                solve_raoult_flash(records.rows(), point.temperature_k, PRESSURE_PA, BULK).unwrap();
            let beta_delta = (point.beta - reference.vapor_fraction.unwrap()).abs();
            let liquid_delta = (0..3)
                .map(|index| {
                    (point.liquid.unwrap()[index] - reference.liquid_mole_fractions.unwrap()[index])
                        .abs()
                })
                .fold(0.0_f64, f64::max);
            let vapor_delta = (0..3)
                .map(|index| {
                    (point.vapor.unwrap()[index] - reference.vapor_mole_fractions.unwrap()[index])
                        .abs()
                })
                .fold(0.0_f64, f64::max);
            (
                metrics.0.max(beta_delta),
                metrics.1.max(liquid_delta),
                metrics.2.max(vapor_delta),
            )
        });
    let max_conservation = forward
        .iter()
        .chain(&reverse)
        .map(|point| point.component_conservation_error)
        .fold(0.0_f64, f64::max);
    let max_mu_mismatch = forward
        .iter()
        .chain(&reverse)
        .map(|point| point.chemical_potential_mismatch_j_mol)
        .fold(0.0_f64, f64::max);
    println!(
        "summary: max|delta beta|={max_beta_delta:.3e} max|delta x|={max_liquid_delta:.3e} max|delta y|={max_vapor_delta:.3e} max conservation={max_conservation:.3e} max|delta mu|={max_mu_mismatch:.3e} J/mol"
    );
    println!(
        "accepted transitions: forward={} reverse={} accepted chatter=0",
        forward
            .iter()
            .map(|point| point.transition_count)
            .sum::<usize>(),
        reverse
            .iter()
            .map(|point| point.transition_count)
            .sum::<usize>(),
    );
}

fn print_lifecycle_table(direction: &str, points: &[LifecyclePoint]) {
    println!("{direction}:");
    println!(
        "  T K      route       oracle     topology  beta      transitions TPD_gas   TPD_liquid"
    );
    for point in points {
        println!(
            "  {:>7.3}  {:<10} {:<10?} {:<9} {:>8.5} {:>11} {:>9} {:>11}",
            point.temperature_k,
            if point.used_accepted_continuation {
                "continued"
            } else {
                "initial"
            },
            point.oracle,
            phase_mask_label(&point.phase_mask),
            point.beta,
            point.transition_count,
            format_tpd(point.gas_minimum_tpd),
            format_tpd(point.liquid_minimum_tpd),
        );
        if let (Some(liquid), Some(vapor)) = (point.liquid, point.vapor) {
            let records = load_nist_ternary_vle_antoine_gauge().unwrap();
            let reference =
                solve_raoult_flash(records.rows(), point.temperature_k, PRESSURE_PA, BULK).unwrap();
            println!("    component        x_RR       x_KiThe      y_RR       y_KiThe");
            for index in 0..3 {
                println!(
                    "    {:<14} {:>9.6}  {:>9.6}  {:>9.6}  {:>9.6}",
                    records.rows()[index].compound,
                    reference.liquid_mole_fractions.unwrap()[index],
                    liquid[index],
                    reference.vapor_mole_fractions.unwrap()[index],
                    vapor[index],
                );
            }
        }
    }
}

fn print_hysteresis_characterization() {
    let (temperature_k, probe_tpd, dg_create, dg_keep) = in_band_gas_temperature();
    let inactive = full_lifecycle_runner(temperature_k, false).solve().unwrap();
    let (_, _, grid) = boundaries_and_grid();
    let mut active_runner = full_lifecycle_runner(grid[2], false);
    let accepted_two_phase = active_runner.solve().unwrap();
    continue_at(&mut active_runner, &accepted_two_phase, temperature_k);
    let retained = active_runner.solve().unwrap();
    println!(
        "hysteresis: T={temperature_k:.9} K gas-TPD={probe_tpd:.6e} J/mol band=({dg_create:.6e}, {dg_keep:.6e}) inactive={:?} active-history={:?}",
        inactive.phase_control_report.final_phase_set.active_mask(),
        retained.phase_control_report.final_phase_set.active_mask(),
    );
}

fn phase_mask_label(mask: &[bool]) -> &'static str {
    match mask {
        [false, true] => "liquid",
        [true, true] => "gas+liquid",
        [true, false] => "gas",
        _ => "invalid",
    }
}

fn format_tpd(value: Option<f64>) -> String {
    value
        .map(|value| format!("{value:.2e}"))
        .unwrap_or_else(|| "-".to_string())
}
