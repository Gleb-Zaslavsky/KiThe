//! Production-runner lifecycle evidence for the frozen ternary Antoine gauge.
//!
//! The Raoult flash used here is an independent test oracle.  The actual
//! equilibrium state is always produced by the regular bounded phase-control
//! runner; this module must never grow a ternary-specific solver.

use nalgebra::DMatrix;

use crate::Thermodynamics::ChemEquilibrium::equilibrium_activity::PhaseActivityModel;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_log_moles::{Phase, R, Solvers};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_problem::{
    EquilibriumConditions, EquilibriumProblem, LogMolesInitialGuess,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_solver_policy::SolverPolicy;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_workflows::PhaseTransitionReason;
use crate::Thermodynamics::ChemEquilibrium::frozen_reference_nist_ternary_vle_antoine_gauge::{
    FROZEN_ANTOINE_REFERENCE_PRESSURE_PA, FrozenAntoineGaugeState, TernaryRaoultFlashPhaseClass,
    common_temperature_interval, gauge_standard_gibbs_j_mol, load_nist_ternary_vle_antoine_gauge,
    load_nist_ternary_vle_interior_rows, molecular_element_matrix, solve_raoult_bubble_temperature,
    solve_raoult_flash, ternary_gauge_gibbs_functions,
};
use crate::Thermodynamics::ChemEquilibrium::prepared_phase_control_runner::{
    PreparedPhaseControlOutcome, PreparedPhaseControlRunner,
};

const PRIMARY_ROW_INDEX: usize = 3;
const TRACE_MOLES: f64 = 1.0e-30;
const REFERENCE_VAPOR_FRACTION: f64 = 0.5;

/// Builds a central `z = beta*y + (1-beta)*x` inventory from an independent
/// bubble calculation.  The selected temperature is deliberately below the
/// bubble point, where the gas-only topology must discover liquid.
fn gas_to_liquid_case() -> ([f64; 3], f64, f64) {
    let records = load_nist_ternary_vle_antoine_gauge().unwrap();
    let row = load_nist_ternary_vle_interior_rows().unwrap().rows()[PRIMARY_ROW_INDEX];
    let liquid = [
        row.liquid_toluene_mole_fraction,
        row.liquid_ethylbenzene_mole_fraction,
        row.liquid_chlorobenzene_mole_fraction,
    ];
    let boundary = solve_raoult_bubble_temperature(records.rows(), liquid, row.pressure_pa)
        .expect("the central frozen liquid composition must have a bounded bubble point");
    let bulk = std::array::from_fn(|index| {
        REFERENCE_VAPOR_FRACTION * boundary.vapor_mole_fractions[index]
            + (1.0 - REFERENCE_VAPOR_FRACTION) * liquid[index]
    });
    (bulk, boundary.temperature_k - 0.5, row.pressure_pa)
}

/// Reuses the exact same conserved inventory above its independent bubble
/// boundary and locates its independent dew-side gas-only state. A bubble
/// point alone is insufficient for a fixed bulk inventory: there is generally
/// a finite two-phase interval before the final liquid droplet disappears.
fn liquid_to_gas_case() -> ([f64; 3], f64, f64) {
    let records = load_nist_ternary_vle_antoine_gauge().unwrap();
    let row = load_nist_ternary_vle_interior_rows().unwrap().rows()[PRIMARY_ROW_INDEX];
    let liquid = [
        row.liquid_toluene_mole_fraction,
        row.liquid_ethylbenzene_mole_fraction,
        row.liquid_chlorobenzene_mole_fraction,
    ];
    let boundary = solve_raoult_bubble_temperature(records.rows(), liquid, row.pressure_pa)
        .expect("the central frozen liquid composition must have a bounded bubble point");
    let bulk = std::array::from_fn(|index| {
        REFERENCE_VAPOR_FRACTION * boundary.vapor_mole_fractions[index]
            + (1.0 - REFERENCE_VAPOR_FRACTION) * liquid[index]
    });
    let (_, upper_temperature_k) = common_temperature_interval(records.rows()).unwrap();
    let temperature_k = (1..=256)
        .map(|step| boundary.temperature_k + step as f64 * 0.05)
        .take_while(|temperature| *temperature < upper_temperature_k)
        .find(|temperature| {
            solve_raoult_flash(records.rows(), *temperature, row.pressure_pa, bulk)
                .map(|flash| flash.phase_class == TernaryRaoultFlashPhaseClass::AllVapor)
                .unwrap_or(false)
        })
        .expect(
            "central frozen inventory must reach a gas-only state in the common Antoine interval",
        );
    (bulk, temperature_k, row.pressure_pa)
}

/// Creates a raw immutable phase-control runner only because the Antoine gauge
/// is test-only evidence rather than a `ResolvedPhaseSystem` from a production
/// repository. The runner itself is the production bounded lifecycle path.
fn ternary_runner(
    bulk: [f64; 3],
    temperature_k: f64,
    pressure_pa: f64,
    gas_initially_active: bool,
) -> PreparedPhaseControlRunner {
    let records = load_nist_ternary_vle_antoine_gauge().unwrap();
    let (_, molecular) = molecular_element_matrix(records.rows()).unwrap();
    let elements = DMatrix::from_fn(6, molecular.ncols(), |row, column| {
        molecular[(row % 3, column)]
    });
    let initial_moles = if gas_initially_active {
        vec![bulk[0], bulk[1], bulk[2], 0.0, 0.0, 0.0]
    } else {
        vec![0.0, 0.0, 0.0, bulk[0], bulk[1], bulk[2]]
    };
    let seed_moles = initial_moles
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
        LogMolesInitialGuess::from_moles(&seed_moles, TRACE_MOLES).unwrap(),
        elements,
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
            pressure_pa,
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

fn assert_two_phase_solution_matches_independent_flash(
    outcome: &PreparedPhaseControlOutcome,
    bulk: [f64; 3],
    temperature_k: f64,
    pressure_pa: f64,
) {
    let records = load_nist_ternary_vle_antoine_gauge().unwrap();
    let reference = solve_raoult_flash(records.rows(), temperature_k, pressure_pa, bulk).unwrap();
    assert_eq!(
        reference.phase_class,
        TernaryRaoultFlashPhaseClass::TwoPhase
    );
    let expected_beta = reference.vapor_fraction.unwrap();
    let expected_liquid = reference.liquid_mole_fractions.unwrap();
    let expected_vapor = reference.vapor_mole_fractions.unwrap();
    let moles = outcome.solution.moles();
    let gas_total = moles[..3].iter().sum::<f64>();
    let liquid_total = moles[3..].iter().sum::<f64>();
    let total = gas_total + liquid_total;
    assert!(gas_total > 1.0e-8 && liquid_total > 1.0e-8);
    let actual_beta = gas_total / total;
    assert!((actual_beta - expected_beta).abs() <= 2.0e-6);

    for component in 0..3 {
        let vapor = moles[component] / gas_total;
        let liquid = moles[component + 3] / liquid_total;
        assert!((vapor - expected_vapor[component]).abs() <= 2.0e-6);
        assert!((liquid - expected_liquid[component]).abs() <= 2.0e-6);
        assert!((moles[component] + moles[component + 3] - bulk[component]).abs() <= 2.0e-8);

        let gas_mu =
            R * temperature_k * (vapor * pressure_pa / FROZEN_ANTOINE_REFERENCE_PRESSURE_PA).ln();
        let liquid_mu = gauge_standard_gibbs_j_mol(
            &records.rows()[component],
            FrozenAntoineGaugeState::Liquid,
            temperature_k,
        )
        .unwrap()
            + R * temperature_k * liquid.ln();
        assert!(
            (gas_mu - liquid_mu).abs() <= 2.0e-4,
            "component {component}: mu_g={gas_mu}, mu_l={liquid_mu}"
        );
    }
    assert!(outcome.acceptance_report.complementarity.satisfied);
}

#[test]
fn production_gas_to_liquid_activation_reaches_accepted_ternary_flash_state() {
    let (bulk, temperature_k, pressure_pa) = gas_to_liquid_case();
    let mut runner = ternary_runner(bulk, temperature_k, pressure_pa, true);
    let outcome = runner
        .solve()
        .expect("gas-only state must activate unstable liquid");

    assert_eq!(
        outcome.phase_control_report.initial_phase_set.active_mask(),
        vec![true, false]
    );
    assert_eq!(
        outcome.phase_control_report.final_phase_set.active_mask(),
        vec![true, true]
    );
    let activation = outcome
        .phase_control_report
        .transitions
        .iter()
        .find(|transition| !transition.activated.is_empty())
        .expect("the outer loop must retain liquid activation evidence");
    assert!(matches!(
        activation.reason,
        PhaseTransitionReason::UnstableInactivePhase { minimum_tpd } if minimum_tpd < 0.0
    ));
    assert_eq!(activation.activated.len(), 1);
    let incipient = activation
        .incipient_composition
        .as_ref()
        .expect("multicomponent TPD activation must retain its minimizer composition");
    assert_eq!(incipient.len(), 3);
    assert!(incipient.iter().all(|fraction| *fraction > 0.0));
    assert!((incipient.iter().sum::<f64>() - 1.0).abs() <= 1.0e-10);

    assert_two_phase_solution_matches_independent_flash(&outcome, bulk, temperature_k, pressure_pa);
}

#[test]
fn production_liquid_to_gas_activation_reaches_accepted_ternary_flash_state() {
    let (bulk, temperature_k, pressure_pa) = gas_to_liquid_case();
    let mut runner = ternary_runner(bulk, temperature_k, pressure_pa, false);
    let outcome = runner
        .solve()
        .expect("liquid-only state must activate unstable gas");

    assert_eq!(
        outcome.phase_control_report.initial_phase_set.active_mask(),
        vec![false, true]
    );
    assert_eq!(
        outcome.phase_control_report.final_phase_set.active_mask(),
        vec![true, true]
    );
    let activation = outcome
        .phase_control_report
        .transitions
        .iter()
        .find(|transition| !transition.activated.is_empty())
        .expect("the outer loop must retain gas activation evidence");
    assert!(matches!(
        activation.reason,
        PhaseTransitionReason::UnstableInactivePhase { minimum_tpd } if minimum_tpd < 0.0
    ));
    assert_eq!(activation.activated.len(), 1);
    let incipient = activation
        .incipient_composition
        .as_ref()
        .expect("multicomponent TPD activation must retain its minimizer composition");
    assert_eq!(incipient.len(), 3);
    assert!(incipient.iter().all(|fraction| *fraction > 0.0));
    assert!((incipient.iter().sum::<f64>() - 1.0).abs() <= 1.0e-10);

    assert_two_phase_solution_matches_independent_flash(&outcome, bulk, temperature_k, pressure_pa);
}

#[test]
fn production_liquid_to_gas_lifecycle_deactivates_the_exhausted_liquid_phase() {
    let (bulk, temperature_k, pressure_pa) = liquid_to_gas_case();
    let records = load_nist_ternary_vle_antoine_gauge().unwrap();
    assert_eq!(
        solve_raoult_flash(records.rows(), temperature_k, pressure_pa, bulk)
            .unwrap()
            .phase_class,
        TernaryRaoultFlashPhaseClass::AllVapor
    );

    let mut runner = ternary_runner(bulk, temperature_k, pressure_pa, false);
    let outcome = runner
        .solve()
        .expect("liquid-only state above the bubble boundary must reach stable gas");
    assert_eq!(
        outcome.phase_control_report.initial_phase_set.active_mask(),
        vec![false, true]
    );
    assert_eq!(
        outcome.phase_control_report.final_phase_set.active_mask(),
        vec![true, false]
    );
    assert!(
        outcome
            .phase_control_report
            .transitions
            .iter()
            .any(|transition| !transition.activated.is_empty()),
        "gas must be activated from a physically liquid-only start"
    );
    assert!(
        outcome
            .phase_control_report
            .transitions
            .iter()
            .any(|transition| !transition.deactivated.is_empty()),
        "the now-unsupported liquid phase must be removed by the outer loop"
    );
    let moles = outcome.solution.moles();
    assert!(moles[3..].iter().sum::<f64>() <= 1.0e-20);
    for component in 0..3 {
        assert!((moles[component] - bulk[component]).abs() <= 2.0e-8);
    }
    assert!(outcome.acceptance_report.complementarity.satisfied);
}

#[test]
fn production_temperature_continuation_uses_only_the_previously_accepted_state() {
    let (bulk, two_phase_temperature_k, pressure_pa) = gas_to_liquid_case();
    let (same_bulk, gas_only_temperature_k, same_pressure_pa) = liquid_to_gas_case();
    assert_eq!(bulk, same_bulk);
    assert_eq!(pressure_pa, same_pressure_pa);

    let mut runner = ternary_runner(bulk, two_phase_temperature_k, pressure_pa, true);
    let first = runner
        .solve()
        .expect("the independent first point must establish a two-phase state");
    assert_eq!(
        first.phase_control_report.final_phase_set.active_mask(),
        vec![true, true]
    );
    assert_two_phase_solution_matches_independent_flash(
        &first,
        bulk,
        two_phase_temperature_k,
        pressure_pa,
    );

    // The runner may continue only from this accepted output.  No rejected
    // inner iterate is visible or eligible to seed the next temperature.
    let accepted_seed = LogMolesInitialGuess::new(first.solution.log_moles().to_vec()).unwrap();
    let accepted_phase_set = first.phase_control_report.final_phase_set.clone();
    let records = load_nist_ternary_vle_antoine_gauge().unwrap();
    runner
        .retarget_numeric(
            EquilibriumConditions::new(
                gas_only_temperature_k,
                pressure_pa,
                FROZEN_ANTOINE_REFERENCE_PRESSURE_PA,
            )
            .unwrap(),
            accepted_seed,
            ternary_gauge_gibbs_functions(records.rows()).unwrap(),
        )
        .unwrap();
    runner
        .set_continuation_phase_set(accepted_phase_set)
        .unwrap();

    let second = runner
        .solve()
        .expect("accepted two-phase continuation must reduce to its stable gas topology");
    assert_eq!(
        second.phase_control_report.initial_phase_set.active_mask(),
        vec![true, true],
        "the second point must begin from the first accepted phase topology"
    );
    assert_eq!(
        second.phase_control_report.final_phase_set.active_mask(),
        vec![true, false]
    );
    assert!(
        second
            .phase_control_report
            .transitions
            .iter()
            .any(|transition| !transition.deactivated.is_empty()),
        "the accepted continuation must allow a real topology reduction"
    );
    let moles = second.solution.moles();
    for component in 0..3 {
        assert!((moles[component] - bulk[component]).abs() <= 2.0e-8);
    }
    assert!(second.acceptance_report.complementarity.satisfied);
}

/// Release-only observation of the completed first lifecycle slice.  The
/// independent flash remains the phase-regime oracle; the table exposes the
/// production runner's topology and accepted phase fraction without turning
/// this characterization into a second acceptance contract.
#[test]
#[ignore = "release-oriented ternary Antoine-gauge lifecycle characterization"]
fn i5_nist_ternary_antoine_gauge_production_lifecycle_characterization() {
    println!("ternary Antoine-gauge production lifecycle characterization");
    println!(
        "  story                         T K      flash      initial final  beta       TPD_gas    TPD_liquid transitions"
    );

    let (bulk, two_phase_temperature_k, pressure_pa) = gas_to_liquid_case();
    let records = load_nist_ternary_vle_antoine_gauge().unwrap();
    let mut gas_runner = ternary_runner(bulk, two_phase_temperature_k, pressure_pa, true);
    let gas_to_liquid = gas_runner.solve().unwrap();
    assert_eq!(
        gas_to_liquid
            .phase_control_report
            .final_phase_set
            .active_mask(),
        vec![true, true]
    );
    assert_two_phase_solution_matches_independent_flash(
        &gas_to_liquid,
        bulk,
        two_phase_temperature_k,
        pressure_pa,
    );
    print_lifecycle_row(
        "gas -> liquid activation",
        two_phase_temperature_k,
        records.rows(),
        bulk,
        pressure_pa,
        &gas_to_liquid,
    );

    let mut liquid_runner = ternary_runner(bulk, two_phase_temperature_k, pressure_pa, false);
    let liquid_to_gas = liquid_runner.solve().unwrap();
    assert_eq!(
        liquid_to_gas
            .phase_control_report
            .final_phase_set
            .active_mask(),
        vec![true, true]
    );
    assert_two_phase_solution_matches_independent_flash(
        &liquid_to_gas,
        bulk,
        two_phase_temperature_k,
        pressure_pa,
    );
    print_lifecycle_row(
        "liquid -> gas activation",
        two_phase_temperature_k,
        records.rows(),
        bulk,
        pressure_pa,
        &liquid_to_gas,
    );

    let (bulk, gas_only_temperature_k, pressure_pa) = liquid_to_gas_case();
    let mut disappearance_runner = ternary_runner(bulk, gas_only_temperature_k, pressure_pa, false);
    let disappearance = disappearance_runner.solve().unwrap();
    assert_eq!(
        disappearance
            .phase_control_report
            .final_phase_set
            .active_mask(),
        vec![true, false]
    );
    assert!(disappearance.solution.moles()[3..].iter().sum::<f64>() <= 1.0e-20);
    assert!(disappearance.acceptance_report.complementarity.satisfied);
    print_lifecycle_row(
        "liquid -> gas disappearance",
        gas_only_temperature_k,
        records.rows(),
        bulk,
        pressure_pa,
        &disappearance,
    );
}

fn print_lifecycle_row(
    story: &str,
    temperature_k: f64,
    records: &[crate::Thermodynamics::ChemEquilibrium::frozen_reference_nist_ternary_vle_antoine_gauge::FrozenAntoineStandardState],
    bulk: [f64; 3],
    pressure_pa: f64,
    outcome: &PreparedPhaseControlOutcome,
) {
    let flash = solve_raoult_flash(records, temperature_k, pressure_pa, bulk).unwrap();
    let moles = outcome.solution.moles();
    let gas_total = moles[..3].iter().sum::<f64>();
    let liquid_total = moles[3..].iter().sum::<f64>();
    let beta = gas_total / (gas_total + liquid_total);
    let minimum_tpds = outcome
        .acceptance_report
        .phase_stability
        .iter()
        .map(|report| report.minimum_tpd)
        .collect::<Vec<_>>();
    println!(
        "  {story:<29} {temperature_k:>7.3}  {:<10?} {:<7} {:<6} {beta:>8.5}  {:>9}  {:>11} {:>11}",
        flash.phase_class,
        phase_mask_label(&outcome.phase_control_report.initial_phase_set.active_mask()),
        phase_mask_label(&outcome.phase_control_report.final_phase_set.active_mask()),
        tpd_cell(minimum_tpds.first().copied().flatten()),
        tpd_cell(minimum_tpds.get(1).copied().flatten()),
        outcome.phase_control_report.transitions.len(),
    );
}

fn phase_mask_label(mask: &[bool]) -> &'static str {
    match mask {
        [true, false] => "gas",
        [false, true] => "liquid",
        [true, true] => "gas+liq",
        _ => "invalid",
    }
}

fn tpd_cell(value: Option<f64>) -> String {
    value
        .map(|value| format!("{value:.2e}"))
        .unwrap_or_else(|| "-".to_string())
}
