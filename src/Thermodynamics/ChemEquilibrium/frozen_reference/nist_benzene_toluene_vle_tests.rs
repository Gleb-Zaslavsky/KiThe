//! I1/I3/I4 contract tests for the frozen NIST benzene/toluene P-x story.
//!
//! The external NIST table is deliberately P-x characterization only. Strict
//! numerical equivalence below is between an independently evaluated ideal
//! Raoult relation from the *local resolved standard states* and canonical
//! TPD minimization using those same states. This keeps an imperfect agreement
//! between NASA records and an experiment visible instead of hiding it behind
//! a tolerance on a different physical source.

use crate::Thermodynamics::ChemEquilibrium::equilibrium_multiphase_domain::{
    MultiphaseEquilibriumLayout, MultiphaseInitialComposition,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_phase_stability::CanonicalPhaseState;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_problem::{
    EquilibriumConditions, TraceSpeciesSeedPolicy,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_workflows::{
    InitialPhaseSet, PhaseSet, PhaseStabilityStatus, compute_phase_stability_reports,
};
use crate::Thermodynamics::ChemEquilibrium::frozen_reference_nist_benzene_toluene_vle::{
    BENZENE_TOLUENE_GAS_PHASE, BENZENE_TOLUENE_LIQUID_PHASE, BENZENE_TOLUENE_TEMPERATURE_K,
    BenzeneTolueneRaoultReference, load_nist_benzene_toluene_vle_dataset, raoult_reference,
    resolve_offline_benzene_toluene_vle,
};
use crate::Thermodynamics::ChemEquilibrium::phase_equilibrium_problem::{
    PhaseEquilibriumBuildRequest, SupportedPhaseModelPolicy, build_phase_equilibrium_problem,
};
use crate::Thermodynamics::ChemEquilibrium::phase_equilibrium_workflow::{
    PhaseControlPolicy, ResolvedPhaseEquilibriumRequest, solve_resolved_pt,
};
use crate::Thermodynamics::ChemEquilibrium::prelude::PhaseIndex;
use crate::Thermodynamics::User_PhaseOrSolution::ResolvedPhaseSystem;
use crate::Thermodynamics::phase_layout::{PhaseComponentId, PhaseId};

const REFERENCE_PRESSURE_PA: f64 = 101_325.0;
const TRACE_MOLES: f64 = 1.0e-30;

#[derive(Debug, Clone, Copy)]
struct LocalRaoultReference {
    bubble_pressure_pa: f64,
    vapor_benzene_mole_fraction: f64,
    vapor_toluene_mole_fraction: f64,
}

/// Computes ideal Raoult equilibrium independently of canonical TPD.
///
/// The relation is obtained directly from the four resolved standard Gibbs
/// values. It is a local I1 oracle, not a fit to the frozen experimental
/// pressures: any NASA-versus-NIST discrepancy remains measurable separately.
fn local_raoult_reference(
    standard_gibbs: &[f64],
    liquid_benzene_mole_fraction: f64,
) -> LocalRaoultReference {
    assert_eq!(standard_gibbs.len(), 4);
    let rt = 8.314_462_618_153_24 * BENZENE_TOLUENE_TEMPERATURE_K;
    let liquid_toluene_mole_fraction = 1.0 - liquid_benzene_mole_fraction;
    let benzene_saturation_pressure =
        REFERENCE_PRESSURE_PA * ((standard_gibbs[2] - standard_gibbs[0]) / rt).exp();
    let toluene_saturation_pressure =
        REFERENCE_PRESSURE_PA * ((standard_gibbs[3] - standard_gibbs[1]) / rt).exp();
    let bubble_pressure_pa = liquid_benzene_mole_fraction * benzene_saturation_pressure
        + liquid_toluene_mole_fraction * toluene_saturation_pressure;
    let vapor_benzene_mole_fraction =
        liquid_benzene_mole_fraction * benzene_saturation_pressure / bubble_pressure_pa;
    let vapor_toluene_mole_fraction =
        liquid_toluene_mole_fraction * toluene_saturation_pressure / bubble_pressure_pa;
    assert!(bubble_pressure_pa.is_finite() && bubble_pressure_pa > 0.0);
    assert!((vapor_benzene_mole_fraction + vapor_toluene_mole_fraction - 1.0).abs() < 1.0e-12);
    LocalRaoultReference {
        bubble_pressure_pa,
        vapor_benzene_mole_fraction,
        vapor_toluene_mole_fraction,
    }
}

/// Extracts local standard Gibbs values through the production bridge before
/// applying the independent Raoult relation. Its pressure is immaterial because
/// standard-state Gibbs energies contain no composition term.
fn local_raoult_from_resolved(
    resolved: &ResolvedPhaseSystem,
    liquid_benzene_mole_fraction: f64,
) -> LocalRaoultReference {
    let layout = MultiphaseEquilibriumLayout::new(resolved.phase_specs().to_vec()).unwrap();
    let probe_composition =
        MultiphaseInitialComposition::from_dense(&layout, vec![0.5, 0.5, TRACE_MOLES, TRACE_MOLES])
            .unwrap();
    let probe = build_phase_equilibrium_problem(
        PhaseEquilibriumBuildRequest::new(
            resolved,
            EquilibriumConditions::new(
                BENZENE_TOLUENE_TEMPERATURE_K,
                REFERENCE_PRESSURE_PA,
                REFERENCE_PRESSURE_PA,
            )
            .unwrap(),
            probe_composition,
            TraceSpeciesSeedPolicy::Absolute { floor: TRACE_MOLES },
            SupportedPhaseModelPolicy::IdealPhaseModelsV1,
        )
        .unwrap(),
    )
    .unwrap();
    let standard_gibbs = probe
        .problem()
        .gibbs()
        .iter()
        .map(|gibbs| gibbs(BENZENE_TOLUENE_TEMPERATURE_K))
        .collect::<Vec<_>>();
    local_raoult_reference(&standard_gibbs, liquid_benzene_mole_fraction)
}

#[test]
fn offline_benzene_toluene_fixture_preserves_state_specific_nasa_provenance() {
    let resolved = resolve_offline_benzene_toluene_vle().unwrap();
    assert!(!resolved.report().nist_fallback_enabled());

    let gas = resolved
        .report()
        .phases()
        .iter()
        .find(|phase| phase.phase().as_option().as_deref() == Some(BENZENE_TOLUENE_GAS_PHASE))
        .unwrap();
    let liquid = resolved
        .report()
        .phases()
        .iter()
        .find(|phase| phase.phase().as_option().as_deref() == Some(BENZENE_TOLUENE_LIQUID_PHASE))
        .unwrap();
    let gas_thermo = gas
        .search()
        .rows()
        .iter()
        .filter(|row| row.property() == "Thermo")
        .collect::<Vec<_>>();
    let liquid_thermo = liquid
        .search()
        .rows()
        .iter()
        .filter(|row| row.property() == "Thermo")
        .collect::<Vec<_>>();
    assert_eq!(gas_thermo.len(), 2);
    assert_eq!(liquid_thermo.len(), 2);
    assert!(gas_thermo.iter().all(|row| row.library() == "NASA_gas"));
    assert!(liquid_thermo.iter().all(|row| row.library() == "NASA_cond"));
}

#[test]
fn canonical_tpd_recovers_local_raoult_binary_liquid_minimum() {
    let dataset = load_nist_benzene_toluene_vle_dataset().unwrap();
    let resolved = resolve_offline_benzene_toluene_vle().unwrap();
    let layout = MultiphaseEquilibriumLayout::new(resolved.phase_specs().to_vec()).unwrap();
    let phase_set = PhaseSet::from_policy(
        &InitialPhaseSet::Explicit {
            active: vec![PhaseIndex::new(0, 2).unwrap()],
            excluded: Vec::new(),
        },
        &[true, false],
    )
    .unwrap();

    // Three interior frozen compositions cover asymmetric and equimolar VLE.
    for frozen in dataset.rows().iter().filter(|row| {
        row.liquid_benzene_mole_fraction > 0.0 && row.liquid_benzene_mole_fraction < 1.0
    }) {
        let external: BenzeneTolueneRaoultReference = raoult_reference(&dataset, frozen).unwrap();
        let local = local_raoult_from_resolved(&resolved, frozen.liquid_benzene_mole_fraction);
        let composition = MultiphaseInitialComposition::from_dense(
            &layout,
            vec![
                local.vapor_benzene_mole_fraction,
                local.vapor_toluene_mole_fraction,
                TRACE_MOLES,
                TRACE_MOLES,
            ],
        )
        .unwrap();
        let bundle = build_phase_equilibrium_problem(
            PhaseEquilibriumBuildRequest::new(
                &resolved,
                EquilibriumConditions::new(
                    BENZENE_TOLUENE_TEMPERATURE_K,
                    local.bubble_pressure_pa,
                    REFERENCE_PRESSURE_PA,
                )
                .unwrap(),
                composition,
                TraceSpeciesSeedPolicy::Absolute { floor: TRACE_MOLES },
                SupportedPhaseModelPolicy::IdealPhaseModelsV1,
            )
            .unwrap(),
        )
        .unwrap();
        let reports = compute_phase_stability_reports(
            bundle.problem().initial_log_moles().as_slice(),
            bundle.problem().gibbs(),
            bundle.problem().phases(),
            &[0, 0, 1, 1],
            bundle.problem().element_composition(),
            BENZENE_TOLUENE_TEMPERATURE_K,
            local.bubble_pressure_pa,
            REFERENCE_PRESSURE_PA,
            &phase_set,
        )
        .unwrap();
        let liquid = reports
            .iter()
            .find(|report| report.phase.index() == 1)
            .unwrap();
        assert_eq!(liquid.status, PhaseStabilityStatus::Evaluated);
        assert!(!liquid.active);
        assert!(liquid.minimum_tpd.unwrap().abs() <= 1.0e-6);
        let incipient = liquid.incipient_composition.as_ref().unwrap();
        assert_eq!(incipient.len(), 2);
        assert!((incipient[0] - frozen.liquid_benzene_mole_fraction).abs() <= 1.0e-8);
        assert!((incipient[1] - (1.0 - frozen.liquid_benzene_mole_fraction)).abs() <= 1.0e-8);
        let minimizer = liquid.minimizer.as_ref().unwrap();
        assert_eq!(minimizer.active_component_count, 2);
        assert!(minimizer.max_abs_constraint_residual <= minimizer.constraint_tolerance);

        // I5 is a characterization target only. It must remain visibly
        // separate from the exact local I1/I3 equivalence asserted above.
        assert!(external.bubble_pressure_pa.is_finite());
    }
}

#[test]
fn bounded_phase_control_activates_the_binary_liquid_only_above_local_bubble_pressure() {
    let resolved = resolve_offline_benzene_toluene_vle().unwrap();
    let layout = MultiphaseEquilibriumLayout::new(resolved.phase_specs().to_vec()).unwrap();
    let liquid_benzene_mole_fraction = 0.5;
    let local = local_raoult_from_resolved(&resolved, liquid_benzene_mole_fraction);
    let initial = MultiphaseInitialComposition::from_dense(
        &layout,
        vec![
            local.vapor_benzene_mole_fraction,
            local.vapor_toluene_mole_fraction,
            0.0,
            0.0,
        ],
    )
    .unwrap();
    let liquid_benzene = PhaseComponentId::new(
        PhaseId::new(Some(BENZENE_TOLUENE_LIQUID_PHASE.to_owned())),
        "C6H6(L)",
    );

    // At fixed overall composition, sub-bubble pressure is vapour-stable and
    // super-bubble pressure must create liquid. This is the physical sign
    // check that prevents an outer-loop activation criterion from being
    // accidentally inverted.
    for (pressure_factor, expect_activation) in [(0.98, false), (1.02, true)] {
        let conditions = EquilibriumConditions::new(
            BENZENE_TOLUENE_TEMPERATURE_K,
            local.bubble_pressure_pa * pressure_factor,
            REFERENCE_PRESSURE_PA,
        )
        .unwrap();
        let solution = solve_resolved_pt(
            ResolvedPhaseEquilibriumRequest::new(&resolved, conditions, initial.clone())
                .with_trace_seed_policy(TraceSpeciesSeedPolicy::Absolute { floor: TRACE_MOLES })
                .with_phase_control_policy(PhaseControlPolicy::default()),
        )
        .unwrap();
        let report = solution.phase_control_report().unwrap();
        assert_eq!(
            report.final_phase_set.active_mask()[1],
            expect_activation,
            "pressure factor {pressure_factor} must decide liquid activity through canonical TPD"
        );
        assert_eq!(solution.phase_control_transitions() > 0, expect_activation);
        let liquid_moles = solution.moles_for(&liquid_benzene).unwrap();
        if expect_activation {
            assert!(liquid_moles > TRACE_MOLES);
            assert!(report.transitions.iter().any(|transition| {
                transition
                    .incipient_composition
                    .as_ref()
                    .is_some_and(|composition| composition.len() == 2)
            }));

            // Validate fully-qualified phase identities rather than collapsing
            // C6H6(g) and C6H6(l) into one bare-substance total. The
            // production residual has separately enforced the two transfer
            // equilibria; recomputing mu here makes that result inspectable.
            let bundle = build_phase_equilibrium_problem(
                PhaseEquilibriumBuildRequest::new(
                    &resolved,
                    conditions,
                    initial.clone(),
                    TraceSpeciesSeedPolicy::Absolute { floor: TRACE_MOLES },
                    SupportedPhaseModelPolicy::IdealPhaseModelsV1,
                )
                .unwrap(),
            )
            .unwrap();
            let state = CanonicalPhaseState::from_log_moles(
                &solution
                    .numerical_component_moles()
                    .iter()
                    .map(|moles| moles.ln())
                    .collect::<Vec<_>>(),
                bundle.problem().gibbs(),
                bundle.problem().phases(),
                &[0, 0, 1, 1],
                bundle.problem().element_composition(),
                BENZENE_TOLUENE_TEMPERATURE_K,
                conditions.pressure(),
                REFERENCE_PRESSURE_PA,
                &[true, true],
                &[true, true],
            )
            .unwrap();
            let chemical_potentials = state.chemical_potentials();
            assert!((chemical_potentials[0] - chemical_potentials[2]).abs() <= 1.0e-4);
            assert!((chemical_potentials[1] - chemical_potentials[3]).abs() <= 1.0e-4);

            // This fixture has no chemical reactions: each named molecular
            // species may move between phases but its gas-plus-liquid amount
            // must remain unchanged.
            assert!(
                ((solution.component_moles()[0] + solution.component_moles()[2])
                    - initial.moles()[0])
                    .abs()
                    <= 1.0e-8
            );
            assert!(
                ((solution.component_moles()[1] + solution.component_moles()[3])
                    - initial.moles()[1])
                    .abs()
                    <= 1.0e-8
            );

            // Element balances are the phase-qualified conservation contract.
            // They cover both C and H despite the duplicate molecular names
            // living in different physical phases.
            for element in 0..bundle.problem().element_composition().ncols() {
                let initial_total = initial
                    .moles()
                    .iter()
                    .enumerate()
                    .map(|(component, &moles)| {
                        bundle.problem().element_composition()[(component, element)] * moles
                    })
                    .sum::<f64>();
                let final_total = solution
                    .component_moles()
                    .iter()
                    .enumerate()
                    .map(|(component, &moles)| {
                        bundle.problem().element_composition()[(component, element)] * moles
                    })
                    .sum::<f64>();
                assert!((final_total - initial_total).abs() <= 1.0e-8);
            }
        } else {
            assert_eq!(liquid_moles, 0.0);
        }
    }
}

/// Release-oriented source characterization, intentionally not a strict
/// NASA-versus-NIST regression. The frozen publication provides P-x, while the
/// local standard states define a separate ideal-mixture prediction.
#[test]
#[ignore = "release-oriented NIST P-x versus local NASA ideal-Raoult characterization"]
fn i5_nist_benzene_toluene_px_characterization() {
    let dataset = load_nist_benzene_toluene_vle_dataset().unwrap();
    let resolved = resolve_offline_benzene_toluene_vle().unwrap();
    println!("NIST ThermoML benzene/toluene P-x characterization at 353.15 K");
    println!("source=P-x only; no experimental vapor composition is inferred");
    println!(
        "{:>7} {:>14} {:>14} {:>14} {:>12}",
        "x_B", "NIST P kPa", "frozen Raoult", "local NASA", "NASA-NIST %"
    );
    for frozen in dataset.rows() {
        let frozen_raoult = raoult_reference(&dataset, frozen).unwrap();
        let local = local_raoult_from_resolved(&resolved, frozen.liquid_benzene_mole_fraction);
        let percent_error = 100.0 * (local.bubble_pressure_pa - frozen.experimental_pressure_pa)
            / frozen.experimental_pressure_pa;
        println!(
            "{:7.3} {:14.4} {:14.4} {:14.4} {:+12.4}",
            frozen.liquid_benzene_mole_fraction,
            frozen.experimental_pressure_pa / 1_000.0,
            frozen_raoult.bubble_pressure_pa / 1_000.0,
            local.bubble_pressure_pa / 1_000.0,
            percent_error,
        );
        assert!(local.bubble_pressure_pa.is_finite() && local.bubble_pressure_pa > 0.0);
    }
}
